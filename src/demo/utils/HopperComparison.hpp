// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause
#ifndef DEME_HOPPER_COMPARISON_HPP
#define DEME_HOPPER_COMPARISON_HPP

#include "DEM/API.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <stdexcept>

namespace hopper {

// Record the selected force geometry independently from optional legacy/Feng comparison diagnostics.
struct Options {
    bool smoke = false, report = false, diagnostics = false, frames = true, fixedCD = false, help = false;
    std::filesystem::path output;
    std::string geometry = "default";

    static Options Parse(int argc, char** argv) {
        Options o;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--smoke-test")
                o.smoke = true;
            else if (arg == "--comparison-report")
                o.report = true;
            else if (arg == "--feng-diagnostics")
                o.diagnostics = o.report = true;
            else if (arg == "--no-frames")
                o.frames = false;
            else if (arg == "--fixed-cd")
                o.fixedCD = true;
            else if (arg == "--help")
                o.help = true;
            else if (arg == "--geometry=default") {
                o.geometry = "default";
            } else if (arg == "--geometry=feng") {
                o.geometry = "feng";
                o.report = true;
            } else if (arg == "--output-dir" && i + 1 < argc)
                o.output = argv[++i];
            else
                throw std::runtime_error("Unknown or incomplete argument: " + arg);
        }
        if (o.output.empty()) {
            o.output = "DemoOutput_HopperSphereMeshedCylinder";
            if (o.smoke)
                o.output += "_smoke";
            if (o.report)
                o.output += o.geometry == "feng" ? "_feng" : (o.diagnostics ? "_feng_diagnostics" : "_default");
        }
        if (o.report && !o.help) {
            for (const auto* name :
                 {"run.csv", "metrics.csv", "coverage.csv", "pair_runs.csv", "initial.csv", "final.csv"}) {
                if (std::filesystem::exists(o.output / name))
                    throw std::runtime_error("Comparison output already exists; select a fresh --output-dir: " +
                                             o.output.string());
            }
        }
        return o;
    }
};

// One mutually exclusive reason per candidate patch. These are diagnostic gates, not certified force eligibility.
inline unsigned int Reason(const deme::MeshMeshFengDiagnostic& row) {
    if (!row.ownersWatertight)
        return 0;
    if (row.ambiguousPairCount > 0)
        return 1;
    if (row.segmentCount == 0)
        return 2;
    if (!row.closurePassed)
        return 3;
    if (!row.hasContactLine)
        return 4;
    return 5;
}

// Collect host-side reports only at existing output boundaries. Pair runs describe consecutive sampled owner-pair
// presence; they deliberately do not use transient patch indices as persistent contact identifiers.
class Report {
    struct Coverage {
        std::array<size_t, 6> reasons{};
        size_t patches = 0, active = 0, activeEligible = 0, comparable = 0;
        double segments = 0, ambiguousPairs = 0, elasticLoad = 0, eligibleLoad = 0;
        double areaRatioSum = 0, normalCosSum = 0, pointDistanceSum = 0;
        size_t forceEligible = 0, usedFeng = 0;
        double usedElasticLoad = 0;
        std::array<size_t, 7> fallback{};
        std::set<std::pair<deme::bodyID_t, deme::bodyID_t>> pairs;
    };
    struct PairRun {
        size_t first, last, samples;
        double start, end;
        unsigned int category;
    };
    Options options;
    std::ofstream metrics, coverage, runs;
    std::vector<unsigned int> families;
    std::vector<float> masses;
    std::vector<float3> moi;
    std::vector<bool> crossed;
    std::map<std::pair<deme::bodyID_t, deme::bodyID_t>, PairRun> live;
    deme::bodyID_t gate;
    size_t sample = 0;
    double dt;
    double cylinderCompliance, flumeCompliance;
    static constexpr double outletZ = -0.04;

    // Fail on report I/O errors rather than leaving apparently successful but truncated comparison artifacts.
    static std::ofstream Open(const std::filesystem::path& path) {
        std::ofstream out;
        out.exceptions(std::ios::failbit | std::ios::badbit);
        out.open(path);
        out << std::setprecision(17);
        return out;
    }

    unsigned int Category(deme::bodyID_t a, deme::bodyID_t b) const {
        const bool ca = families.at(a) == 100, cb = families.at(b) == 100;
        if (ca && cb)
            return 0;
        if ((ca && b == gate) || (cb && a == gate))
            return 2;
        if ((ca && families.at(b) == 10) || (cb && families.at(a) == 10))
            return 1;
        return 3;
    }

    static const char* CategoryName(unsigned int category) {
        return std::array<const char*, 4>{"cylinder_cylinder", "cylinder_hopper", "cylinder_gate", "other"}[category];
    }

    void WriteRun(const std::pair<deme::bodyID_t, deme::bodyID_t>& pair, const PairRun& run, bool censored) {
        runs << pair.first << ',' << pair.second << ',' << CategoryName(run.category) << ',' << run.first << ','
             << run.last << ',' << run.samples << ',' << run.start << ',' << run.end << ',' << run.end - run.start
             << ',' << censored << '\n';
    }

  public:
    // Metadata records the actual force flavor and diagnostic switch independently, plus timing and particle setup.
    Report(deme::DEMSolver& sim,
           const Options& o,
           deme::bodyID_t gateID,
           size_t spheres,
           size_t cylinders,
           unsigned int sampleSteps,
           double settling,
           double discharge,
           const deme::DEMMaterial& cylinderMaterial,
           const deme::DEMMaterial& flumeMaterial)
        : options(o), gate(gateID), dt(sim.GetTimeStepSize()) {
        const double cylinderE = cylinderMaterial.mat_prop.at("E"), cylinderNu = cylinderMaterial.mat_prop.at("nu");
        const double flumeE = flumeMaterial.mat_prop.at("E"), flumeNu = flumeMaterial.mat_prop.at("nu");
        cylinderCompliance = (1 - cylinderNu * cylinderNu) / cylinderE;
        flumeCompliance = (1 - flumeNu * flumeNu) / flumeE;
        const auto count = static_cast<deme::bodyID_t>(sim.GetNumOwners());
        families = sim.GetOwnerFamily(0, count);
        masses = sim.GetOwnerMass(0, count);
        moi = sim.GetOwnerMOI(0, count);
        crossed.resize(count, false);
        auto metadata = Open(o.output / "run.csv");
        metadata << "schema,geometry,force_model,diagnostics,smoke,fixed_cd,fixed_cd_steps,dt,sample_steps,settling_s,"
                    "discharge_s,"
                    "spheres,cylinders,outlet_z,cylinder_E,cylinder_nu,flume_E,flume_nu\n"
                 << "2," << o.geometry << ",frictional_hertzian," << (o.diagnostics || o.geometry == "feng") << ','
                 << o.smoke << ',' << o.fixedCD << ',' << (o.fixedCD ? 10 : 0) << ',' << dt << ',' << sampleSteps << ','
                 << settling << ',' << discharge << ',' << spheres << ',' << cylinders << ',' << outletZ << ','
                 << cylinderE << ',' << cylinderNu << ',' << flumeE << ',' << flumeNu << '\n';
        metrics = Open(o.output / "metrics.csv");
        metrics << "sample,time,phase,species,count,below_plane,ever_below_plane,discharged_mass,remaining_mass,mean_z,"
                   "translational_ke,rotational_ke,axis_x2,axis_y2,axis_z2,max_speed,wall_seconds\n";
        coverage = Open(o.output / "coverage.csv");
        coverage << "sample,force_time,phase,category,patches,nonwatertight,ambiguous,no_segment,open_boundary,"
                    "degenerate,eligible,"
                    "segments,ambiguous_pairs,active_patches,active_eligible,legacy_elastic_load,eligible_elastic_load,"
                    "comparable,area_ratio_sum,normal_cos_sum,point_distance_sum,owner_pairs,new_owner_pairs,"
                    "force_eligible,used_feng,used_elastic_load,fallback_inactive,fallback_solid,fallback_grouping,"
                    "fallback_gate,fallback_boundary,fallback_normal\n";
        runs = Open(o.output / "pair_runs.csv");
        runs << "owner_a,owner_b,category,first_sample,last_sample,observations,first_force_time,last_force_time,"
                "observed_span,right_censored\n";
    }

    // Owner positions and velocities are post-integration. Angular velocities and MOI are both owner-local.
    // Discharge is the first OBSERVED center below the plane, so it is a sampled passage count, not a flux integral.
    void State(deme::DEMSolver& sim,
               double time,
               const std::string& phase,
               double wall,
               const char* snapshot = nullptr) {
        const auto count = static_cast<deme::bodyID_t>(families.size());
        const auto p = sim.GetOwnerPosition(0, count), v = sim.GetOwnerVelocity(0, count),
                   w = sim.GetOwnerAngVel(0, count);
        const auto q = sim.GetOwnerOriQ(0, count);
        std::ofstream state;
        if (snapshot) {
            state = Open(options.output / snapshot);
            state << "owner,species,x,y,z,vx,vy,vz,wx,wy,wz,qx,qy,qz,qw,mass,ix,iy,iz\n";
        }
        for (const auto family : {99u, 100u}) {
            size_t n = 0, below = 0, ever = 0;
            double totalMass = 0, dischargedMass = 0, z = 0, linear = 0, rotational = 0, vmax = 0;
            double3 axis2 = make_double3(0, 0, 0);
            const char* species = family == 99 ? "sphere" : "cylinder";
            for (size_t i = 0; i < count; ++i) {
                if (families[i] != family)
                    continue;
                for (double value : {double(p[i].x), double(p[i].y), double(p[i].z), double(v[i].x), double(v[i].y),
                                     double(v[i].z), double(w[i].x), double(w[i].y), double(w[i].z), double(q[i].x),
                                     double(q[i].y), double(q[i].z), double(q[i].w)}) {
                    if (!std::isfinite(value))
                        throw std::runtime_error("Non-finite hopper comparison state");
                }
                ++n;
                below += p[i].z < outletZ;
                crossed[i] = crossed[i] || p[i].z < outletZ;
                ever += crossed[i];
                totalMass += masses[i];
                dischargedMass += crossed[i] ? masses[i] : 0;
                z += p[i].z;
                vmax = std::max(vmax, double(length(v[i])));
                linear += .5 * masses[i] * dot(v[i], v[i]);
                rotational += .5 * (moi[i].x * double(w[i].x) * w[i].x + moi[i].y * double(w[i].y) * w[i].y +
                                    moi[i].z * double(w[i].z) * w[i].z);
                if (family == 100) {
                    // Rotate the cylinder's local X axis into world space using DEME's (x,y,z,w) quaternion.
                    const auto qi = q[i];
                    const float3 axis = make_float3(1 - 2 * (qi.y * qi.y + qi.z * qi.z),
                                                    2 * (qi.x * qi.y + qi.w * qi.z), 2 * (qi.x * qi.z - qi.w * qi.y));
                    axis2 += make_double3(axis.x * axis.x, axis.y * axis.y, axis.z * axis.z);
                }
                if (snapshot)
                    state << i << ',' << species << ',' << p[i].x << ',' << p[i].y << ',' << p[i].z << ',' << v[i].x
                          << ',' << v[i].y << ',' << v[i].z << ',' << w[i].x << ',' << w[i].y << ',' << w[i].z << ','
                          << q[i].x << ',' << q[i].y << ',' << q[i].z << ',' << q[i].w << ',' << masses[i] << ','
                          << moi[i].x << ',' << moi[i].y << ',' << moi[i].z << '\n';
            }
            metrics << sample << ',' << time << ',' << phase << ',' << species << ',' << n << ',' << below << ','
                    << ever << ',' << dischargedMass << ',' << totalMass - dischargedMass << ',' << (n ? z / n : 0)
                    << ',' << linear << ',' << rotational << ',' << (n ? axis2.x / n : 0) << ','
                    << (n ? axis2.y / n : 0) << ',' << (n ? axis2.z / n : 0) << ',' << vmax << ',' << wall << '\n';
        }
        ++sample;
    }

    // Weight coverage by the elastic normal term of this demo's unchanged Hertz law. Damping, tangential force,
    // and rolling resistance are deliberately excluded; this proxy must not be reported as measured total force.
    void Geometry(deme::DEMSolver& sim, double time, const std::string& phase) {
        if (!options.diagnostics && options.geometry != "feng")
            return;
        std::array<Coverage, 4> totals;
        std::set<std::pair<deme::bodyID_t, deme::bodyID_t>> present;
        for (const auto& r : sim.GetMeshMeshFengDiagnostics()) {
            const unsigned int category = Category(r.ownerA, r.ownerB);
            auto& c = totals[category];
            ++c.patches;
            ++c.reasons[Reason(r)];
            c.forceEligible += r.fengEligible;
            c.usedFeng += r.usedFeng;
            ++c.fallback[static_cast<unsigned int>(r.fallbackReason)];
            c.segments += r.segmentCount;
            c.ambiguousPairs += r.ambiguousPairCount;
            const auto pair = std::minmax(r.ownerA, r.ownerB);
            c.pairs.insert(pair);
            present.insert(pair);
            if (r.legacyArea > 0 && r.legacyPenetration > 0) {
                ++c.active;
                const auto compliance = [&](deme::bodyID_t owner) {
                    return families.at(owner) == 100 ? cylinderCompliance : flumeCompliance;
                };
                const double effectiveE = 1.0 / (compliance(r.ownerA) + compliance(r.ownerB));
                const double load = (4.0 / 3.0) * effectiveE * std::sqrt(r.legacyArea / deme::PI) * r.legacyPenetration;
                c.elasticLoad += load;
                c.usedElasticLoad += r.usedFeng ? load : 0;
                if (r.hasContactLine) {
                    ++c.activeEligible;
                    c.eligibleLoad += load;
                    ++c.comparable;
                    c.areaRatioSum += r.fengArea / r.legacyArea;
                    c.normalCosSum += dot(r.fengNormal, to_double3(r.legacyNormal));
                    c.pointDistanceSum += length(r.fengContactPoint - r.legacyContactPoint);
                }
            }
        }
        for (unsigned int i = 0; i < totals.size(); ++i) {
            const auto& c = totals[i];
            size_t newPairs = 0;
            for (const auto& pair : c.pairs) {
                auto it = live.find(pair);
                if (it == live.end()) {
                    ++newPairs;
                    live.emplace(pair, PairRun{sample - 1, sample - 1, 1, time - dt, time - dt, i});
                } else {
                    it->second.last = sample - 1;
                    ++it->second.samples;
                    it->second.end = time - dt;
                }
            }
            coverage << sample - 1 << ',' << time - dt << ',' << phase << ',' << CategoryName(i) << ',' << c.patches;
            for (auto reason : c.reasons)
                coverage << ',' << reason;
            coverage << ',' << c.segments << ',' << c.ambiguousPairs << ',' << c.active << ',' << c.activeEligible
                     << ',' << c.elasticLoad << ',' << c.eligibleLoad << ',' << c.comparable << ',' << c.areaRatioSum
                     << ',' << c.normalCosSum << ',' << c.pointDistanceSum << ',' << c.pairs.size() << ',' << newPairs
                     << ',' << c.forceEligible << ',' << c.usedFeng << ',' << c.usedElasticLoad;
            for (size_t reason = 1; reason < c.fallback.size(); ++reason)
                coverage << ',' << c.fallback[reason];
            coverage << '\n';
        }
        for (auto it = live.begin(); it != live.end();) {
            if (!present.count(it->first)) {
                WriteRun(it->first, it->second, false);
                it = live.erase(it);
            } else
                ++it;
        }
    }

    // Mark owner-pair runs still present at the final observation as right-censored.
    void Finish() {
        for (const auto& item : live)
            WriteRun(item.first, item.second, true);
        metrics.flush();
        coverage.flush();
        runs.flush();
    }
};
}  // namespace hopper
#endif
