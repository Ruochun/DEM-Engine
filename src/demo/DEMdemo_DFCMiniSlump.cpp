// C++ reproduction of pyDEME_DFCMiniSlump.py for investigating repeated-run stalls.
// Units, material parameters, per-step API calls and output synchronization match Python.
#include "DEM/API.h"
#include "DFCMiniSlump/ConeTransform.hpp"
#include "DFCMiniSlump/Helpers.hpp"
#include <iomanip>
#include <numeric>

using namespace deme;
namespace fs = std::filesystem;

int main(int argc, char** argv) {
    try {
        std::cout << std::unitbuf;
        auto args = minislump::parse(argc, argv);
        // Never recursively delete an arbitrary user path. Repetitions use separate directories.
        fs::create_directories(args.output);
        if (!fs::is_empty(args.output))
            throw std::runtime_error("Output directory must be empty: " + args.output.string());
        const fs::path data = DFC_MINISLUMP_DATA_DIR;
        constexpr double pi = 3.14159265358979323846;
        const double specimen_volume =
            pi * std::pow(101.6 / 2, 2) * 50.8 * 2 / 3 - pi * std::pow(69.85 / 2, 2) * 50.8 / 3;
        const double rho = 845 * (1 + .21 + 1.2422) * 1e-12;
        const double va0 = (1 - std::pow(.5, .5)) * (1 - 845.0 / 3150 - .21 * 845 / 1000 - .03);
        DEMSolver sim;
        sim.SetVerbosity("INFO");
        sim.SetOutputFormat("CSV");
        sim.SetOutputContent(std::vector<std::string>{"ABSV", "XYZ"});
        sim.SetMeshOutputFormat("VTK");
        sim.SetContactOutputContent(std::vector<std::string>{"OWNER", "CNT_WILDCARD"});
        if (args.diagnostic) {
            sim.SetContactOutputFormat("CSV");
            sim.SetContactOutputContent(std::vector<std::string>{"OWNER", "FORCE", "POINT", "NORMAL", "CNT_WILDCARD"});
        }
        sim.SetJitifyClumpTemplates(false);
        sim.SetJitifyMassProperties(false);
        sim.InstructBoxDomainDimension(1000, 1000, 1000);
        sim.SetErrorOutAvgContacts(250);
        sim.SetCudaDebugSync();
        sim.SetGravitationalAcceleration(make_float3(0, -9810, 0));
        auto mat = sim.LoadMaterial({{"ENm_mat", .021},
                                     {"ENa_mat", 100},
                                     {"mortar_layer_mat", 2.15},
                                     {"alpha_mat", .25},
                                     {"beta_mat", .5},
                                     {"np_mat", 1},
                                     {"sgmTmax_mat", 3e-4},
                                     {"sgmTau0_mat", 2e-5},
                                     {"kappa0_mat", 100},
                                     {"eta_inf_mat", 1.25e-5},
                                     {"flocbeta_mat", .01},
                                     {"flocTcr_mat", 100},
                                     {"flocm_mat", 1},
                                     {"lambda_init_mat", 1e-8}});
        auto floor = sim.AddBCPlane(make_float3(0, 0, 0), make_float3(0, 1, 0), mat);
        floor->SetFamily(10);
        auto bottom = sim.Track(floor);
        if (!args.no_cone) {
            auto cone = sim.AddWavefrontMeshObject((data / "minislump-Cut4.obj").string(), mat);
            cone->SetMass(1);
            cone->SetInitPos(make_float3(CONE_POSITION[0], CONE_POSITION[1], CONE_POSITION[2]));
            cone->SetInitQuat(
                make_float4(CONE_QUATERNION[0], CONE_QUATERNION[1], CONE_QUATERNION[2], CONE_QUATERNION[3]));
            cone->SetFamily(11);
            sim.Track(cone);
        }
        sim.SetFamilyFixed(10);
        sim.SetFamilyFixed(11);
        auto model = sim.ReadContactForceModel((data / "DFCModel.cu").string());
        model->SetMustHaveMatProp({});
        model->SetMustPairwiseMatProp({});
        model->SetPerContactWildcards({"contact_info_step_time", "contact_info_lambda", "contact_info_strain_x",
                                       "contact_info_strain_y", "contact_info_strain_z"});
        model->SetPerOwnerWildcards({"ParticleThixoAvg"});
        model->DefineCustomModelPrerequisites(minislump::read(data / "prerequisites.cuh"));

        // Normalize coated sphere masses without subtracting overlaps, and rotate -90 degrees about X.
        std::istringstream input(minislump::read(args.particles));
        std::vector<float3> positions;
        std::vector<double> radii, volumes;
        std::string line;
        while (std::getline(input, line)) {
            if (line.rfind("//", 0) == 0)
                continue;
            std::istringstream row(line);
            long id;
            double x, y, z, diameter;
            std::string extra;
            if (!(row >> id >> x >> y >> z >> diameter) || (row >> extra))
                continue;
            double radius = diameter / 2 + 2.15;
            if (!std::isfinite(radius) || radius <= 0 || !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
                throw std::runtime_error("Invalid particle geometry");
            positions.push_back(make_float3(x, z, -y));
            radii.push_back(radius);
            volumes.push_back(4.0 / 3 * pi * std::pow(radius, 3));
        }
        if (radii.empty() || radii.size() > std::numeric_limits<bodyID_t>::max())
            throw std::runtime_error("Invalid particle count");
        double volume_sum = std::accumulate(volumes.begin(), volumes.end(), 0.0);
        double effective_density = rho * specimen_volume / volume_sum;
        std::vector<std::shared_ptr<DEMClumpTemplate>> templates;
        for (size_t i = 0; i < radii.size(); ++i)
            templates.push_back(sim.LoadSphereType(effective_density * volumes[i], radii[i], mat));
        auto particles = sim.Track(sim.AddClumps(templates, positions));
        std::cout << "Loaded " << radii.size() << " particles; effective density=" << effective_density << '\n';
        auto plate = sim.AddBCPlane(make_float3(0, args.plate, 0), make_float3(0, -1, 0), mat);
        plate->SetFamily(11);
        auto total_ke = sim.CreateInspector("clump_kinetic_energy");
        double dt = args.diagnostic ? 1e-10 : (args.sensitivity ? args.sensitivity : 1e-8);
        sim.SetInitTimeStep(dt);
        sim.SetMaxVelocity(400000);
        sim.SetErrorOutVelocity(500000);
        sim.SetInitBinNumTarget(args.bins);
        sim.DisableAdaptiveBinSize();
        std::cout << "Starting Initialize()\n";
        sim.Initialize();
        std::cout << "Finished Initialize()\n";
        const auto count = static_cast<bodyID_t>(radii.size());
        const auto masses = sim.GetOwnerMass(0, count);
        const double mass = std::accumulate(masses.begin(), masses.end(), 0.0), weight = mass * 9810;
        const double bottom_mass = bottom->Mass();
        if (!std::isfinite(bottom_mass) || bottom_mass <= 0 || weight <= 0)
            throw std::runtime_error("Invalid mass/weight");
        auto translational_ke = [&]() {
            auto v = sim.GetOwnerVelocity(0, count);
            double ke = 0;
            for (size_t i = 0; i < v.size(); ++i)
                ke += .5 * masses[i] * (double(v[i].x) * v[i].x + double(v[i].y) * v[i].y + double(v[i].z) * v[i].z);
            return ke;
        };
        // Capture initial forces and optionally the first and refreshed steps after contact suppression.
        if (args.diagnostic) {
            auto owners = [&](const std::string& name) {
                auto n = static_cast<bodyID_t>(sim.GetNumOwners());
                auto p = sim.GetOwnerPosition(0, n), v = sim.GetOwnerVelocity(0, n), a = sim.GetOwnerAcc(0, n);
                auto m = sim.GetOwnerMass(0, n);
                auto f = sim.GetOwnerFamily(0, n);
                auto file = minislump::output(args.output / name);
                file << "owner,family,mass,x,y,z,vx,vy,vz,ax,ay,az\n";
                for (bodyID_t i = 0; i < n; ++i)
                    file << i << ',' << unsigned(f[i]) << ',' << m[i] << ',' << p[i].x << ',' << p[i].y << ',' << p[i].z
                         << ',' << v[i].x << ',' << v[i].y << ',' << v[i].z << ',' << a[i].x << ',' << a[i].y << ','
                         << a[i].z << '\n';
            };
            owners("owners_initial.csv");
            sim.WriteMeshFile((args.output / "cone_initial.vtk").string());
            sim.WriteSphereFile((args.output / "particles_initial.csv").string());
            sim.DoDynamicsThenSync(dt / 2);
            sim.WriteContactFileIncludingPotentialPairs((args.output / "contacts_initial.csv").string());
            owners("owners_after_one_step.csv");
            if (args.release) {
                sim.DoDynamicsThenSync(0);
                sim.DisableContactBetweenFamilies(11, 0);
                sim.DoDynamicsThenSync(dt / 2);
                sim.WriteContactFileIncludingPotentialPairs((args.output / "contacts_release_first_step.csv").string());
                sim.DoDynamicsThenSync(0);
                sim.DoDynamicsThenSync(dt / 2);
                sim.WriteContactFileIncludingPotentialPairs((args.output / "contacts_release_refreshed.csv").string());
            }
            std::cout << "Diagnostic complete; destroying solver\n";
            return 0;
        }
        auto manifest = minislump::output(args.output / "frame_times.csv");
        auto height = minislump::output(args.output / "MiniSlump_height_vs_time.csv");
        auto stats = minislump::output(args.output / "MiniSlump_stats_vs_time.csv");
        auto hk = minislump::output(args.output / "height_KE_vs_time.csv");
        auto reaction = minislump::output(args.output / "bottom_reaction_vs_time.csv");
        // Retain trajectories until the end, as in Python, to preserve output timing.
        std::ostringstream trajectories;
        trajectories.precision(17);
        trajectories << "time,particle_id,x,y,z\n";
        manifest << "frame,step_index,completed_steps,time,dt_used,phase\n";
        height << "time,top_height,top_particle_x,top_particle_y,top_particle_z,average_top_layer_height,top_layer_"
                  "particle_count\n";
        stats << "time,va0,KE,KE_total\n";
        hk << "frame,time,phase,top_height,average_top_layer_height,top_layer_particle_count,KE,KE_total,top_radius,"
              "bottom_radius,bottom_layer_particle_count\n";
        reaction << "frame,time,phase,plate_contact_Fx_N,plate_contact_Fy_N,plate_contact_Fz_N,bottom_reaction_N,total_"
                    "particle_mass,gravity_y_mm_s2,particle_weight_N,R_over_W\n";
        long long frame = 0, step = 0, completed = 0;
        double t = 0, last_dt = 0, out_time = 0;
        auto write_frame = [&](const std::string& phase) {
            sim.DoDynamicsThenSync(0);
            auto acc = bottom->ContactAcc();
            double fx = bottom_mass * acc.x, fy = bottom_mass * acc.y, fz = bottom_mass * acc.z;
            reaction << frame << ',' << t << ',' << phase << ',' << fx << ',' << fy << ',' << fz << ',' << -fy << ','
                     << mass << ",-9810," << weight << ',' << -fy / weight << std::endl;
            std::cout << "Time " << t << ", " << phase << ", frame " << frame << ", last dt " << last_dt << ", next dt "
                      << dt << '\n';
            std::ostringstream suffix;
            suffix << std::setw(6) << std::setfill('0') << frame;
            sim.WriteSphereFile((args.output / ("MiniSlump_va_M" + suffix.str() + ".csv")).string());
            if (!args.no_cone)
                sim.WriteMeshFile((args.output / ("MiniSlump_mesh" + suffix.str() + ".vtk")).string());
            manifest << frame << ',' << step << ',' << completed << ',' << t << ',' << last_dt << ',' << phase
                     << std::endl;
            auto p = particles->Positions();
            for (size_t i = 0; i < std::min(size_t(50), p.size()); ++i)
                trajectories << t << ',' << i << ',' << p[i].x << ',' << p[i].y << ',' << p[i].z << '\n';
            double ke = translational_ke();
            auto s = minislump::surface(p, radii);
            auto top = p[s.top_id];
            height << t << ',' << s.top << ',' << top.x << ',' << top.y << ',' << top.z << ',' << s.average << ','
                   << s.top_count << std::endl;
            double kinetic = total_ke->GetValue();
            stats << t << ',' << va0 << ',' << ke << ',' << kinetic << std::endl;
            hk << frame << ',' << t << ',' << phase << ',' << s.top << ',' << s.average << ',' << s.top_count << ','
               << ke << ',' << kinetic << ',' << s.top_radius << ',' << s.bottom_radius << ',' << s.bottom_count
               << std::endl;
            std::cout << "Contacts " << sim.GetNumContacts() << ", translational KE " << ke << '\n';
        };
        auto normal_frame = [&](const std::string& phase) {
            if (args.frame_time ? t >= out_time : step % args.frame_steps == 0) {
                write_frame(phase);
                ++frame;
                if (args.frame_time)
                    out_time += args.frame_time;
            }
        };
        if (args.sensitivity) {
            auto steps = std::llround(.1 / dt), interval = std::llround(.001 / dt);
            sim.DoDynamicsThenSync(0);
            write_frame("fixed-step settling");
            ++frame;
            for (step = 1; step <= steps; ++step) {
                sim.DoDynamics(dt);
                t = step * dt;
                completed = step;
                last_dt = dt;
                if (step % interval == 0) {
                    sim.DoDynamicsThenSync(0);
                    write_frame("fixed-step settling");
                    ++frame;
                }
            }
        } else {
            if (args.settle)
                while (true) {
                    if (step % 20 == 0 && dt < 5e-5) {
                        dt = std::min(dt * 1.05, 5e-5);
                        sim.UpdateStepSize(dt);
                    }
                    if (dt > 1.5 - t) {
                        dt = 1.5 - t;
                        sim.UpdateStepSize(dt);
                    }
                    last_dt = dt;
                    sim.DoDynamics(dt);
                    t += dt;
                    ++completed;
                    if (t >= 1.5) {
                        double ke = translational_ke();
                        auto release = minislump::output(args.output / "settling_release.json");
                        release << "{\"time\":" << t << ",\"completed_steps\":" << completed
                                << ",\"step_index\":" << step
                                << ",\"reason\":\"fixed_settling_time\",\"translational_KE\":" << ke
                                << ",\"settling_end_time\":1.5,\"sgmTmax\":0.0003,\"lambda_init\":1e-8}\n";
                        break;
                    }
                    normal_frame("settling");
                    ++step;
                    if (step % 1000 == 0)
                        std::cout << "Progress settling: step=" << step << ", time=" << t << '\n';
                }
            double phase_end = t;
            dt = 1e-8;
            sim.DoDynamicsThenSync(0);
            sim.UpdateStepSize(dt);
            sim.DisableContactBetweenFamilies(11, 0);
            while (t < phase_end + args.duration) {
                last_dt = dt;
                sim.DoDynamics(dt);
                t += dt;
                ++completed;
                // Python/Chrono update AFTER integration on nine out of ten step indices.
                if (step % 10 != 0 && dt < 5e-6) {
                    dt = std::min(dt * 1.01, 5e-6);
                    sim.UpdateStepSize(dt);
                }
                normal_frame("slumping");
                ++step;
                if (step % 1000 == 0)
                    std::cout << "Progress slumping: step=" << step << ", time=" << t << '\n';
            }
        }
        auto final = minislump::output(args.output / args.out_name);
        final << "# X,Y,Z\n";
        for (auto p : particles->Positions())
            final << p.x << ',' << p.y << ',' << p.z << '\n';
        auto trajectory_file = minislump::output(args.output / "particle_trajectories.csv");
        trajectory_file << trajectories.str();
        std::cout << "Simulation complete; destroying solver\n";
    } catch (const std::exception& e) {
        std::cerr << "MiniSlump: " << e.what() << '\n';
        return 1;
    }
    std::cout << "Solver destroyed; exiting\n";
    return 0;
}
