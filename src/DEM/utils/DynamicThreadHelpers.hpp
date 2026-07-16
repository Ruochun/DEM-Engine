//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

#ifndef DEME_DYNAMIC_THREAD_HELPERS_HPP
#define DEME_DYNAMIC_THREAD_HELPERS_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <DEM/Defines.h>
#include <DEM/VariableTypes.h>

namespace deme {

// Host-side snapshot of dT's old contact arrays, used only for DEBUG-verbosity lost-contact diagnostics.
struct LostContactDebugSnapshot {
    size_t nPatchContacts = 0;
    size_t nPrimitiveContacts = 0;
    std::vector<bodyID_t> idPatchA;
    std::vector<bodyID_t> idPatchB;
    std::vector<contact_t> contactTypePatch;
    std::vector<bodyID_t> contactPatchIsland;
    std::vector<bodyID_t> idPrimitiveA;
    std::vector<bodyID_t> idPrimitiveB;
    std::vector<contact_t> contactTypePrimitive;
    std::vector<contactPairs_t> geomToPatchMap;
    std::vector<float3> primitivePenetrationStorage;
    std::vector<float3> primitiveAreaStorage;
};

// Internal estimator for tracking the drift/cost trade-off:
// J(d) ~= a + b/d + c*(d/dmax) + e*(d/dmax)^2, with forgetting for a non-stationary baseline.
struct DriftRLS {
    static constexpr int N = 4;  // a, b, c, e
    double theta[N] = {0.0, 0.0, 0.0, 0.0};
    double P[N][N] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
    double lambda = 0.999;
    double sigma2_ema = 1e-6;
    bool initialized = false;

    void reset(double p0 = 1e2) {
        for (int i = 0; i < N; i++) {
            theta[i] = 0.0;
            for (int j = 0; j < N; j++) {
                P[i][j] = (i == j) ? p0 : 0.0;
            }
        }
        sigma2_ema = 1e-6;
        initialized = true;
    }

    static double dot(const double* a, const double* b) {
        double s = 0.0;
        for (int i = 0; i < N; i++) {
            s += a[i] * b[i];
        }
        return s;
    }

    double huberWeight(double r) const {
        const double sigma = std::sqrt(std::max(1e-12, sigma2_ema));
        const double t = 2.5 * sigma;
        const double ar = std::abs(r);
        return (ar <= t) ? 1.0 : t / ar;
    }

    void update(unsigned int d, double d0, double y) {
        if (!initialized) {
            reset();
        }
        if (!std::isfinite(y)) {
            return;
        }

        const double dd = static_cast<double>(std::max(1u, d));
        d0 = std::max(1.0, d0);
        const double x = dd / (dd + d0);
        const double phi[N] = {1.0, d0 / dd, x, x * x};

        const double y_hat = dot(theta, phi);
        if (!std::isfinite(y_hat)) {
            reset();
            return;
        }
        const double r = y - y_hat;

        constexpr double beta = 0.05;
        sigma2_ema = (1.0 - beta) * sigma2_ema + beta * (r * r);
        if (!std::isfinite(sigma2_ema) || sigma2_ema < 0.0) {
            sigma2_ema = 1e-6;
        }

        const double s = std::sqrt(huberWeight(r));
        double phiw[N];
        for (int i = 0; i < N; i++) {
            phiw[i] = s * phi[i];
        }
        const double yw = s * y;

        double Pphi[N] = {0.0, 0.0, 0.0, 0.0};
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                Pphi[i] += P[i][j] * phiw[j];
            }
        }

        const double denom = lambda + dot(phiw, Pphi);
        if (!std::isfinite(denom) || denom <= 1e-18) {
            reset();
            return;
        }

        double K[N];
        for (int i = 0; i < N; i++) {
            K[i] = Pphi[i] / denom;
        }

        const double errw = yw - dot(theta, phiw);
        if (!std::isfinite(errw)) {
            reset();
            return;
        }
        for (int i = 0; i < N; i++) {
            theta[i] += K[i] * errw;
        }
        if (theta[1] < 0.0) {
            theta[1] = 0.0;
        }
        if (theta[2] < 0.0) {
            theta[2] = 0.0;
        }
        if (theta[3] < 0.0) {
            theta[3] = 0.0;
        }
        for (int i = 0; i < N; i++) {
            if (!std::isfinite(theta[i])) {
                reset();
                return;
            }
        }

        double newP[N][N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double ssum = P[i][j];
                for (int k = 0; k < N; k++) {
                    ssum -= (K[i] * phiw[k]) * P[k][j];
                }
                newP[i][j] = ssum / lambda;
            }
        }
        constexpr double P_ABS_MAX = 1e12;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (!std::isfinite(newP[i][j]) || std::abs(newP[i][j]) > P_ABS_MAX) {
                    reset();
                    return;
                }
                P[i][j] = newP[i][j];
            }
        }
    }

    double predict(unsigned int d, double d0) const {
        const double dd = static_cast<double>(std::max(1u, d));
        d0 = std::max(1.0, d0);
        const double x = dd / (dd + d0);
        const double phi[N] = {1.0, d0 / dd, x, x * x};
        return theta[0] * phi[0] + theta[1] * phi[1] + theta[2] * phi[2] + theta[3] * phi[3];
    }
};

// dT-side state used by the adaptive future-drift scheduler. It is deliberately kept outside DEMDynamicThread so
// scheduler tuning state does not bloat the main worker-thread declaration.
struct FutureDriftRegulator {
    double last_total_time = 0.0;
    double debug_cum_time = 0.0;
    double last_debug_cum_time = 0.0;

    uint64_t last_step_sample = 0;
    bool has_last_step_sample = false;

    // calibrateParams may be called multiple times; only the first call after a kT update should advance the
    // timing baseline and update the tuner.
    bool receive_pending = false;
    uint64_t pending_recv_stamp = 0;
    double pending_total_time = 0.0;

    static constexpr int COST_WINDOW = 100;
    double cost_window[COST_WINDOW] = {0.0};
    unsigned int drift_window[COST_WINDOW] = {0u};
    int window_size = 0;
    int window_pos = 0;

    unsigned int last_proposed = 0;
    unsigned int last_sent_proposed = 0;
    unsigned int last_sent_true = 0;
    unsigned int last_sent_wait = 0;
    unsigned int last_wait_cmd = 0;

    double lag_ema = 0.0;
    bool lag_ema_initialized = false;

    uint64_t next_send_step = 0;
    unsigned int next_send_wait = 0;
    bool pending_send = false;

    unsigned int last_observed_kinematic_lag_steps = 0;

    DriftRLS drift_rls;
    uint64_t drift_rls_samples = 0;
    double drift_scale_ema = 0.0;
    bool drift_scale_initialized = false;
    double cost_scale_ema = 0.0;
    bool cost_scale_initialized = false;

    void Clear() { *this = FutureDriftRegulator{}; }
};

inline unsigned clamp_drift_u(unsigned v, unsigned maxv) {
    return (v < 1u) ? 1u : (v > maxv ? maxv : v);
}

inline unsigned clamp_wait_i(int v, unsigned maxv) {
    if (v <= 0) {
        return 0u;
    }
    const unsigned u = (unsigned)v;
    return (u > maxv) ? maxv : u;
}

inline unsigned apply_wait_policy_u(unsigned w,
                                    double lag_pred,
                                    double upper_ratio,
                                    double lower_ratio,
                                    unsigned maxv) {
    w = (w > maxv) ? maxv : w;
    double total = lag_pred + (double)w;
    if (total < 1.0) {
        total = 1.0;
    }
    if (upper_ratio <= 0.0) {
        w = 0u;
    } else if (upper_ratio < 1.0) {
        const unsigned uw = clamp_wait_i((int)std::ceil(total * upper_ratio), maxv);
        if (w > uw) {
            w = uw;
        }
    }
    if (lower_ratio > 0.0) {
        const unsigned lw = clamp_wait_i((int)std::floor(total * lower_ratio), maxv);
        if (w <= lw) {
            w = 0u;
        }
    }
    return w;
}

inline void ema_asym(double& ema, bool& init, double x, double a_up, double a_dn, double minv) {
    if (!init) {
        init = true;
        ema = x;
    } else {
        const double a = (x > ema) ? a_up : a_dn;
        ema = (1.0 - a) * ema + a * x;
    }
    if (!std::isfinite(ema) || ema < minv) {
        ema = std::max(minv, x);
    }
}

inline void ring_push(FutureDriftRegulator& r, double cost, unsigned obs) {
    const int W = FutureDriftRegulator::COST_WINDOW;
    const int i = r.window_pos;
    r.cost_window[i] = cost;
    r.drift_window[i] = obs;
    if (r.window_size < W) {
        r.window_size++;
    }
    r.window_pos = (i + 1) % W;
}

inline double drift_ref_quantile(const FutureDriftRegulator& r, double floor_ref) {
    constexpr int WIN = 30;
    const int n = std::min(r.window_size, WIN);
    if (n <= 0) {
        return floor_ref;
    }
    std::array<unsigned, WIN> a;
    int idx = (r.window_pos > 0) ? (r.window_pos - 1) : (FutureDriftRegulator::COST_WINDOW - 1);
    for (int i = 0; i < n; ++i) {
        a[i] = r.drift_window[idx];
        idx = (idx > 0) ? (idx - 1) : (FutureDriftRegulator::COST_WINDOW - 1);
    }
    const int q = n / 5;
    std::nth_element(a.begin(), a.begin() + q, a.begin() + n);
    double qv = (double)a[q];
    if (r.drift_scale_initialized) {
        qv = std::min(qv, r.drift_scale_ema);
    }
    return std::max(floor_ref, qv);
}

inline bool rls_is_bad(const DriftRLS& rls, unsigned obs, double drift_ref, double scale) {
    const double yhat = rls.predict(obs, drift_ref);
    if (!std::isfinite(yhat) || std::abs(yhat) > 1000.0 * scale) {
        return true;
    }
    const double clip = std::max(1e-3, 1000.0 * scale);
    for (int i = 0; i < DriftRLS::N; ++i) {
        const double t = rls.theta[i];
        if (!std::isfinite(t) || std::abs(t) > clip) {
            return true;
        }
    }
    return false;
}

}  // namespace deme

#endif
