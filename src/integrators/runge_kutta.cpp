//
// Created by Mykyta Khomiakov on 25/11/2026.
//

#include "integrators/runge_kutta.h++"

#include <span>
#include <stdexcept>
#include <utility>

namespace integrators {
    namespace {
        /// @brief out = x + k y for every body (a trial state: r_n + h k_r or v_n + h k_v)
        /// @details out must not be a view of x or y.
        void combine(soa::Vec3Span out, soa::ConstVec3Span x, double k, soa::ConstVec3Span y) {
            soa::per_axis([k](std::span<double> o, std::span<const double> xc, std::span<const double> yc) {
                for (std::size_t i = 0; i < o.size(); ++i) {
                    o[i] = xc[i] + k * yc[i];
                }
            }, out, x, y);
        }
    }

    /// @brief Sizes every buffer for the given system and copies its masses into stage
    /// @param bodies The system of stellar bodies
    void RKWorkspace::prepare(const Object& bodies) {
        const std::size_t N = bodies.size();
        if (stage.size() != N) {
            stage.resize(N);
        }
        stage.mass = bodies.mass;   // same length, so the SoA invariant of stage is kept
        k_r.resize(N);
        sum_r.resize(N);
        sum_v.resize(N);
    }

    /// @brief One second-order Runge-Kutta step (explicit midpoint method)
    /// @param bodies The system of stellar bodies; on entry ax, ay, az must hold a(r_n)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size (may be negative)
    /// @param force Force routine, see Force
    /// @param ws Scratch storage, see RKWorkspace
    void rk2_step(Object& bodies, double G, double e, double h, const Force& force, RKWorkspace& ws) {
        ws.prepare(bodies);     // views are taken after this, since it may resize the buffers

        const Object& state = bodies;
        const soa::ConstVec3Span r = state.displacement();     // r_n
        const soa::ConstVec3Span v = state.velocity();         // v_n = k1_r
        const soa::ConstVec3Span a = state.acceleration();     // a_n = k1_v

        // Stage 2: half a step with the stage 1 slopes, and the slopes there (midpoint)
        combine(ws.stage.displacement(), r, 0.5 * h, v);       // r_mid = r_n + h/2 k1_r
        combine(ws.k_r.view(), v, 0.5 * h, a);                 // k2_r  = v_n + h/2 k1_v
        force(G, ws.stage, e);                                 // k2_v  = a(r_mid) into stage ax, ay, az

        // Advance with the midpoint slopes (r first: it only reads k2_r, not v)
        soa::axpy(bodies.displacement(), h, ws.k_r.view());                    // r += h k2_r
        soa::axpy(bodies.velocity(), h, std::as_const(ws.stage).acceleration());   // v += h k2_v

        force(G, bodies, e);    // a_{n+1}, the k1_v of the next step
    }

    /// @brief One RK2 step with its own scratch storage
    void rk2_step(Object& bodies, double G, double e, double h, const Force& force) {
        RKWorkspace ws;
        rk2_step(bodies, G, e, h, force, ws);
    }

    /// @brief Advances a system by a number of RK2 steps
    /// @param bodies The system of stellar bodies (advanced in place)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size
    /// @param steps Number of steps
    /// @param force Force routine, see Force
    void rk2_run(Object& bodies, double G, double e, double h, std::size_t steps, const Force& force) {
        if (!force) {
            throw std::invalid_argument("rk2_run: no force routine given");
        }

        RKWorkspace ws;
        force(G, bodies, e);    // a_0
        for (std::size_t n = 0; n < steps; ++n) {
            rk2_step(bodies, G, e, h, force, ws);
        }
    }

    /// @brief One classical fourth-order Runge-Kutta step
    /// @param bodies The system of stellar bodies; on entry ax, ay, az must hold a(r_n)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size (may be negative)
    /// @param force Force routine, see Force
    /// @param ws Scratch storage, see RKWorkspace
    /// @details The bodies keep r_n and v_n until the very end, so every stage reads the state at the
    ///          start of the step. k_r always holds the latest displacement slope and the stage
    ///          accelerations the latest velocity slope; each is consumed (added to the weighted sums and
    ///          used to build the next trial state) before it is overwritten.
    void rk4_step(Object& bodies, double G, double e, double h, const Force& force, RKWorkspace& ws) {
        ws.prepare(bodies);     // views are taken after this, since it may resize the buffers

        const Object& state = bodies;
        const soa::ConstVec3Span r = state.displacement();     // r_n
        const soa::ConstVec3Span v = state.velocity();         // v_n = k1_r
        const soa::ConstVec3Span a = state.acceleration();     // a_n = k1_v

        const soa::Vec3Span r_trial = ws.stage.displacement();
        const soa::ConstVec3Span k_v = std::as_const(ws.stage).acceleration();     // filled by force(stage)
        const soa::Vec3Span k_r = ws.k_r.view();
        const soa::Vec3Span sum_r = ws.sum_r.view();
        const soa::Vec3Span sum_v = ws.sum_v.view();

        // Stage 1 is free: k1 = (v_n, a_n)
        soa::copy(v, sum_r);
        soa::copy(a, sum_v);

        // Stage 2: half a step with the stage 1 slopes
        combine(r_trial, r, 0.5 * h, v);        // r_n + h/2 k1_r
        combine(k_r, v, 0.5 * h, a);            // k2_r = v_n + h/2 k1_v
        force(G, ws.stage, e);                  // k2_v
        soa::axpy(sum_r, 2.0, k_r);
        soa::axpy(sum_v, 2.0, k_v);

        // Stage 3: half a step again, now with the stage 2 slopes
        combine(r_trial, r, 0.5 * h, k_r);      // r_n + h/2 k2_r (read k2_r before overwriting it)
        combine(k_r, v, 0.5 * h, k_v);          // k3_r = v_n + h/2 k2_v
        force(G, ws.stage, e);                  // k3_v
        soa::axpy(sum_r, 2.0, k_r);
        soa::axpy(sum_v, 2.0, k_v);

        // Stage 4: a FULL step with the stage 3 slopes
        combine(r_trial, r, h, k_r);            // r_n + h k3_r
        combine(k_r, v, h, k_v);                // k4_r = v_n + h k3_v
        force(G, ws.stage, e);                  // k4_v
        soa::axpy(sum_r, 1.0, k_r);
        soa::axpy(sum_v, 1.0, k_v);

        // Advance with the weighted average of the four slopes (weights 1, 2, 2, 1)
        soa::axpy(bodies.displacement(), h / 6.0, sum_r);
        soa::axpy(bodies.velocity(), h / 6.0, sum_v);

        force(G, bodies, e);    // a_{n+1}, the k1_v of the next step
    }

    /// @brief One RK4 step with its own scratch storage
    void rk4_step(Object& bodies, double G, double e, double h, const Force& force) {
        RKWorkspace ws;
        rk4_step(bodies, G, e, h, force, ws);
    }

    /// @brief Advances a system by a number of RK4 steps
    /// @param bodies The system of stellar bodies (advanced in place)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size
    /// @param steps Number of steps
    /// @param force Force routine, see Force
    void rk4_run(Object& bodies, double G, double e, double h, std::size_t steps, const Force& force) {
        if (!force) {
            throw std::invalid_argument("rk4_run: no force routine given");
        }

        RKWorkspace ws;
        force(G, bodies, e);    // a_0
        for (std::size_t n = 0; n < steps; ++n) {
            rk4_step(bodies, G, e, h, force, ws);
        }
    }
}
