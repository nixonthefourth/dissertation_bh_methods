//
// Created by Mykyta Khomiakov on 25/11/2026.
//

#ifndef MAIN_RUNGE_KUTTA_H
#define MAIN_RUNGE_KUTTA_H

#include <cstddef>

#include "helpers/object.h++"
#include "helpers/soa.h++"
#include "integrators/force.h++"

namespace integrators {

    /// @brief Scratch storage of the Runge-Kutta steps (allocate once, reuse every step)
    /// @details Unlike the kick-drift schemes, a Runge-Kutta stage needs the accelerations of a trial
    ///          configuration, which is not the state the bodies are in. The notebook's acc_stage swaps
    ///          the trial displacements into the bodies temporarily; here they are written into a
    ///          separate system, stage, which holds the masses of the real system. The force routine is
    ///          called on stage, so the real system is never touched until the step is complete, and any
    ///          Force (direct summation or Barnes-Hut) works without modification.
    ///              - stage          trial displacements and masses; the force writes the stage slopes
    ///                               k_v into its ax, ay, az (its velocities are not used)
    ///              - k_r            displacement slope of the latest stage
    ///              - sum_r, sum_v   weighted sums of the slopes (RK4 only)
    ///          The step functions call prepare() themselves, so a workspace can be passed in any state.
    struct RKWorkspace {
        Object stage;
        soa::Vec3Array k_r;
        soa::Vec3Array sum_r, sum_v;

        /// @brief Sizes every buffer for the given system and copies its masses into stage
        void prepare(const Object& bodies);
    };

    /// @brief One second-order Runge-Kutta step (explicit midpoint method; Runge 1895)
    /// @param bodies The system of stellar bodies; on entry ax, ay, az must hold a(r_n)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size (may be negative)
    /// @param force Force routine, see Force
    /// @param ws Scratch storage, see RKWorkspace
    /// @details Stage 1 is free: k1_r = v_n, k1_v = a_n. Stage 2 is the midpoint:
    ///              k2_r = v_n + h/2 a_n,   k2_v = a(r_n + h/2 v_n),
    ///          and the step is r_{n+1} = r_n + h k2_r, v_{n+1} = v_n + h k2_v.
    ///          Non-symplectic, global error O(h^2). Two force evaluations per step (the midpoint and
    ///          a(r_{n+1}), which is the k1_v of the next step); on exit ax, ay, az hold a(r_{n+1}).
    void rk2_step(Object& bodies, double G, double e, double h, const Force& force, RKWorkspace& ws);

    /// @brief One RK2 step with its own scratch storage (convenient for a single step; allocates)
    void rk2_step(Object& bodies, double G, double e, double h, const Force& force);

    /// @brief Advances a system by a number of RK2 steps
    /// @param bodies The system of stellar bodies (advanced in place)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size
    /// @param steps Number of steps
    /// @param force Force routine, see Force
    /// @details Computes a(r_0) first and allocates the scratch storage once. Total cost: 2 steps + 1
    ///          force evaluations.
    void rk2_run(Object& bodies, double G, double e, double h, std::size_t steps, const Force& force);

    /// @brief One classical fourth-order Runge-Kutta step (Kutta 1901)
    /// @param bodies The system of stellar bodies; on entry ax, ay, az must hold a(r_n)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size (may be negative)
    /// @param force Force routine, see Force
    /// @param ws Scratch storage, see RKWorkspace
    /// @details k1 = (v_n, a_n),
    ///          k2 = (v_n + h/2 k1_v, a(r_n + h/2 k1_r)),
    ///          k3 = (v_n + h/2 k2_v, a(r_n + h/2 k2_r)),
    ///          k4 = (v_n + h   k3_v, a(r_n + h   k3_r)),
    ///          r_{n+1} = r_n + h/6 (k1_r + 2 k2_r + 2 k3_r + k4_r), and likewise for v.
    ///          Non-symplectic, global error O(h^4). Four force evaluations per step (three stages and
    ///          a(r_{n+1}), which is the k1_v of the next step); on exit ax, ay, az hold a(r_{n+1}).
    void rk4_step(Object& bodies, double G, double e, double h, const Force& force, RKWorkspace& ws);

    /// @brief One RK4 step with its own scratch storage (convenient for a single step; allocates)
    void rk4_step(Object& bodies, double G, double e, double h, const Force& force);

    /// @brief Advances a system by a number of RK4 steps
    /// @param bodies The system of stellar bodies (advanced in place)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size
    /// @param steps Number of steps
    /// @param force Force routine, see Force
    /// @details Computes a(r_0) first and allocates the scratch storage once. Total cost: 4 steps + 1
    ///          force evaluations.
    void rk4_run(Object& bodies, double G, double e, double h, std::size_t steps, const Force& force);
}

#endif //MAIN_RUNGE_KUTTA_H
