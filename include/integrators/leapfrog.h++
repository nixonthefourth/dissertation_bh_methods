//
// Created by Mykyta Khomiakov on 25/11/2026.
//

#ifndef MAIN_LEAPFROG_H
#define MAIN_LEAPFROG_H

#include <cstddef>

#include "helpers/object.h++"
#include "integrators/force.h++"

namespace integrators {

    /// @brief One second-order Leapfrog step in kick-drift-kick form
    /// @param bodies The system of stellar bodies; on entry ax, ay, az must hold a(r_n)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size (may be negative)
    /// @param force Force routine, see Force
    /// @details v_{n+1/2} = v_n + h/2 a_n,   r_{n+1} = r_n + h v_{n+1/2},   a_{n+1} = a(r_{n+1}),
    ///          v_{n+1} = v_{n+1/2} + h/2 a_{n+1}.
    ///          One force evaluation per step; on exit ax, ay, az hold a(r_{n+1}), ready for the next step.
    ///          The kick-drift-kick Leapfrog and the Velocity Verlet step are the same scheme, so this
    ///          gives exactly the result of vv_step (it is kept separately as the L2 entry of the notebook).
    void leapfrog_step(Object& bodies, double G, double e, double h, const Force& force);

    /// @brief Advances a system by a number of second-order Leapfrog steps
    /// @param bodies The system of stellar bodies (advanced in place)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size
    /// @param steps Number of steps
    /// @param force Force routine, see Force
    /// @details Computes a(r_0) first, so that the stored accelerations need not be valid on entry.
    ///          Total cost: steps + 1 force evaluations.
    void leapfrog_run(Object& bodies, double G, double e, double h, std::size_t steps, const Force& force);
}

#endif //MAIN_LEAPFROG_H
