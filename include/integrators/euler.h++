//
// Created by Mykyta Khomiakov on 25/11/2026.
//

#ifndef MAIN_EULER_H
#define MAIN_EULER_H

#include <cstddef>

#include "helpers/object.h++"
#include "integrators/force.h++"

namespace integrators {

    /// @brief One first-order Forward (explicit) Euler step
    /// @param bodies The system of stellar bodies; on entry ax, ay, az must hold a(r_n)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size (may be negative)
    /// @param force Force routine, see Force
    /// @details r_{n+1} = r_n + h v_n,   v_{n+1} = v_n + h a_n,   a_{n+1} = a(r_{n+1}).
    ///          The displacement is updated with the OLD velocity v_n, so it has to be updated before the
    ///          velocity; with the opposite order the method becomes the symplectic Euler method.
    ///          Non-symplectic, global error O(h). One force evaluation per step; on exit ax, ay, az hold
    ///          a(r_{n+1}), so that the stored accelerations always belong to the stored displacements.
    void euler_step(Object& bodies, double G, double e, double h, const Force& force);

    /// @brief Advances a system by a number of Forward Euler steps
    /// @param bodies The system of stellar bodies (advanced in place)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size
    /// @param steps Number of steps
    /// @param force Force routine, see Force
    /// @details Computes a(r_0) first, so that the stored accelerations need not be valid on entry.
    ///          Total cost: steps + 1 force evaluations.
    void euler_run(Object& bodies, double G, double e, double h, std::size_t steps, const Force& force);
}

#endif //MAIN_EULER_H
