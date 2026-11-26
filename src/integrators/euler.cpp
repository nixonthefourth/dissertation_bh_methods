//
// Created by Mykyta Khomiakov on 25/11/2026.
//

#include "integrators/euler.h++"

#include <stdexcept>

namespace integrators {

    /// @brief One first-order Forward (explicit) Euler step
    /// @param bodies The system of stellar bodies; on entry ax, ay, az must hold a(r_n)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size (may be negative)
    /// @param force Force routine, see Force
    void euler_step(Object& bodies, double G, double e, double h, const Force& force) {
        soa::axpy(bodies.displacement(), h, bodies.velocity());         // r += h v_n (old velocity)
        soa::axpy(bodies.velocity(), h, bodies.acceleration());         // v += h a_n (old acceleration)
        force(G, bodies, e);    // a_{n+1} into ax, ay, az
    }

    /// @brief Advances a system by a number of Forward Euler steps
    /// @param bodies The system of stellar bodies (advanced in place)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size
    /// @param steps Number of steps
    /// @param force Force routine, see Force
    void euler_run(Object& bodies, double G, double e, double h, std::size_t steps, const Force& force) {
        if (!force) {
            throw std::invalid_argument("euler_run: no force routine given");
        }

        force(G, bodies, e);    // a_0
        for (std::size_t n = 0; n < steps; ++n) {
            euler_step(bodies, G, e, h, force);
        }
    }
}
