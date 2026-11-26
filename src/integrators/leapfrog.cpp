//
// Created by Mykyta Khomiakov on 25/11/2026.
//

#include "integrators/leapfrog.h++"

#include <stdexcept>

namespace integrators {

    /// @brief One second-order Leapfrog step in kick-drift-kick form
    /// @param bodies The system of stellar bodies; on entry ax, ay, az must hold a(r_n)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size (may be negative)
    /// @param force Force routine, see Force
    void leapfrog_step(Object& bodies, double G, double e, double h, const Force& force) {
        soa::axpy(bodies.velocity(), 0.5 * h, bodies.acceleration());       // kick:  v += h/2 a_n
        soa::axpy(bodies.displacement(), h, bodies.velocity());             // drift: r += h v_{n+1/2}
        force(G, bodies, e);    // a_{n+1} into ax, ay, az
        soa::axpy(bodies.velocity(), 0.5 * h, bodies.acceleration());       // kick:  v += h/2 a_{n+1}
    }

    /// @brief Advances a system by a number of second-order Leapfrog steps
    /// @param bodies The system of stellar bodies (advanced in place)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size
    /// @param steps Number of steps
    /// @param force Force routine, see Force
    void leapfrog_run(Object& bodies, double G, double e, double h, std::size_t steps, const Force& force) {
        // Edge Case Handling
        if (!force) {
            throw std::invalid_argument("leapfrog_run: no force routine given");
        }

        // Runtime
        force(G, bodies, e);    // a_0
        for (std::size_t n = 0; n < steps; ++n) {
            leapfrog_step(bodies, G, e, h, force);
        }
    }
}
