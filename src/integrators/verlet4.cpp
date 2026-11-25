//
// Created by Mykyta Khomiakov on 24/11/2026.
//

#include "integrators/verlet4.h++"

#include <stdexcept>

namespace integrators {

    /// @brief One second-order Velocity Verlet step in kick-drift-kick form
    /// @param bodies The system of stellar bodies; on entry ax, ay, az must hold a(r_n)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size (may be negative)
    /// @param force Force routine, see Force
    void vv_step(Object& bodies, double G, double e, double h, const Force& force) {
        soa::axpy(bodies.velocity(), 0.5 * h, bodies.acceleration());       // kick:  v += h/2 a_n
        soa::axpy(bodies.displacement(), h, bodies.velocity());             // drift: r += h v_{n+1/2}
        force(G, bodies, e);    // a_{n+1} into ax, ay, az
        soa::axpy(bodies.velocity(), 0.5 * h, bodies.acceleration());       // kick:  v += h/2 a_{n+1}
    }

    /// @brief One fourth-order Velocity Verlet step (Yoshida composition)
    /// @param bodies The system of stellar bodies; on entry ax, ay, az must hold a(r_n)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size
    /// @param force Force routine, see Force
    void vv4_step(Object& bodies, double G, double e, double h, const Force& force) {
        vv_step(bodies, G, e, yoshida::w1 * h, force);
        vv_step(bodies, G, e, yoshida::w0 * h, force);
        vv_step(bodies, G, e, yoshida::w1 * h, force);
    }

    /// @brief Advances a system by a number of fourth-order Velocity Verlet steps
    /// @param bodies The system of stellar bodies (advanced in place)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size
    /// @param steps Number of steps
    /// @param force Force routine, see Force
    void vv4_run(Object& bodies, double G, double e, double h, std::size_t steps, const Force& force) {
        if (!force) {
            throw std::invalid_argument("vv4_run: no force routine given");
        }

        force(G, bodies, e);    // a_0
        for (std::size_t n = 0; n < steps; ++n) {
            vv4_step(bodies, G, e, h, force);
        }
    }
}
