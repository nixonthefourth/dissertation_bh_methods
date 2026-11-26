//
// Created by Mykyta Khomiakov on 25/11/2026.
//

#include "integrators/forest_ruth.h++"

#include <stdexcept>

namespace integrators {

    /// @brief One fourth-order Forest-Ruth step
    /// @param bodies The system of stellar bodies; on entry ax, ay, az must hold a(r_n)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size (may be negative)
    /// @param force Force routine, see Force
    void fr4_step(Object& bodies, double G, double e, double h, const Force& force) {
        for (std::size_t k = 0; k < forest_ruth::c.size(); ++k) {
            soa::axpy(bodies.velocity(), forest_ruth::c[k] * h, bodies.acceleration());    // kick:  v += c_k h a

            // Stage 4 has d_4 = 0: no drift, so the current accelerations stay valid
            if (forest_ruth::d[k] == 0.0) {
                continue;
            }

            soa::axpy(bodies.displacement(), forest_ruth::d[k] * h, bodies.velocity());    // drift: r += d_k h v
            force(G, bodies, e);    // a at the new displacements
        }
    }

    /// @brief Advances a system by a number of Forest-Ruth steps
    /// @param bodies The system of stellar bodies (advanced in place)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size
    /// @param steps Number of steps
    /// @param force Force routine, see Force
    void fr4_run(Object& bodies, double G, double e, double h, std::size_t steps, const Force& force) {
        if (!force) {
            throw std::invalid_argument("fr4_run: no force routine given");
        }

        force(G, bodies, e);    // a_0
        for (std::size_t n = 0; n < steps; ++n) {
            fr4_step(bodies, G, e, h, force);
        }
    }
}
