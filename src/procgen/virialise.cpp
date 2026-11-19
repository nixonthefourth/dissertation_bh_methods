//
// Created by Mykyta Khomiakov on 19/11/2026.
//

#include "procgen/virialise.h++"

#include <cmath>
#include <stdexcept>

#include "procgen/virialise.h++"
#include "helpers/verificators.h++"

namespace procgen {

    /// @brief Rescales all velocities so that the system has the virial ratio Q = T / |V|
    /// @param bodies The system of stellar bodies
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param q_vir Target virial ratio (default: equilibrium)
    /// @details T and V are those of Eq. (10) and Eq. (11). Q = 0.5 is exact virial equilibrium
    ///          (2T + V = 0); Q < 0.5 makes the system collapse. Random sampling only satisfies the
    ///          target ratio on average, and the rescaling removes the sampling noise. T scales with the
    ///          square of the velocities, so multiplying every v by sqrt(Q |V| / T) gives T = Q |V|
    ///          exactly. Applied after the centre-of-mass correction, P = 0 is preserved.
    ///          Used after init_plummer and by the random-sphere generator.
    void virialise(Object& bodies, double G, double e, double q_vir) {
        const double t = verify::t_abs(bodies);
        const double v = verify::v_abs(bodies, G, e);

        if (t <= 0.0) {
            throw std::invalid_argument("virialise: kinetic energy is zero, velocities cannot be rescaled");
        }

        soa::scale(bodies.velocity(), std::sqrt(q_vir * std::abs(v) / t));
    }
}
