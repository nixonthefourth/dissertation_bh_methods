//
// Created by Mykyta Khomiakov on 19/11/2026.
//

#include "procgen/plummer.h++"

#include <cmath>
#include <numbers>
#include <stdexcept>

#include "helpers/minimal.h++"

namespace procgen {
    /// @brief Default Plummer scale length
    /// @param G Constant of gravitational attraction
    /// @param M Total mass of the cluster
    /// @details The total energy of a Plummer sphere is E = -3 pi G M^2 / (64 a), so E = -1/4 gives
    ///          a = 3 pi G M^2 / 16 (about 0.589 for G = M = 1).
    /// @return Returns the scale length a
    double plummer_scale_length(double G, double M) {
        return 3.0 * std::numbers::pi * G * M * M / 16.0;
    }

    /// @brief Plummer sphere sampler (Aarseth, Henon & Wielen 1974)
    /// @param N Number of bodies
    /// @param G Constant of gravitational attraction
    /// @param M Total mass of the cluster
    /// @param a Plummer scale length
    /// @param r_cut Radius cut-off in units of a (removes the very rare, extremely distant bodies)
    /// @details The model and the sampling method are described in the Plummer Sphere Explanation
    ///          (Eqs. (21) to (27)). Draws N bodies of mass M / N whose positions follow the Plummer
    ///          density profile
    ///              rho(r) = 3 M / (4 pi a^3) * (1 + r^2 / a^2)^(-5/2)      (Eq. (22))
    ///          and whose speeds follow the matching isotropic distribution function
    ///          f(E) ~ (-E)^(7/2) (Eq. (24)). The system is returned in the centre-of-mass frame
    ///          (P_0 = 0, Eq. (13)). The random numbers are drawn in the same order as in the notebook.
    /// @return Returns the system of bodies
    Object init_plummer(std::size_t N, double G, double M, double a, double r_cut) {
        if (N == 0) {
            throw std::invalid_argument("init_plummer: N must be positive");
        }
        if (a <= 0.0 || r_cut <= 0.0) {
            throw std::invalid_argument("init_plummer: a and r_cut must be positive");
        }

        Object bodies(N);
        const double m = M / static_cast<double>(N);

        for (std::size_t k = 0; k < N; ++k) {
            // Radius: invert the cumulative mass profile M(r) / M = X  ->  r = a / sqrt(X^(-2/3) - 1)
            // (Eqs. (23) and (27))
            double r;

            while (true) {
                const double X = minimal::random_uniform();
                if (X > 0.0) {
                    r = a / std::sqrt(std::pow(X, -2.0 / 3.0) - 1.0);
                    if (r < r_cut * a) {
                        break;
                    }
                }
            }

            // Speed: v = q * v_esc, with q drawn from g(q) = q^2 (1 - q^2)^(7/2) by rejection (Eq. (26)).
            // The maximum of g on [0, 1] is about 0.092, so a box of height 0.1 encloses it.
            double q;
            while (true) {
                q = minimal::random_uniform();
                const double y = 0.1 * minimal::random_uniform();
                if (y < q * q * std::pow(1.0 - q * q, 3.5)) {
                    break;
                }
            }
            const double v_esc = std::sqrt(2.0 * G * M / std::sqrt(r * r + a * a));    // local escape speed (Eq. (25))

            // Separate statements: the order of evaluation of function arguments is unspecified in
            // C++, so the position direction is drawn first explicitly, as in the notebook
            const Vec3 position = r * minimal::random_direction();
            const Vec3 velocity = q * v_esc * minimal::random_direction();

            bodies.add(m, position, velocity);
        }

        // Move to the centre-of-mass frame, so that the cluster does not drift and P_0 = 0 (up to rounding)
        minimal::to_com_frame(bodies, M);

        return bodies;
    }
}
