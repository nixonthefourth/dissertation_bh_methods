//
// Created by Mykyta Khomiakov on 19/11/2026.
//

#include "../include/helpers/verificators.h++"

#include <cmath>

namespace verify {
    /// @brief Kinetic energy T = sum_i m_i |v_i|^2 / 2 (Eq. (10))
    /// @param bodies The system of stellar bodies
    /// @details The notebook computes |p_i|^2 / (2 m_i) with p_i = m_i v_i, which is the same
    ///          quantity; m_i |v_i|^2 / 2 avoids the division (and a NaN for a zero mass).
    /// @return Returns T
    double t_abs(const Object& bodies) {
        return 0.5 * soa::weighted_norm_sq(bodies.velocity(), bodies.mass);
    }

    /// @brief Softened potential energy V = -G sum_{i<j} m_i m_j / sqrt(|r_ij|^2 + e^2) (Eq. (11))
    /// @param bodies The system of stellar bodies
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor (it has to be the same one as in acc, otherwise H is not the
    ///          energy that the integrator is actually conserving)
    /// @return Returns V
    double v_abs(const Object& bodies, double G, double e) {
        const std::size_t N = bodies.size();
        const double e2 = e * e;
        const soa::ConstVec3Span r = bodies.displacement();

        double v_tot = 0.0;     // potential energy total (without the -G factor)

        // Sum over every UNIQUE pair (i < j): a pair holds one potential energy, not two
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = i + 1; j < N; ++j) {
                v_tot += bodies.mass[i] * bodies.mass[j] / std::sqrt(norm_sq(r[i] - r[j]) + e2);
            }
        }

        return -G * v_tot;
    }

    ///@brief Calculate the absolute Hamiltonian value H = T + V (Eq. (12))
    ///@param v_abs System's potential energy (calculated previously)
    ///@param t_abs System's kinetic energy (calculated previously)
    ///@return Returns the total Hamiltonian
    double h_abs(const double v_abs, const double t_abs) { return v_abs + t_abs; }

    ///@brief Convert the Hamiltonian into a relative energy error H_rel (Eq. (15) with X = H)
    ///@param h_0 – initial Hamiltonian value at the start of the simulation clock
    ///@param  h_n – current Hamiltonian value at the given moment in time
    ///@return Returns the relative Hamiltonian value
    double h_rel(const double h_0, const double h_n) { return std::abs(h_n - h_0) / std::abs(h_0); }

    /// @brief Total linear momentum P = sum_i m_i v_i, a 3-vector (Eq. (13))
    /// @param bodies The system of stellar bodies
    /// @return Returns P
    Vec3 p_abs(const Object& bodies) {
        return soa::weighted_sum(bodies.velocity(), bodies.mass);
    }

    /// @brief Characteristic linear momentum of the system, sum_i m_i |v_i|
    /// @param bodies The system of stellar bodies
    /// @details Used as the reference P_ref of Eq. (15) when the total momentum is (close to) zero.
    /// @return Returns the scalar momentum scale
    double p_scale(const Object& bodies) {
        const soa::ConstVec3Span v = bodies.velocity();

        double scale = 0.0;
        for (std::size_t i = 0; i < bodies.size(); ++i) {
            scale += bodies.mass[i] * norm(v[i]);
        }

        return scale;
    }

    /// @brief Total angular momentum L = sum_i m_i (r_i x v_i), a 3-vector (Eq. (14))
    /// @param bodies The system of stellar bodies
    /// @details Cross product, NOT an element-wise product.
    /// @return Returns L
    Vec3 l_abs(const Object& bodies) {
        return soa::weighted_cross_sum(bodies.displacement(), bodies.velocity(), bodies.mass);
    }

    /// @brief Characteristic angular momentum of the system, sum_i m_i |r_i x v_i|
    /// @param bodies The system of stellar bodies
    /// @details Used as the reference L_ref of Eq. (15) when the total angular momentum is (close to) zero.
    /// @return Returns the scalar angular momentum scale
    double l_scale(const Object& bodies) {
        const soa::ConstVec3Span r = bodies.displacement();
        const soa::ConstVec3Span v = bodies.velocity();

        double scale = 0.0;
        for (std::size_t i = 0; i < bodies.size(); ++i) {
            scale += bodies.mass[i] * norm(cross(r[i], v[i]));
        }

        return scale;
    }

    /// @brief Relative error of a conserved VECTOR (Eq. (15) with X = P or X = L)
    /// @param x_0 Initial vector at the start of the simulation
    /// @param x_n Current vector at the time of the reading
    /// @param scale Characteristic size (p_scale or l_scale at t = 0), used instead of |x_0| if |x_0| ~ 0
    /// @details The length of the change divided by the length of the initial vector, or by the
    ///          characteristic size if the initial vector is zero (|x_0| < 1e-12 * scale).
    /// @return Returns the relative error
    double vec_rel(const Vec3& x_0, const Vec3& x_n, const double scale) {
        double reference = norm(x_0);
        if (reference < 1e-12 * scale) {    // X_0 = 0: the relative error is undefined, use the scale
            reference = scale;
        }

        return norm(x_n - x_0) / reference;
    }

    /// @brief Snapshot of all conserved quantities at one instant: returns (H, P, L)
    /// @param bodies The system of stellar bodies
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor (the same one as in acc)
    /// @return Returns the ConservationState {H, P, L}
    ConservationState conservation_state(const Object& bodies, const double G, const double e) {
        return {h_abs(v_abs(bodies, G, e), t_abs(bodies)), p_abs(bodies), l_abs(bodies)};
    }
}
