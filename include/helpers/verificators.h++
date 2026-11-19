//
// Created by Mykyta Khomiakov on 19/11/2026.
//

#ifndef MAIN_VERIFICATORS_H
#define MAIN_VERIFICATORS_H

#include <cstddef>
#include <cstdint>
#include <optional>
#include <random>

#include "Vec3.h++"
#include "soa.h++"
#include "object.h++"

namespace verify {
    /// @brief Kinetic energy T = sum_i m_i |v_i|^2 / 2 (Eq. (10))
    /// @param bodies The system of stellar bodies
    /// @details The notebook computes |p_i|^2 / (2 m_i) with p_i = m_i v_i, which is the same
    ///          quantity; m_i |v_i|^2 / 2 avoids the division (and a NaN for a zero mass).
    /// @return Returns T
    double t_abs(const Object& bodies);

    /// @brief Softened potential energy V = -G sum_{i<j} m_i m_j / sqrt(|r_ij|^2 + e^2) (Eq. (11))
    /// @param bodies The system of stellar bodies
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor (it has to be the same one as in acc, otherwise H is not the
    ///          energy that the integrator is actually conserving)
    /// @return Returns V
    double v_abs(const Object& bodies, double G, double e);

    ///@brief Calculate the absolute Hamiltonian value H = T + V (Eq. (12))
    ///@param v_abs System's potential energy (calculated previously)
    ///@param t_abs System's kinetic energy (calculated previously)
    ///@return Returns the total Hamiltonian
    double h_abs(double v_abs, double t_abs);

    ///@brief Convert the Hamiltonian into a relative energy error H_rel (Eq. (15) with X = H)
    ///@param h_0 – initial Hamiltonian value at the start of the simulation clock
    ///@param  h_n – current Hamiltonian value at the given moment in time
    ///@return Returns the relative Hamiltonian value
    double h_rel(const double h_0, const double h_n);

    /// @brief Total linear momentum P = sum_i m_i v_i, a 3-vector (Eq. (13))
    /// @param bodies The system of stellar bodies
    /// @return Returns P
    Vec3 p_abs(const Object& bodies);

    /// @brief Characteristic linear momentum of the system, sum_i m_i |v_i|
    /// @param bodies The system of stellar bodies
    /// @details Used as the reference P_ref of Eq. (15) when the total momentum is (close to) zero.
    /// @return Returns the scalar momentum scale
    double p_scale(const Object& bodies);

    /// @brief Total angular momentum L = sum_i m_i (r_i x v_i), a 3-vector (Eq. (14))
    /// @param bodies The system of stellar bodies
    /// @details Cross product, NOT an element-wise product.
    /// @return Returns L
    Vec3 l_abs(const Object& bodies);

    /// @brief Characteristic angular momentum of the system, sum_i m_i |r_i x v_i|
    /// @param bodies The system of stellar bodies
    /// @details Used as the reference L_ref of Eq. (15) when the total angular momentum is (close to) zero.
    /// @return Returns the scalar angular momentum scale
    double l_scale(const Object& bodies);

    /// @brief Relative error of a conserved VECTOR (Eq. (15) with X = P or X = L)
    /// @param x_0 Initial vector at the start of the simulation
    /// @param x_n Current vector at the time of the reading
    /// @param scale Characteristic size (p_scale or l_scale at t = 0), used instead of |x_0| if |x_0| ~ 0
    /// @details The length of the change divided by the length of the initial vector, or by the
    ///          characteristic size if the initial vector is zero (|x_0| < 1e-12 * scale).
    /// @return Returns the relative error
    double vec_rel(const Vec3& x_0, const Vec3& x_n, double scale);

    /// @brief Relative linear momentum error P_rel (Eq. (15) with X = P), the same function as vec_rel
    inline double p_rel(const Vec3& p_0, const Vec3& p_n, const double scale) { return vec_rel(p_0, p_n, scale); }

    /// @brief Relative angular momentum error L_rel (Eq. (15) with X = L), the same function as vec_rel
    inline double l_rel(const Vec3& l_0, const Vec3& l_n, const double scale) { return vec_rel(l_0, l_n, scale); }

    /// @brief Snapshot of all conserved quantities at one instant
    struct ConservationState {
        double H;   // Hamiltonian, Eq. (12)
        Vec3 P;     // total linear momentum, Eq. (13)
        Vec3 L;     // total angular momentum, Eq. (14)
    };

    /// @brief Snapshot of all conserved quantities at one instant: returns (H, P, L)
    /// @param bodies The system of stellar bodies
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor (the same one as in acc)
    /// @return Returns the ConservationState {H, P, L}
    ConservationState conservation_state(const Object& bodies, double G, double e);
}

#endif //MAIN_VERIFICATORS_H
