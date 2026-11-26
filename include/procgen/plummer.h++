//
// Created by Mykyta Khomiakov on 19/11/2026.
//

#ifndef MAIN_PLUMMER_H
#define MAIN_PLUMMER_H

#include <cstddef>

#include "helpers/object.h++"

/// @brief Procedural generation of initial conditions
namespace procgen {

    /// @brief Default Plummer scale length
    /// @param G Constant of gravitational attraction
    /// @param M Total mass of the cluster
    /// @details The total energy of a Plummer sphere is E = -3 pi G M^2 / (64 a), so E = -1/4 gives
    ///          a = 3 pi G M^2 / 16 (about 0.589 for G = M = 1).
    /// @return Returns the scale length a
    double plummer_scale_length(double G, double M);

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
    Object init_plummer(std::size_t N, double G, double M, double a, double r_cut = 10.0);

    /// @brief Uniform random sphere with unequal masses (the notebook's init_random_sphere)
    /// @param N Number of bodies
    /// @param G Constant of gravitational attraction
    /// @param M Total mass of the system
    /// @param R Radius of the sphere that the bodies are scattered in
    /// @param e Softening factor (needed for the potential energy of the virial rescaling)
    /// @param delta Relative spread of the masses, m_i uniform in [1 - delta, 1 + delta] before
    ///        normalisation (0 <= delta < 1, so that every mass is positive)
    /// @param q_vir Target virial ratio Q = T / |V|
    /// @details Unlike init_plummer, the bodies have unequal masses and are distributed uniformly inside
    ///          a sphere of radius R rather than following the Plummer profile. The raw masses are drawn
    ///          first and normalised so that sum_i m_i = M exactly. Each body then receives a radius
    ///          r = R X^(1/3) (uniform in volume), a random speed in [0, 1) and isotropic directions for
    ///          its displacement and velocity. The system is moved to the centre-of-mass frame and its
    ///          velocities are rescaled to the virial ratio q_vir. The random numbers are drawn in the
    ///          same order as in the notebook.
    /// @return Returns the system of bodies
    Object init_random_sphere(std::size_t N, double G, double M, double R, double e, double delta, double q_vir);
}

#endif //MAIN_PLUMMER_H
