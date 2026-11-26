//
// Created by Mykyta Khomiakov on 25/11/2026.
//

#ifndef MAIN_FOREST_RUTH_H
#define MAIN_FOREST_RUTH_H

#include <array>
#include <cstddef>

#include "helpers/object.h++"
#include "integrators/force.h++"

namespace integrators {

    /// @brief Coefficients of the Forest-Ruth fourth-order scheme (Forest & Ruth 1990)
    /// @details chi is the real root of 48 chi^3 + 24 chi^2 - 1 = 0 (called x by Forest & Ruth),
    ///              chi = (2^(1/3) + 2^(-1/3) - 1) / 6.
    ///          A step consists of four kick-drift stages with the kick coefficients c_1..c_4 and the
    ///          drift coefficients d_1..d_4; each set sums to 1, and d_4 = 0.
    namespace forest_ruth {
        inline constexpr double cbrt2 = 1.2599210498948731647672106;        // 2^(1/3)
        inline constexpr double inv_cbrt2 = 0.7937005259840997373758528;    // 2^(-1/3)
        inline constexpr double chi = (cbrt2 + inv_cbrt2 - 1.0) / 6.0;      //  0.1756035959798288

        inline constexpr std::array<double, 4> c{chi + 0.5, -chi, -chi, chi + 0.5};                    // kicks
        inline constexpr std::array<double, 4> d{2.0 * chi + 1.0, -4.0 * chi - 1.0, 2.0 * chi + 1.0, 0.0};  // drifts
    }

    /// @brief One fourth-order Forest-Ruth step
    /// @param bodies The system of stellar bodies; on entry ax, ay, az must hold a(r_n)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size (may be negative)
    /// @param force Force routine, see Force
    /// @details For k = 1..4:  v += c_k h a,  r += d_k h v,  a = a(r).
    ///          Stage 4 has d_4 = 0: the bodies do not move, so its drift and force evaluation are
    ///          skipped and the current accelerations are reused by the first kick of the next step.
    ///          Symplectic and time-reversible, global error O(h^4). Three force evaluations per step;
    ///          on exit ax, ay, az hold a(r_{n+1}).
    void fr4_step(Object& bodies, double G, double e, double h, const Force& force);

    /// @brief Advances a system by a number of Forest-Ruth steps
    /// @param bodies The system of stellar bodies (advanced in place)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size
    /// @param steps Number of steps
    /// @param force Force routine, see Force
    /// @details Computes a(r_0) first, so that the stored accelerations need not be valid on entry.
    ///          Total cost: 3 steps + 1 force evaluations.
    void fr4_run(Object& bodies, double G, double e, double h, std::size_t steps, const Force& force);
}

#endif //MAIN_FOREST_RUTH_H
