//
// Created by Mykyta Khomiakov on 24/11/2026.
//

#ifndef MAIN_VERLET4_H
#define MAIN_VERLET4_H

#include <cstddef>
#include <functional>

#include "helpers/object.h++"
#include "helpers/soa.h++"

/// @brief Numerical Integration Schemes (symplectic and non-symplectic)
namespace integrators {

    /// @brief Force routine with the call form force(G, bodies, e)
    /// @details It must write the accelerations of the current displacements into the system's own
    ///          ax, ay, az. Both direct summation and the Barnes-Hut tree fit this form:
    ///              [](double G, Object& b, double e) { physics::acc(G, b, e); }    (a lambda, since acc is overloaded)
    ///              physics::acc_bh_theta(0.5)
    ///          The force is evaluated a few times per step only, so the cost of the std::function call
    ///          is negligible next to the O(N^2) or O(N log N) force evaluation itself.
    using Force = std::function<void(double G, Object& bodies, double e)>;

    /// @brief Coefficients of Yoshida's fourth-order triple-jump composition (Yoshida 1990)
    /// @details A step of size h is made of three Velocity Verlet sub-steps of sizes w1 h, w0 h, w1 h with
    ///              w1 = 1 / (2 - 2^(1/3)),    w0 = -2^(1/3) / (2 - 2^(1/3)).
    ///          They satisfy 2 w1 + w0 = 1 (the sub-steps add up to h) and 2 w1^3 + w0^3 = 0 (the
    ///          h^3 error terms of the second-order sub-steps cancel), and the sequence is symmetric, so
    ///          the composition is time-reversible and its error is O(h^4). w0 < 0: the middle sub-step
    ///          runs backwards in time.
    namespace yoshida {
        inline constexpr double cbrt2 = 1.2599210498948731647672106;    // 2^(1/3)
        inline constexpr double w1 = 1.0 / (2.0 - cbrt2);              //  1.3512071919596578
        inline constexpr double w0 = -cbrt2 / (2.0 - cbrt2);           // -1.7024143839193153
    }

    /// @brief One second-order Velocity Verlet step in kick-drift-kick form
    /// @param bodies The system of stellar bodies; on entry ax, ay, az must hold a(r_n)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size (may be negative)
    /// @param force Force routine, see Force
    /// @details v_{n+1/2} = v_n + h/2 a_n,   r_{n+1} = r_n + h v_{n+1/2},   a_{n+1} = a(r_{n+1}),
    ///          v_{n+1} = v_{n+1/2} + h/2 a_{n+1}.
    ///          Algebraically identical to r_{n+1} = r_n + h v_n + h^2/2 a_n, v_{n+1} = v_n + h/2 (a_n + a_{n+1}),
    ///          but it needs no scratch copy of a_n. One force evaluation per step; on exit ax, ay, az
    ///          hold a(r_{n+1}), ready for the next step.
    void vv_step(Object& bodies, double G, double e, double h, const Force& force);

    /// @brief One fourth-order Velocity Verlet step (Yoshida composition)
    /// @param bodies The system of stellar bodies; on entry ax, ay, az must hold a(r_n)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size
    /// @param force Force routine, see Force
    /// @details vv_step(w1 h) followed by vv_step(w0 h) and vv_step(w1 h). Each sub-step ends with the
    ///          acceleration that the next one starts from, so a full step costs three force evaluations.
    ///          The method is symplectic and time-reversible like Velocity Verlet, with a global error
    ///          of O(h^4). On exit ax, ay, az hold a(r_{n+1}).
    void vv4_step(Object& bodies, double G, double e, double h, const Force& force);

    /// @brief Advances a system by a number of fourth-order Velocity Verlet steps
    /// @param bodies The system of stellar bodies (advanced in place)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size
    /// @param steps Number of steps
    /// @param force Force routine, see Force
    /// @details Computes a(r_0) first, so that the stored accelerations need not be valid on entry
    ///          (e.g. straight after Object::add or minimal::copy_bodies), then makes the steps.
    ///          Total cost: 3 steps + 1 force evaluations.
    void vv4_run(Object& bodies, double G, double e, double h, std::size_t steps, const Force& force);
}

#endif //MAIN_VERLET4_H
