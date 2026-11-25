//
// Created by Mykyta Khomiakov on 25/11/2026.
//

#ifndef MAIN_LEAPFROG_H
#define MAIN_LEAPFROG_H
#include <functional>

#include "helpers/object.h++"

namespace integrators {
    /// @brief Force routine with the call form force(G, bodies, e)
    /// @details It must write the accelerations of the current displacements into the system's own
    ///          ax, ay, az. Both direct summation and the Barnes-Hut tree fit this form:
    ///              [](double G, Object& b, double e) { physics::acc(G, b, e); }    (a lambda, since acc is overloaded)
    ///              physics::acc_bh_theta(0.5)
    ///          The force is evaluated a few times per step only, so the cost of the std::function call
    ///          is negligible next to the O(N^2) or O(N log N) force evaluation itself.
    using Force = std::function<void(double G, Object& bodies, double e)>;

    /// @brief One second-order Leapfrog step in kick-drift-kick form
    /// @param bodies The system of stellar bodies; on entry ax, ay, az must hold a(r_n)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size (may be negative)
    /// @param force Force routine, see Force
    void leapfrog_step(Object& bodies, double G, double e, double h, const Force& force);

    /// @brief Advances a system by a number of second-order Yoshida steps
    /// @param bodies The system of stellar bodies (advanced in place)
    /// @param G Constant of gravitational attraction
    /// @param e Softening factor
    /// @param h Step size
    /// @param steps Number of steps
    /// @param force Force routine, see Force
    void leapfrog_run(Object& bodies, double G, double e, double h, std::size_t steps, const Force& force);
}

#endif //MAIN_LEAPFROG_H
