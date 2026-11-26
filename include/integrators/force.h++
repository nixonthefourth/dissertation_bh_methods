//
// Created by Mykyta Khomiakov on 25/11/2026.
//

#ifndef MAIN_FORCE_H
#define MAIN_FORCE_H

#include <functional>

#include "helpers/object.h++"

/// @brief Numerical Integration Schemes (symplectic and non-symplectic)
namespace integrators {

    /// @brief Force routine with the call form force(G, bodies, e), shared by every integrator
    /// @details It must write the accelerations of the current displacements into the system's own
    ///          ax, ay, az. Both direct summation and the Barnes-Hut tree fit this form:
    ///              [](double G, Object& b, double e) { physics::acc(G, b, e); }    (a lambda, since acc is overloaded)
    ///              physics::acc_bh_theta(0.5)
    ///          The force is evaluated a few times per step only, so the cost of the std::function call
    ///          is negligible next to the O(N^2) or O(N log N) force evaluation itself.
    using Force = std::function<void(double G, Object& bodies, double e)>;
}

#endif //MAIN_FORCE_H
