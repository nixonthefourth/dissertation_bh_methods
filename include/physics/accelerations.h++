//
// Created by Mykyta Khomiakov on 19/11/2026.
//

#ifndef MAIN_ACCELERATIONS_H
#define MAIN_ACCELERATIONS_H

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <optional>
#include <random>

#include "../include/helpers/Vec3.h++"
#include "../include/helpers/soa.h++"
#include "../include/helpers/object.h++"

namespace physics {
    /// @brief Computes the softened gravitational acceleration of every body by direct summation
    /// @param G Constant of gravitational attraction
    /// @param bodies The system of stellar bodies (displacements and masses are read)
    /// @param e Softening factor to prevent division by zero if bodies get too close
    /// @param out Output view of size N; out[i] = -G sum_{j != i} m_j r_ij / (|r_ij|^2 + e^2)^(3/2)
    /// @details Each unique pair (i < j) is evaluated once and applied to both bodies (Newton's third
    ///          law), which halves the work of the notebook's double loop. The result agrees with the
    ///          notebook's acc up to floating-point rounding (the summation order differs).
    ///          This is the hot loop, so it works on the raw component arrays directly; out must not
    ///          be a view of the displacements or masses.
    void acc(double G, const Object& bodies, double e, soa::Vec3Span out);

    /// @brief Computes the accelerations and stores them in the system's ax, ay, az
    /// @param G Constant of gravitational attraction
    /// @param bodies The system of stellar bodies
    /// @param e Softening factor
    void acc(double G, Object& bodies, double e);
}

#endif //MAIN_ACCELERATIONS_H
