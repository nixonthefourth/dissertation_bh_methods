//
// Created by Mykyta Khomiakov on 19/11/2026.
//

#ifndef MAIN_ACCELERATIONS_H
#define MAIN_ACCELERATIONS_H

#include <span>

#include "helpers/Vec3.h++"
#include "helpers/soa.h++"
#include "helpers/object.h++"

/// @brief Force evaluation
namespace physics {
    /// @brief Softened gravitational acceleration of every body by direct summation (SoA kernel)
    /// @param G Constant of gravitational attraction
    /// @param r Displacements of the N bodies, as three arrays x, y, z
    /// @param m Masses of the N bodies
    /// @param e Softening factor to prevent division by zero if bodies get too close
    /// @param out Output view of size N; out[i] = -G sum_{j != i} m_j r_ij / (|r_ij|^2 + e^2)^(3/2) (Eq. (4))
    /// @details This is the core of the direct summation and the only place where the pair loop is
    ///          written. It works on bare SoA arrays and does not need an Object, so it can evaluate
    ///          the accelerations of any set of positions, e.g. a trial configuration held in a
    ///          soa::Vec3Array scratch buffer, without copying it into a system first.
    ///          Each unique pair (i < j) is evaluated once and applied to both bodies (Newton's third
    ///          law), which halves the work of the notebook's double loop and conserves the total
    ///          momentum to rounding. The result agrees with the notebook's acc up to floating-point
    ///          rounding (the summation order differs).
    ///          Throws std::invalid_argument if the sizes do not match, or if out overlaps r or m
    ///          (out is zeroed before the inputs are read, so an overlap would corrupt the result).
    void acc(double G, soa::ConstVec3Span r, std::span<const double> m, double e, soa::Vec3Span out);

    /// @brief Computes the accelerations of a system into a separate output view
    /// @param G Constant of gravitational attraction
    /// @param bodies The system of stellar bodies (displacements and masses are read)
    /// @param e Softening factor
    /// @param out Output view of size N, which must not be a view of the displacements or masses
    /// @details Equivalent to acc(G, bodies.displacement(), bodies.mass, e, out).
    void acc(double G, const Object& bodies, double e, soa::Vec3Span out);

    /// @brief Computes the accelerations and stores them in the system's ax, ay, az
    /// @param G Constant of gravitational attraction
    /// @param bodies The system of stellar bodies
    /// @param e Softening factor
    /// @details This is the call form force(G, bodies, e) expected by the integrators.
    void acc(double G, Object& bodies, double e);
}

#endif //MAIN_ACCELERATIONS_H
