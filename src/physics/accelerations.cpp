//
// Created by Mykyta Khomiakov on 19/11/2026.
//

#include "physics/accelerations.h++"

#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>

namespace physics {
    namespace {
        /// @brief True if the memory of the two arrays overlaps (std::less gives a total order on pointers)
        bool overlaps(std::span<const double> a, std::span<const double> b) {
            if (a.empty() || b.empty()) {
                return false;
            }
            const std::less<const double*> before;
            return before(a.data(), b.data() + b.size()) && before(b.data(), a.data() + a.size());
        }

        /// @brief True if any array of the output view overlaps any of the input arrays
        bool aliases(soa::ConstVec3Span out, soa::ConstVec3Span r, std::span<const double> m) {
            for (const auto o : {out.x, out.y, out.z}) {
                for (const auto in : {r.x, r.y, r.z, m}) {
                    if (overlaps(o, in)) {
                        return true;
                    }
                }
            }
            return false;
        }
    }

    /// @brief Softened gravitational acceleration of every body by direct summation (SoA kernel)
    /// @param G Constant of gravitational attraction
    /// @param r Displacements of the N bodies, as three arrays x, y, z
    /// @param m Masses of the N bodies
    /// @param e Softening factor to prevent division by zero if bodies get too close
    /// @param out Output view of size N; out[i] = -G sum_{j != i} m_j r_ij / (|r_ij|^2 + e^2)^(3/2) (Eq. (4))
    /// @details Each unique pair (i < j) is evaluated once and applied to both bodies (Newton's third
    ///          law). This is the hot loop, so it works on raw pointers to the component arrays: the
    ///          compiler then sees plain contiguous double arrays and can keep the body i values in
    ///          registers for the whole inner loop.
    void acc(double G, soa::ConstVec3Span r, std::span<const double> m, double e, soa::Vec3Span out) {
        const std::size_t N = r.size();
        if (m.size() != N) {
            throw std::invalid_argument("acc: there must be one mass per displacement");
        }
        if (out.size() != N) {
            throw std::invalid_argument("acc: output view must have one entry per body");
        }
        if (aliases(out, r, m)) {
            throw std::invalid_argument("acc: output view must not overlap the displacements or masses");
        }

        const double e2 = e * e;

        // Inputs
        const double* x = r.x.data();
        const double* y = r.y.data();
        const double* z = r.z.data();
        const double* mass = m.data();

        // Output
        double* out_x = out.x.data();
        double* out_y = out.y.data();
        double* out_z = out.z.data();

        soa::fill(out, {});

        // Sigma operation over every unique pair
        for (std::size_t i = 0; i < N; ++i) {
            const double x_i = x[i], y_i = y[i], z_i = z[i], m_i = mass[i];
            double a_x = 0.0, a_y = 0.0, a_z = 0.0;     // accumulated locally, written back once per i

            for (std::size_t j = i + 1; j < N; ++j) {
                // Displacement from j to i
                const double r_x = x_i - x[j];
                const double r_y = y_i - y[j];
                const double r_z = z_i - z[j];

                // Softened distance squared, and 1 / (|r_ij|^2 + e^2)^(3/2) without std::pow
                const double s = r_x * r_x + r_y * r_y + r_z * r_z + e2;
                const double inv_denominator = 1.0 / (s * std::sqrt(s));

                // Pull of j on i
                const double k_i = mass[j] * inv_denominator;
                a_x += k_i * r_x;
                a_y += k_i * r_y;
                a_z += k_i * r_z;

                // Equal and opposite pull of i on j
                const double k_j = m_i * inv_denominator;
                out_x[j] -= k_j * r_x;
                out_y[j] -= k_j * r_y;
                out_z[j] -= k_j * r_z;
            }

            out_x[i] += a_x;
            out_y[i] += a_y;
            out_z[i] += a_z;
        }

        // Apply -G to the net accelerations
        soa::scale(out, -G);
    }

    /// @brief Computes the accelerations of a system into a separate output view
    /// @param G Constant of gravitational attraction
    /// @param bodies The system of stellar bodies (displacements and masses are read)
    /// @param e Softening factor
    /// @param out Output view of size N, which must not be a view of the displacements or masses
    void acc(double G, const Object& bodies, double e, soa::Vec3Span out) {
        acc(G, bodies.displacement(), bodies.mass, e, out);
    }

    /// @brief Computes the accelerations and stores them in the system's ax, ay, az
    /// @param G Constant of gravitational attraction
    /// @param bodies The system of stellar bodies
    /// @param e Softening factor
    void acc(double G, Object& bodies, double e) {
        acc(G, std::as_const(bodies), e, bodies.acceleration());
    }
}
