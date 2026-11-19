//
// Created by Mykyta Khomiakov on 19/11/2026.
//

#ifndef MAIN_SOA_H
#define MAIN_SOA_H

#include <cassert>
#include <cstddef>
#include <span>
#include <vector>

#include "Vec3.h++"

/// @brief Helpers for processing vector quantities stored as Structure of Arrays (SoA)
/// @details A vector quantity of N bodies (displacement, velocity, acceleration) is stored as three
///          separate arrays x[0..N), y[0..N), z[0..N). Writing three loops by hand for every operation
///          is easy to get wrong (e.g. a stray "vy" in the z loop), so the helpers below take care of
///          the three axes in one place:
///            - Vec3Span / ConstVec3Span bundle the three arrays of one quantity into a single view,
///            - operator[] / set / add / sub read and write one body as a Vec3,
///            - per_axis() runs the same code once per axis,
///            - bulk operations (axpy, scale, fill, translate, weighted sums, ...) are built on per_axis.
///          Views do not own memory: they stay valid as long as the arrays are not resized.
namespace soa {

    /// @brief Read-only view of a vector quantity stored as three arrays x, y, z
    struct ConstVec3Span {
        std::span<const double> x, y, z;

        ConstVec3Span(std::span<const double> x, std::span<const double> y, std::span<const double> z)
            : x(x), y(y), z(z) {
            assert(x.size() == y.size() && y.size() == z.size());
        }

        /// @brief Number of bodies
        [[nodiscard]] std::size_t size() const { return x.size(); }

        /// @brief Reads body i as a Vec3
        Vec3 operator[](std::size_t i) const { return {x[i], y[i], z[i]}; }
    };

    /// @brief Mutable view of a vector quantity stored as three arrays x, y, z
    struct Vec3Span {
        std::span<double> x, y, z;

        Vec3Span(std::span<double> x, std::span<double> y, std::span<double> z)
            : x(x), y(y), z(z) {
            assert(x.size() == y.size() && y.size() == z.size());
        }

        /// @brief A mutable view can always be used where a read-only one is expected
        operator ConstVec3Span() const { return {x, y, z}; }

        /// @brief Number of bodies
        [[nodiscard]] std::size_t size() const { return x.size(); }

        /// @brief Reads body i as a Vec3 (a copy: write back with set, add or sub)
        Vec3 operator[](std::size_t i) const { return {x[i], y[i], z[i]}; }

        /// @brief Overwrites body i with v
        void set(std::size_t i, const Vec3& v) const { x[i] = v.x; y[i] = v.y; z[i] = v.z; }

        /// @brief Adds v to body i
        void add(std::size_t i, const Vec3& v) const { x[i] += v.x; y[i] += v.y; z[i] += v.z; }

        /// @brief Subtracts v from body i
        void sub(std::size_t i, const Vec3& v) const { x[i] -= v.x; y[i] -= v.y; z[i] -= v.z; }
    };

    /// @brief Owning SoA buffer of N vectors, for scratch storage such as the old accelerations of
    ///        the Velocity Verlet step (allocate once, reuse every step)
    struct Vec3Array {
        std::vector<double> x, y, z;

        Vec3Array() = default;
        explicit Vec3Array(std::size_t n) : x(n, 0.0), y(n, 0.0), z(n, 0.0) {}

        [[nodiscard]] std::size_t size() const { return x.size(); }
        void resize(std::size_t n) { x.resize(n); y.resize(n); z.resize(n); }

        /// @brief Mutable view of the buffer
        Vec3Span view() { return {x, y, z}; }

        /// @brief Read-only view of the buffer
        [[nodiscard]] ConstVec3Span view() const { return {x, y, z}; }
    };

    /// @brief Runs f once per axis
    /// @details f is called three times: with the x arrays of all spans, then the y arrays, then the
    ///          z arrays. Code written inside f therefore deals with ONE axis only and cannot mix them.
    ///          Example (r += h v for all bodies):
    ///              soa::per_axis([h](std::span<double> r, std::span<const double> v) {
    ///                  for (std::size_t i = 0; i < r.size(); ++i) r[i] += h * v[i];
    ///              }, displacement, ConstVec3Span(velocity));
    template <class F, class... Spans>
    void per_axis(F&& f, const Spans&... spans) {
        f(spans.x...);
        f(spans.y...);
        f(spans.z...);
    }

    /// @brief y += k * x for every body (the update used by r += v h and v += a h)
    inline void axpy(Vec3Span y, double k, ConstVec3Span x) {
        assert(y.size() == x.size());
        per_axis([k](std::span<double> yc, std::span<const double> xc) {
            for (std::size_t i = 0; i < yc.size(); ++i) {
                yc[i] += k * xc[i];
            }
        }, y, x);
    }

    /// @brief v *= k for every body
    inline void scale(Vec3Span v, double k) {
        per_axis([k](std::span<double> c) {
            for (double& value : c) {
                value *= k;
            }
        }, v);
    }

    /// @brief Sets every body to the same vector (e.g. fill(acceleration, {}) for zero)
    inline void fill(Vec3Span v, const Vec3& value) {
        for (std::size_t i = 0; i < v.size(); ++i) {
            v.set(i, value);
        }
    }

    /// @brief Adds the same offset to every body (e.g. translate(displacement, -r_cm))
    inline void translate(Vec3Span v, const Vec3& offset) {
        for (std::size_t i = 0; i < v.size(); ++i) {
            v.add(i, offset);
        }
    }

    /// @brief dst = src for every body (the views must have the same size)
    inline void copy(ConstVec3Span src, Vec3Span dst) {
        assert(src.size() == dst.size());
        per_axis([](std::span<const double> s, std::span<double> d) {
            for (std::size_t i = 0; i < s.size(); ++i) {
                d[i] = s[i];
            }
        }, src, dst);
    }

    /// @brief sum_i v_i
    inline Vec3 sum(ConstVec3Span v) {
        Vec3 total{};
        for (std::size_t i = 0; i < v.size(); ++i) {
            total += v[i];
        }
        return total;
    }

    /// @brief sum_i w_i v_i (e.g. total momentum with w = masses and v = velocities)
    inline Vec3 weighted_sum(ConstVec3Span v, std::span<const double> w) {
        assert(v.size() == w.size());
        Vec3 total{};
        for (std::size_t i = 0; i < v.size(); ++i) {
            total += w[i] * v[i];
        }
        return total;
    }

    /// @brief sum_i w_i |v_i|^2 (e.g. twice the kinetic energy with w = masses and v = velocities)
    inline double weighted_norm_sq(ConstVec3Span v, std::span<const double> w) {
        assert(v.size() == w.size());
        double total = 0.0;
        for (std::size_t i = 0; i < v.size(); ++i) {
            total += w[i] * norm_sq(v[i]);
        }
        return total;
    }

    /// @brief sum_i w_i (a_i x b_i) (e.g. total angular momentum with w = masses, a = r, b = v)
    inline Vec3 weighted_cross_sum(ConstVec3Span a, ConstVec3Span b, std::span<const double> w) {
        assert(a.size() == b.size() && b.size() == w.size());
        Vec3 total{};
        for (std::size_t i = 0; i < a.size(); ++i) {
            total += w[i] * cross(a[i], b[i]);
        }
        return total;
    }
}

#endif //MAIN_SOA_H
