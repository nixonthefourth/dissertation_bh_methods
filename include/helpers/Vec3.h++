//
// Created by Mykyta Khomiakov on 19/11/2026.
//

#ifndef MAIN_VEC3_H
#define MAIN_VEC3_H

#include <cmath>
#include <ostream>
#include <stdexcept>

/// @brief Three-dimensional vector of doubles
/// @details Holds the x (horizontal), y (vertical) and z (depth) components. Used for displacement,
///          velocity and acceleration, so that the same code works for planar (z = 0) and fully
///          three-dimensional problems.
struct Vec3 {
    // Definition
    double x, y, z;

    // Constructors
    constexpr Vec3() : x(0.0), y(0.0), z(0.0) {}
    constexpr Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    /// @brief Scalar Multiplication
    /// @param scalar Constant for Multiplication
    /// @return Returns the multiplied Vec3 struct implicitly
    constexpr Vec3 operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

    /// @brief In-place scalar multiplication
    /// @param scalar Constant for Multiplication
    /// @return Returns the multiplied Vec3 struct in-place with the *= operator
    constexpr Vec3& operator*=(double scalar) {
        x *= scalar;
        y *= scalar;
        z *= scalar;

        return *this;
    }

    /// @brief Scalar Division
    /// @param scalar Constant for Division
    /// @details Checks for division by zero first, consistently with operator/=
    /// @return Returns the divided Vec3 struct implicitly
    Vec3 operator/(double scalar) const {
        if (scalar == 0.0) {
            throw std::invalid_argument("Vec3 Division by zero!");
        }

        return {x / scalar, y / scalar, z / scalar};
    }

    /// @brief Scalar In-Place Division
    /// @param scalar Constant for Division
    /// @details Checks for division by zero first before dividing the values
    /// @return Returns the divided Vec3 struct
    Vec3& operator/=(double scalar) {
        if (scalar == 0.0) {
            throw std::invalid_argument("Vec3 Division by zero!");
        }

        x /= scalar;
        y /= scalar;
        z /= scalar;

        return *this;
    }

    /// @brief Unary negation
    /// @return Returns (-x, -y, -z)
    constexpr Vec3 operator-() const {
        return {-x, -y, -z};
    }

    /// @brief Allows subtraction between two Vec3 structures
    /// @param other Other Vec3 struct taken by reference
    /// @return Returns the implicit Vec3 subtraction
    constexpr Vec3 operator-(const Vec3& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }

    /// @brief Allows for in-place subtraction of two Vec3 structures
    /// @param other Simply the other Vec3 struct taken by reference
    /// @return Returns the in-place Vec3 subtraction
    constexpr Vec3& operator-=(const Vec3& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;

        return *this;
    }

    /// @brief Allows for implicit addition of two Vec3 structs
    /// @param other Another Vec3 structure to add
    /// @return Returns the addition of two Vec3 structures implicitly
    constexpr Vec3 operator+(const Vec3& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }

    /// @brief Allows the inplace addition of two Vec3 structs
    /// @param other Another Vec3 structure
    /// @return Returns the addition of two Vec3 structures inplace
    constexpr Vec3& operator+=(const Vec3& other) {
        x += other.x;
        y += other.y;
        z += other.z;

        return *this;
    }

    /// @brief Component-wise equality (exact, mainly for tests)
    constexpr bool operator==(const Vec3& other) const = default;

    /// @brief Reverse operand positioning through accessing `*` operator.
    /// @param scalar Scalar multiplication constant.
    /// @param v Value by reference of the Vec3 structure.
    /// @return Returns multiplied value of Vec3 as `k * (x, y, z)`.
    friend constexpr Vec3 operator*(double scalar, const Vec3& v) {
        return {v.x * scalar, v.y * scalar, v.z * scalar};
    }

    /// @brief Stream output for Vec3, e.g. std::cout << position.
    /// @return Returns the stream with the vector written as (x, y, z)
    friend std::ostream& operator<<(std::ostream& os, const Vec3& v) {
        return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    }
};

/// @brief Dot product a . b (equivalent of np.dot)
constexpr double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/// @brief Cross product a x b (equivalent of np.cross)
constexpr Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}

/// @brief Squared Euclidean length |v|^2 (no square root, cheaper when only comparing lengths)
constexpr double norm_sq(const Vec3& v) {
    return dot(v, v);
}

/// @brief Euclidean length |v| (equivalent of np.linalg.norm)
inline double norm(const Vec3& v) {
    return std::sqrt(norm_sq(v));
}

#endif //MAIN_VEC3_H
