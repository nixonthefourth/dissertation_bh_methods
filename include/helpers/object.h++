//
// Created by Mykyta Khomiakov on 19/11/2026.
//

#ifndef MAIN_OBJECT_H
#define MAIN_OBJECT_H

#include <array>
#include <cstddef>
#include <vector>

#include "Vec3.h++"
#include "soa.h++"

/// @brief The whole system of bodies, stored as a fully split Structure of Arrays (SoA)
/// @details One Object holds ALL bodies of a simulation. Every component has its own array, and body i
///          is the i-th entry of each of them:
///              - displacement  dx[i], dy[i], dz[i]
///              - velocity      vx[i], vy[i], vz[i]
///              - acceleration  ax[i], ay[i], az[i]
///              - mass          mass[i]
///
///          This replaces the notebook's list of Body instances (an Array of Structures), so that each
///          component sits contiguously in memory and the force and integrator loops stream through it.
///
///          To avoid handling the three axes by hand, views are used:
///              - displacement(), velocity(), acceleration()  as  soa::Vec3Span (mutable)
///              - the same on a const Object                  as  soa::ConstVec3Span (read-only)
///          e.g. obj.displacement()[i] reads r_i as a Vec3, obj.velocity().set(i, v) writes v_i, and
///          soa::axpy(obj.displacement(), h, obj.velocity()) performs r += h v for all bodies.
///
///          All ten arrays must always have the same length, so bodies are only added with add() and
///          the size is only changed with resize(). Views become invalid after either of these calls.
class Object {
    public:
        std::vector<double> dx, dy, dz;     // displacement components
        std::vector<double> vx, vy, vz;     // velocity components
        std::vector<double> ax, ay, az;     // acceleration components, filled in by acc()
        std::vector<double> mass;           // masses (scalars)

        Object() = default;

        /// @brief Creates an empty system with memory reserved for n bodies
        explicit Object(std::size_t n) { reserve(n); }

        /// @brief Number of bodies in the system
        [[nodiscard]] std::size_t size() const { return mass.size(); }

        /// @brief True if the system holds no bodies
        [[nodiscard]] bool empty() const { return mass.empty(); }

        /// @brief Reserves memory for n bodies in every array (avoids reallocations while adding)
        void reserve(std::size_t n) {
            for (auto* component : components()) {
                component->reserve(n);
            }
        }

        /// @brief Resizes every array to n bodies (new bodies are zero, including their mass)
        void resize(std::size_t n) {
            for (auto* component : components()) {
                component->resize(n, 0.0);
            }
        }

        /// @brief Adds one body (the equivalent of the notebook's Body(m, r, v))
        /// @param m Mass (scalar)
        /// @param r Displacement vector
        /// @param v Velocity vector
        /// @details The acceleration starts at zero and is filled in by the integrator on its first call.
        void add(double m, const Vec3& r, const Vec3& v) {
            mass.push_back(m);
            dx.push_back(r.x); dy.push_back(r.y); dz.push_back(r.z);
            vx.push_back(v.x); vy.push_back(v.y); vz.push_back(v.z);
            ax.push_back(0.0); ay.push_back(0.0); az.push_back(0.0);
        }

        // Views of the vector quantities
        soa::Vec3Span displacement() { return {dx, dy, dz}; }
        soa::Vec3Span velocity() { return {vx, vy, vz}; }
        soa::Vec3Span acceleration() { return {ax, ay, az}; }

        [[nodiscard]] soa::ConstVec3Span displacement() const { return {dx, dy, dz}; }
        [[nodiscard]] soa::ConstVec3Span velocity() const { return {vx, vy, vz}; }
        [[nodiscard]] soa::ConstVec3Span acceleration() const { return {ax, ay, az}; }

        /// @brief Total mass of the system, sum_i m_i
        [[nodiscard]] double total_mass() const {
            double M = 0.0;
            for (const double m : mass) {
                M += m;
            }
            return M;
        }

        /// @brief True if all ten arrays have the same length (the SoA invariant)
        [[nodiscard]] bool is_consistent() const {
            const std::size_t n = mass.size();
            return dx.size() == n && dy.size() == n && dz.size() == n &&
                   vx.size() == n && vy.size() == n && vz.size() == n &&
                   ax.size() == n && ay.size() == n && az.size() == n;
        }

    private:
        /// @brief Pointers to all ten arrays, so that reserve and resize cannot miss one
        std::array<std::vector<double>*, 10> components() {
            return {&dx, &dy, &dz, &vx, &vy, &vz, &ax, &ay, &az, &mass};
        }
};

#endif //MAIN_OBJECT_H
