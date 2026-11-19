//
// Created by Mykyta Khomiakov on 19/11/2026.
//

#include "helpers/minimal.h++"

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace minimal {
    /// @brief Shared random number generator
    /// @details A single generator for the whole program (function-local static, created on first use),
    ///          seeded with default_seed so that every run is reproducible. A random_device seed would
    ///          give different initial conditions on every run, which breaks the comparisons.
    /// @return Returns a reference to the generator
    std::mt19937_64& rng() {
        static std::mt19937_64 gen(default_seed);
        return gen;
    }

    /// @brief Re-seeds the shared random number generator
    /// @param seed New seed
    void set_seed(std::uint64_t seed) {
        rng().seed(seed);
    }

    /// @brief Uniform random number in [0, 1), the equivalent of rnd.random()
    /// @return Returns the random number
    double random_uniform() {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        return distribution(rng());
    }

    /// @brief Random unit vector distributed uniformly over the surface of a sphere
    /// @return Returns the random vector in the type of Vec3
    Vec3 random_direction() {
        // Two separate statements fix the order of the draws (cos_t first, then phi)
        const double cos_t = 2.0 * random_uniform() - 1.0;                 // cos(polar angle) uniform in [-1, 1]
        const double phi = 2.0 * std::numbers::pi * random_uniform();      // azimuth uniform in [0, 2 pi)
        const double sin_t = std::sqrt(1.0 - cos_t * cos_t);

        return {sin_t * std::cos(phi), sin_t * std::sin(phi), cos_t};
    }

    /// @brief Independent copy of a system, so that every run starts from identical conditions
    /// @param bodies The system of stellar bodies
    /// @details Matches the notebook, where Body(m, r, v) starts with a = 0: masses, displacements and
    ///          velocities are copied, and the accelerations are reset to zero.
    /// @return Returns the copied system
    Object copy_bodies(const Object& bodies) {
        Object copy = bodies;
        soa::fill(copy.acceleration(), {});
        return copy;
    }

    /// @brief Moves a system to its centre-of-mass frame (in place)
    /// @param bodies The system of stellar bodies
    /// @param M Total mass of the system. Defaults to the sum of the masses; the generators pass the
    ///          nominal total mass instead, which differs from the sum only by rounding.
    /// @details Subtracts the mass-weighted mean displacement and velocity, so that the system does not
    ///          drift and its total linear momentum is zero (up to rounding).
    void to_com_frame(Object& bodies, std::optional<double> M) {
        const double M_tot = M.value_or(bodies.total_mass());
        if (M_tot == 0.0) {
            throw std::invalid_argument("to_com_frame: total mass is zero");
        }

        const Vec3 r_cm = soa::weighted_sum(bodies.displacement(), bodies.mass) / M_tot;
        const Vec3 v_cm = soa::weighted_sum(bodies.velocity(), bodies.mass) / M_tot;

        soa::translate(bodies.displacement(), -r_cm);
        soa::translate(bodies.velocity(), -v_cm);
    }

    /// @brief Euclidean distance between two bodies
    /// @param bodies The system of stellar bodies
    /// @param i Index of body i
    /// @param j Index of body j
    /// @return Returns |r_i - r_j|
    double euclidean_distance(const Object& bodies, std::size_t i, std::size_t j) {
        const soa::ConstVec3Span r = bodies.displacement();
        return norm(r[i] - r[j]);
    }

    /// @brief Randomised initial speed
    /// @param bodies The system of stellar bodies
    /// @param i Index of the body whose mass sets the speed scale
    /// @param G Constant of gravitational attraction
    /// @param distance Previously found Euclidean distance
    /// @details The random factor in [0, 1) gives random orbit eccentricities below the circular
    ///          speed sqrt(G m / d), so that no body is ejected.
    /// @return Returns the scalar speed
    double vel_init(const Object& bodies, std::size_t i, double G, double distance) {
        return random_uniform() * std::sqrt((G * bodies.mass[i]) / distance);
    }

    /// @brief Gets a vector perpendicular to v
    /// @param v Input vector
    /// @details (-y, x, 0) is perpendicular to (x, y, z); if v lies on the z axis, (0, 1, 0) is used.
    ///          The result is NOT normalised.
    /// @return Returns the perpendicular vector
    Vec3 get_perpendicular(const Vec3& v) {
        if (v.x != 0.0 || v.y != 0.0) {
            return {-v.y, v.x, 0.0};
        }

        return {0.0, 1.0, 0.0};
    }

    /// @brief Converts a scalar speed into a velocity vector
    /// @param bodies The system of stellar bodies
    /// @param i Index of body i
    /// @param j Index of body j
    /// @param v Scalar speed
    /// @return Returns v times the unit vector perpendicular to r_j - r_i
    Vec3 vel_direction(const Object& bodies, std::size_t i, std::size_t j, double v) {
        const soa::ConstVec3Span r = bodies.displacement();
        const Vec3 perp_vector = get_perpendicular(r[j] - r[i]);

        return v * (perp_vector / norm(perp_vector));   // normalise by its own length
    }
}
