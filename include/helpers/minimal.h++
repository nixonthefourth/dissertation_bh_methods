//
// Created by Mykyta Khomiakov on 19/11/2026.
//

#ifndef MAIN_MINIMAL_H
#define MAIN_MINIMAL_H

#include <cstddef>
#include <cstdint>
#include <optional>
#include <random>

#include "Vec3.h++"
#include "soa.h++"
#include "object.h++"

/// @brief Helper functions shared by the experiments, the initial-condition generators and the
///        force evaluation (port of the notebook's helper cell).
namespace minimal {

    /// @brief Default randomisation seed, the same one as in the notebook
    inline constexpr std::uint64_t default_seed = 256;

    /// @brief Shared random number generator used by every helper that draws random numbers
    /// @details Seeded with default_seed, so that runs are reproducible (as with rnd.seed(seed)).
    std::mt19937_64& rng();

    /// @brief Re-seeds the shared random number generator
    void set_seed(std::uint64_t seed);

    /// @brief Uniform random number in [0, 1), the equivalent of rnd.random()
    double random_uniform();

    /// @brief Random unit vector, uniformly distributed over the surface of a sphere
    Vec3 random_direction();

    /// @brief Independent copy of a system, so that every run starts from identical conditions
    Object copy_bodies(const Object& bodies);

    /// @brief Moves a system to its centre-of-mass frame (in place)
    void to_com_frame(Object& bodies, std::optional<double> M = std::nullopt);

    /// @brief Euclidean distance between bodies i and j
    double euclidean_distance(const Object& bodies, std::size_t i, std::size_t j);

    /// @brief Randomised initial speed of body i from its mass and a distance
    double vel_init(const Object& bodies, std::size_t i, double G, double distance);

    /// @brief A vector perpendicular to v
    Vec3 get_perpendicular(const Vec3& v);

    /// @brief Converts a scalar speed into a velocity vector perpendicular to r_j - r_i
    Vec3 vel_direction(const Object& bodies, std::size_t i, std::size_t j, double v);
}

#endif //MAIN_MINIMAL_H
