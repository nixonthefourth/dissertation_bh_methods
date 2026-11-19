//
// Created by Mykyta Khomiakov on 19/11/2026.
//

#ifndef MAIN_PLUMMER_H
#define MAIN_PLUMMER_H

#include <cstddef>

#include "helpers/object.h++"

/// @brief Procedural generation of initial conditions
namespace procgen {

    /// @brief Default Plummer scale length a = 3 pi G M^2 / 16, which gives the total energy
    ///        E = -1/4 (Henon units, used when the experiment file has "a": null)
    double plummer_scale_length(double G, double M);

    /// @brief Plummer sphere sampler (Aarseth, Henon & Wielen 1974), used by Experiment 2
    Object init_plummer(std::size_t N, double G, double M, double a, double r_cut = 10.0);
}

#endif //MAIN_PLUMMER_H
