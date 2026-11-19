//
// Created by Mykyta Khomiakov on 19/11/2026.
//

#ifndef MAIN_VIRIALISE_H
#define MAIN_VIRIALISE_H

#include "helpers/object.h++"

namespace procgen {

    /// @brief Rescales all velocities so that the system has the virial ratio Q = T / |V|
    void virialise(Object& bodies, double G, double e, double q_vir = 0.5);
}

#endif //MAIN_VIRIALISE_H
