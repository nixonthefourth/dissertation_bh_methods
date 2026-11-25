cmake_minimum_required(VERSION 4.4)
project(main)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(
        main main.cpp

        src/helpers/minimal.cpp
        include/helpers/Vec3.h++
        include/helpers/soa.h++
        include/helpers/object.h++
        include/helpers/minimal.h++

        src/procgen/plummer.cpp
        src/procgen/virialise.cpp
        include/procgen/plummer.h++
        include/procgen/virialise.h++
        src/physics/accelerations.cpp
        include/physics/accelerations.h++
        src/helpers/verificators.cpp
        include/helpers/verificators.h++
        include/helpers/verificators.h++
        src/physics/barnes_hut.cpp
        include/physics/barnes_hut.h++
        include/integrators/verlet4.h++
        src/integrators/verlet4.cpp
        src/integrators/leapfrog.cpp
        include/integrators/leapfrog.h++
)

# Lets sources include headers as "helpers/minimal.h++" or "procgen/plummer.h++"
target_include_directories(main PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
