# N-Body Simulation: C++ Framework

This repository holds the C++ implementation of the dissertation project *Comparative Reliability Performance of N-Body Simulation under Hamiltonian Framework Applying Barnes-Hut Trees*. It is the target of the transfer announced in the notebook prototype: every feature is first validated in Python and then ported here.

This document describes every module of the port and shows, with short examples, how each of its functions is used. Section 4 gives a complete program that touches every module; Sections 5 to 9 then treat the modules one at a time, from the vector type up to the Barnes-Hut octree, and Section 10 shows how they combine into a time integrator. Equation numbers, such as Eq. (19), refer to the numbering of the notebook.

All examples in this document were compiled with GCC 11.4 (`-std=c++2b -O2`) and run; the values quoted in the comments are those that were printed.

---

## Contents

1. [Status of the port](#1-status-of-the-port)
2. [Folder layout](#2-folder-layout)
3. [Building](#3-building)
4. [Quick start](#4-quick-start)
5. [Vectors: `Vec3`](#5-vectors-vec3)
6. [Structure of Arrays: `soa` and `Object`](#6-structure-of-arrays-soa-and-object)
7. [Shared helpers: `minimal`](#7-shared-helpers-minimal)
8. [Initial conditions: `procgen`](#8-initial-conditions-procgen)
9. [Conserved quantities: `verify`](#9-conserved-quantities-verify)
10. [Force evaluation: direct summation](#10-force-evaluation-direct-summation)
11. [Force evaluation: the Barnes-Hut octree](#11-force-evaluation-the-barnes-hut-octree)
12. [Combining the modules: a time integrator](#12-combining-the-modules-a-time-integrator)
13. [Correspondence with the notebook](#13-correspondence-with-the-notebook)
14. [Validation of the Barnes-Hut port](#14-validation-of-the-barnes-hut-port)
15. [Accuracy and cost of the Barnes-Hut method](#15-accuracy-and-cost-of-the-barnes-hut-method)
16. [Known limitations](#16-known-limitations)

---

## 1. Status of the port

| Component | Notebook | C++ | State |
|---|---|---|---|
| Vector and SoA helpers | `Body`, NumPy arrays | `Vec3`, `soa`, `Object` | ported |
| Shared helpers | helper cell | `minimal` | ported |
| Initial conditions | `init_plummer`, `virialise` | `procgen` | ported |
| | `init_random_sphere` | | not yet ported |
| Conserved quantities | conservation cell | `verify` | ported |
| Direct summation | `acc` | `physics::acc` | ported |
| Barnes-Hut octree | `acc_bh` and the `_bh` functions | `physics::Nodes`, `physics::acc_bh` | ported |
| Integrators | six `run_*` functions | | not yet ported |
| JSON experiment pipeline | loader, logger, plotter, runtime | | not yet ported |

`main.cpp` is still a placeholder; the modules above are compiled into the executable but are not yet driven by it. Section 12 shows how an integrator can be written on top of the existing modules in the meantime.

---

## 2. Folder layout

```
main/
├── CMakeLists.txt
├── main.cpp
├── include/
│   ├── helpers/    Vec3.h++  soa.h++  object.h++  minimal.h++  verificators.h++
│   ├── procgen/    plummer.h++  virialise.h++
│   └── physics/    accelerations.h++  barnes_hut.h++
└── src/
    ├── helpers/    minimal.cpp  verificators.cpp
    ├── procgen/    plummer.cpp  virialise.cpp
    └── physics/    accelerations.cpp  barnes_hut.cpp
```

Headers use the extension `.h++` and are included relative to `include/` (for example `#include "physics/barnes_hut.h++"`). `Vec3.h++`, `soa.h++` and `object.h++` are header-only; the remaining modules have a source file in `src/`.

| Header | Namespace | Contents |
|---|---|---|
| `helpers/Vec3.h++` | global | `Vec3`, `dot`, `cross`, `norm_sq`, `norm` |
| `helpers/soa.h++` | `soa` | SoA views, owning buffer and bulk operations |
| `helpers/object.h++` | global | `Object`, the system of bodies |
| `helpers/minimal.h++` | `minimal` | random numbers, copying, centre-of-mass frame, velocity set-up |
| `helpers/verificators.h++` | `verify` | energy, momentum and angular momentum, and their relative errors |
| `procgen/plummer.h++` | `procgen` | Plummer sphere sampler |
| `procgen/virialise.h++` | `procgen` | rescaling to a given virial ratio |
| `physics/accelerations.h++` | `physics` | direct summation |
| `physics/barnes_hut.h++` | `physics` | Barnes-Hut octree and its diagnostics |

---

## 3. Building

The project requires a C++23 compiler and CMake 4.4 or newer (see `CMakeLists.txt`). With CMake:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

The build type matters for timing: the figures in Section 15 were measured with optimisation (`-O2`), so a Debug build is not representative of them.

To compile a program of one's own against the framework without CMake, it suffices to add `include/` to the include path and to compile the six source files with it:

```
g++ -std=c++23 -O2 -Iinclude my_program.cpp src/*/*.cpp -o my_program
```

(With GCC 11, `-std=c++23` is spelled `-std=c++2b`.) Within the CMake project, a new source file is added to the list in `add_executable`.

---

## 4. Quick start

The following program generates a Plummer sphere in virial equilibrium, records its conserved quantities, and compares the Barnes-Hut forces with direct summation.

```cpp
#include <algorithm>
#include <iostream>

#include "helpers/object.h++"
#include "helpers/verificators.h++"
#include "procgen/plummer.h++"
#include "procgen/virialise.h++"
#include "physics/accelerations.h++"
#include "physics/barnes_hut.h++"

int main() {
    const double G = 1.0, M = 1.0, e = 0.05;      // Hénon units, softening 0.05
    const std::size_t N = 1000;

    // 1. Initial conditions: a Plummer sphere in virial equilibrium
    const double a = procgen::plummer_scale_length(G, M);
    Object bodies = procgen::init_plummer(N, G, M, a);
    procgen::virialise(bodies, G, e);              // Q = T / |V| = 0.5

    // 2. Conserved quantities at t = 0
    const verify::ConservationState s0 = verify::conservation_state(bodies, G, e);
    std::cout << "H_0 = " << s0.H << ", P_0 = " << s0.P << ", L_0 = " << s0.L << '\n';

    // 3. Forces by direct summation and by Barnes-Hut, compared body by body
    soa::Vec3Array a_direct(N), a_tree(N);
    physics::acc(G, bodies, e, a_direct.view());
    physics::acc_bh(G, bodies, e, 0.5, a_tree.view());

    const std::vector<double> err = physics::force_error(a_tree.view(), a_direct.view());
    std::cout << "largest relative force error: " << *std::ranges::max_element(err) << '\n';
}
```

Output (with the default seed):

```
H_0 = -0.264435, P_0 = (-3.26272e-17, -1.17049e-17, -1.17082e-18), L_0 = (-0.0136341, -0.0179452, -0.00113097)
largest relative force error: 0.0466041
```

The total momentum is zero up to rounding because `init_plummer` returns the system in its centre-of-mass frame.

---

## 5. Vectors: `Vec3`

Declared in `include/helpers/Vec3.h++`. `Vec3` is a plain aggregate of three `double` components `x`, `y` and `z`, used for single displacements, velocities and accelerations. Planar problems simply keep `z = 0`. Almost everything is `constexpr`.

| Operation | Meaning |
|---|---|
| `Vec3()`, `Vec3(x, y, z)` | zero vector; vector with the given components |
| `a + b`, `a - b`, `-a` | addition, subtraction, negation |
| `a * k`, `k * a`, `a / k` | scalar multiplication and division |
| `a += b`, `a -= b`, `a *= k`, `a /= k` | in-place forms |
| `a == b` | exact component-wise comparison (intended for tests) |
| `dot(a, b)`, `cross(a, b)` | scalar and vector product (`np.dot`, `np.cross`) |
| `norm_sq(a)`, `norm(a)` | squared length and length (`np.linalg.norm`) |
| `std::cout << a` | prints `(x, y, z)` |

Division by zero throws `std::invalid_argument` rather than producing infinities.

```cpp
#include "helpers/Vec3.h++"

Vec3 r{1.0, 2.0, 2.0};
const Vec3 v{0.0, 1.0, 0.0};

const Vec3 sum  = r + v;          // (1, 3, 2)
const Vec3 half = 0.5 * r;        // (0.5, 1, 1); r * 0.5 is equivalent
r += v;                           // (1, 3, 2)
r /= 2.0;                         // (0.5, 1.5, 1)

const double d  = dot(r, v);      // 1.5
const Vec3   c  = cross(r, v);    // (-1, 0, 0.5)
const double l2 = norm_sq(r);     // 3.5, no square root: cheaper when only comparing lengths
const double l  = norm(r);        // 1.87083

std::cout << -r << '\n';          // (-0.5, -1.5, -1)

constexpr Vec3 k = cross(Vec3{1, 0, 0}, Vec3{0, 1, 0});
static_assert(k == Vec3{0, 0, 1});   // evaluated at compile time
```

---

## 6. Structure of Arrays: `soa` and `Object`

The notebook stores a list of `Body` objects, each with its own position, velocity and acceleration (an Array of Structures). The C++ port instead stores every component of every body in its own contiguous array (a Structure of Arrays), so that the force and integrator loops stream through memory. The namespace `soa` supplies the tools that make this layout as convenient to use as a list of `Vec3`, and `Object` uses them to hold a whole system.

### 6.1 Views and buffers

Declared in `include/helpers/soa.h++`.

| Type | Description |
|---|---|
| `soa::Vec3Span` | Mutable view of three arrays `x`, `y`, `z` of equal length. `v[i]` reads body `i` as a `Vec3` (a copy); `v.set(i, u)`, `v.add(i, u)` and `v.sub(i, u)` write it. Converts implicitly to `ConstVec3Span`. |
| `soa::ConstVec3Span` | Read-only counterpart; only `v[i]` and `size()`. |
| `soa::Vec3Array` | Owning buffer of `n` vectors (three `std::vector<double>`). `view()` returns a `Vec3Span` (or a `ConstVec3Span` on a `const` buffer); `resize(n)` grows all three arrays together. |

Views do not own memory and remain valid only as long as the underlying arrays are not resized.

```cpp
#include "helpers/soa.h++"

soa::Vec3Array r(3);                  // three vectors, all zero
soa::Vec3Span rs = r.view();

rs.set(0, {1.0, 0.0, 0.0});
rs.set(1, {0.0, 2.0, 0.0});
rs.add(1, {0.0, 0.5, 0.0});           // r_1 += (0, 0.5, 0)
const Vec3 r1 = rs[1];                // (0, 2.5, 0), a copy: write back with set, add or sub
const std::size_t n = rs.size();      // 3

// A view can also wrap arrays that are owned elsewhere
std::vector<double> x{1, 2}, y{3, 4}, z{5, 6};
soa::Vec3Span w(x, y, z);
w.set(0, {0, 0, 0});                  // x[0] = y[0] = z[0] = 0
```

A `Vec3Array` must be sized before it is written to, either through the constructor `Vec3Array(n)` or through `resize(n)`. It is the intended type for scratch storage that is allocated once and reused at every step, such as the old accelerations of the Velocity Verlet step (Section 12).

### 6.2 Bulk operations and reductions

Every operation below acts on all bodies at once and handles the three axes internally.

| Function | Effect |
|---|---|
| `axpy(y, k, x)` | $y_i \mathrel{+}= k\,x_i$ (the update behind $r \mathrel{+}= h\,v$ and $v \mathrel{+}= h\,a$) |
| `scale(v, k)` | $v_i \mathrel{*}= k$ |
| `fill(v, u)` | $v_i = u$ (e.g. `fill(a, {})` for zero) |
| `translate(v, u)` | $v_i \mathrel{+}= u$ (e.g. `translate(r, -r_cm)`) |
| `copy(src, dst)` | $dst_i = src_i$ |
| `sum(v)` | $\sum_i v_i$ |
| `weighted_sum(v, w)` | $\sum_i w_i v_i$ (total momentum with $w=m$) |
| `weighted_norm_sq(v, w)` | $\sum_i w_i \lvert v_i\rvert^2$ (twice the kinetic energy with $w=m$) |
| `weighted_cross_sum(a, b, w)` | $\sum_i w_i\,(a_i\times b_i)$ (angular momentum with $w=m$, $a=r$, $b=v$) |
| `per_axis(f, spans...)` | Calls `f` three times: with the `x` arrays of all spans, then the `y` arrays, then the `z` arrays. |

```cpp
soa::Vec3Array r(3), v(3);
soa::fill(v.view(), {0.0, 1.0, 0.0});         // every v_i = (0, 1, 0)

const double h = 0.1;
soa::axpy(r.view(), h, v.view());             // r += h v
soa::scale(v.view(), 2.0);                    // v *= 2
soa::translate(r.view(), {-1.0, 0.0, 0.0});   // r -= (1, 0, 0)

soa::Vec3Array backup(3);
soa::copy(r.view(), backup.view());           // backup = r

const std::vector<double> m{1.0, 2.0, 3.0};
const Vec3   P  = soa::weighted_sum(v.view(), m);                   // total momentum
const double T2 = soa::weighted_norm_sq(v.view(), m);               // 2 T
const Vec3   L  = soa::weighted_cross_sum(r.view(), v.view(), m);   // total angular momentum
```

`per_axis` is the building block for operations that are not provided. The function object deals with one axis only and therefore cannot mix components by mistake (for example a stray `vy` in the `z` loop):

```cpp
// r -= h v, written by hand
soa::per_axis([h](std::span<double> rc, std::span<const double> vc) {
    for (std::size_t i = 0; i < rc.size(); ++i) {
        rc[i] -= h * vc[i];
    }
}, r.view(), v.view());
```

The function object may take its parameters as `std::span<double>` or `std::span<const double>`; a mutable span converts to a read-only one, so declaring an input as `const` documents that it is not modified.

### 6.3 The system of bodies: `Object`

Declared in `include/helpers/object.h++`. One `Object` holds all bodies of a simulation in ten arrays of equal length: `dx, dy, dz` (displacement), `vx, vy, vz` (velocity), `ax, ay, az` (acceleration) and `mass`. Body `i` is the `i`-th entry of each array.

| Member | Description |
|---|---|
| `Object()`, `Object(n)` | empty system; empty system with memory reserved for `n` bodies |
| `add(m, r, v)` | appends one body; its acceleration starts at zero (the notebook's `Body(m, r, v)`) |
| `size()`, `empty()` | number of bodies; whether there are none |
| `reserve(n)`, `resize(n)` | act on all ten arrays together; new bodies from `resize` are zero, mass included |
| `displacement()`, `velocity()`, `acceleration()` | `Vec3Span` views, or `ConstVec3Span` views on a `const Object` |
| `total_mass()` | $\sum_i m_i$ |
| `is_consistent()` | whether all ten arrays have the same length (the SoA invariant) |

```cpp
#include "helpers/object.h++"

Object sys(2);                                      // memory for two bodies
sys.add(1.0, {-0.5, 0.0, 0.0}, {0.0, -0.5, 0.0});   // mass, displacement, velocity
sys.add(1.0, { 0.5, 0.0, 0.0}, {0.0,  0.5, 0.0});

sys.size();                                         // 2
sys.total_mass();                                   // 2

// Three equivalent levels of access
sys.dx[0] += 0.1;                                   // one component of one body
const Vec3 r0 = sys.displacement()[0];              // one body as a Vec3: (-0.4, 0, 0)
sys.velocity().set(1, {0.0, 0.6, 0.0});             // write one body
soa::axpy(sys.displacement(), 0.01, sys.velocity()); // all bodies: r += h v

// A const Object gives read-only views
const Object& c = sys;
soa::ConstVec3Span a = c.acceleration();
```

Bodies should only be added with `add` and the size only changed with `resize`, so that the ten arrays cannot fall out of step; views taken before either call are invalidated by it. Writing directly into the public arrays is permitted as long as their lengths are left unchanged.

---

## 7. Shared helpers: `minimal`

Declared in `include/helpers/minimal.h++` (port of the notebook's helper cell).

### 7.1 Random numbers

All randomness in the framework is drawn from a single generator, `minimal::rng()`, a `std::mt19937_64` seeded with `minimal::default_seed` (256, as in the notebook). Every run of a program therefore produces the same initial conditions, which is what the comparisons between integrators require.

| Function | Description |
|---|---|
| `rng()` | reference to the shared generator, for use with any `<random>` distribution |
| `set_seed(seed)` | re-seeds the shared generator |
| `random_uniform()` | uniform number in $[0, 1)$ (`rnd.random()`) |
| `random_direction()` | unit vector uniformly distributed over the sphere |

```cpp
#include "helpers/minimal.h++"

minimal::set_seed(minimal::default_seed);     // restart the sequence (or any other seed)
const double u = minimal::random_uniform();   // 0.315187 for the default seed
const Vec3 n = minimal::random_direction();   // norm(n) == 1

// Any other distribution can draw from the same generator
std::normal_distribution<double> gauss(0.0, 1.0);
const double g = gauss(minimal::rng());
```

Because the generator is shared, the order of the draws matters: generating a second system, or calling `random_uniform` in between, changes all subsequent numbers. Calling `set_seed` immediately before a generator restores reproducibility.

### 7.2 Systems and frames

| Function | Description |
|---|---|
| `copy_bodies(bodies)` | independent copy with the accelerations reset to zero, so that every run starts from identical conditions |
| `to_com_frame(bodies, M = sum of masses)` | subtracts the centre-of-mass displacement and velocity in place, so that the system does not drift; throws if the total mass is zero |
| `euclidean_distance(bodies, i, j)` | $\lvert r_i - r_j\rvert$ |

```cpp
Object start = procgen::init_plummer(500, G, M, a);
Object run_a = minimal::copy_bodies(start);   // each run gets its own copy
Object run_b = minimal::copy_bodies(start);

minimal::to_com_frame(run_a);                 // P = 0 and centre of mass at the origin
minimal::to_com_frame(run_b, 1.0);            // with the nominal total mass instead of the sum
```

### 7.3 Setting up orbits

These three functions build an initial velocity for a body orbiting another one.

| Function | Description |
|---|---|
| `vel_init(bodies, i, G, d)` | random speed in $[0, \sqrt{G m_i / d})$, i.e. below the circular speed about body `i`, so that no body is ejected |
| `get_perpendicular(v)` | a vector perpendicular to `v` (not normalised): $(-v_y, v_x, 0)$, or $(0, 1, 0)$ if `v` lies on the $z$ axis |
| `vel_direction(bodies, i, j, v)` | velocity of magnitude `v` perpendicular to $r_j - r_i$ |

```cpp
const double G = 1.0;
Object pair;
pair.add(1.0,  {0.0, 0.0, 0.0}, {});          // central body
pair.add(1e-3, {1.0, 0.0, 0.0}, {});          // satellite, at rest for now

const double d = minimal::euclidean_distance(pair, 0, 1);   // 1
const double v = minimal::vel_init(pair, 0, G, d);          // random, below sqrt(G m_0 / d) = 1
pair.velocity().set(1, minimal::vel_direction(pair, 0, 1, v));   // tangential velocity

minimal::to_com_frame(pair);                  // remove the net drift of the pair
```

---

## 8. Initial conditions: `procgen`

Declared in `include/procgen/plummer.h++` and `include/procgen/virialise.h++`.

| Function | Description |
|---|---|
| `plummer_scale_length(G, M)` | $a = 3\pi G M^2/16$, which gives the total energy $E=-1/4$ (Hénon units); about 0.589 for $G=M=1$ |
| `init_plummer(N, G, M, a, r_cut = 10)` | $N$ bodies of mass $M/N$ sampled from the Plummer sphere (Aarseth, Hénon & Wielen, 1974; Eqs. (21) to (27)), returned in the centre-of-mass frame. Bodies beyond `r_cut` $\cdot\,a$ are redrawn. Throws `std::invalid_argument` if `N` is zero or `a`, `r_cut` are not positive. |
| `virialise(bodies, G, e, q_vir = 0.5)` | rescales all velocities so that $T/\lvert V\rvert$ equals `q_vir` exactly; $0.5$ is virial equilibrium, smaller values give a collapsing system. Throws if the kinetic energy is zero. |

```cpp
#include "procgen/plummer.h++"
#include "procgen/virialise.h++"

const double G = 1.0, M = 1.0, e = 0.05;
const double a = procgen::plummer_scale_length(G, M);         // 0.589049

Object cluster = procgen::init_plummer(500, G, M, a);         // default cut-off 10 a
Object compact = procgen::init_plummer(500, G, M, a, 3.0);    // no body beyond 3 a

procgen::virialise(cluster, G, e);          // T / |V| = 0.5: equilibrium
procgen::virialise(compact, G, e, 0.25);    // T / |V| = 0.25: cold, collapses
```

Sampling only meets the virial ratio on average; `virialise` removes the sampling noise. It should be given the same softening as the force routine, since $V$ depends on it. Because it scales the velocities uniformly, it preserves the zero total momentum established by `init_plummer`.

---

## 9. Conserved quantities: `verify`

Declared in `include/helpers/verificators.h++`. These functions measure how well a run conserves the Hamiltonian, the linear momentum and the angular momentum. They always use exact pair sums, also when the forces come from the Barnes-Hut tree.

| Function | Description |
|---|---|
| `t_abs(bodies)` | kinetic energy $T=\sum_i m_i\lvert v_i\rvert^2/2$ (Eq. (10)) |
| `v_abs(bodies, G, e)` | softened potential energy $V=-G\sum_{i<j}m_im_j/\sqrt{\lvert r_{ij}\rvert^2+e^2}$ (Eq. (11)), $\mathcal{O}(N^2)$ |
| `h_abs(V, T)` | Hamiltonian $H=T+V$ (Eq. (12)); note the order of the arguments |
| `p_abs(bodies)`, `l_abs(bodies)` | total linear and angular momentum, as `Vec3` (Eqs. (13), (14)) |
| `p_scale(bodies)`, `l_scale(bodies)` | characteristic sizes $\sum_i m_i\lvert v_i\rvert$ and $\sum_i m_i\lvert r_i\times v_i\rvert$ |
| `h_rel(H_0, H_n)` | relative energy error $\lvert H_n-H_0\rvert/\lvert H_0\rvert$ (Eq. (15)) |
| `vec_rel(X_0, X_n, scale)` | relative error of a conserved vector, $\lvert X_n-X_0\rvert/\lvert X_0\rvert$, with `scale` in place of $\lvert X_0\rvert$ when $\lvert X_0\rvert<10^{-12}\,$`scale` |
| `p_rel`, `l_rel` | `vec_rel` under the names used for $P$ and $L$ |
| `conservation_state(bodies, G, e)` | `ConservationState{H, P, L}` in one call |

The softening `e` must be the one used by the force routine; otherwise $H$ is not the quantity that the integrator actually conserves.

```cpp
#include "helpers/verificators.h++"

// Individual quantities
const double T = verify::t_abs(bodies);
const double V = verify::v_abs(bodies, G, e);
const double H = verify::h_abs(V, T);
const Vec3   P = verify::p_abs(bodies);
const Vec3   L = verify::l_abs(bodies);

// Monitoring a run: snapshot at t = 0, reference scales, snapshot later
const verify::ConservationState s0 = verify::conservation_state(bodies, G, e);
const double p_ref = verify::p_scale(bodies);
const double l_ref = verify::l_scale(bodies);

/* ... advance the system ... */

const verify::ConservationState sn = verify::conservation_state(bodies, G, e);
const double H_rel = verify::h_rel(s0.H, sn.H);
const double P_rel = verify::p_rel(s0.P, sn.P, p_ref);
const double L_rel = verify::l_rel(s0.L, sn.L, l_ref);
```

The scales matter for systems in their centre-of-mass frame. There $P_0$ is zero up to rounding, and dividing by $\lvert P_0\rvert$ would turn rounding noise into an arbitrarily large "relative" error; `vec_rel` detects this case and divides by the characteristic momentum instead. The scales should be taken at $t=0$ and kept fixed for the whole run.

---

## 10. Force evaluation: direct summation

Declared in `include/physics/accelerations.h++`.

| Function | Description |
|---|---|
| `acc(G, const Object& bodies, e, soa::Vec3Span out)` | Writes $a_i=-G\sum_{j\neq i}m_j\,r_{ij}/(\lvert r_{ij}\rvert^2+e^2)^{3/2}$ (Eq. (4)) into `out`, which must have one entry per body and must not alias the displacements or masses. Throws `std::invalid_argument` on a size mismatch. |
| `acc(G, Object& bodies, e)` | Writes the accelerations into the system's own `ax`, `ay`, `az`. |

Each unique pair is evaluated once and applied to both bodies (Newton's third law), which halves the work of the notebook's double loop and makes the total momentum conserved to rounding.

```cpp
#include "physics/accelerations.h++"

Object b;
b.add(1.0, {0, 0, 0}, {});
b.add(1.0, {1, 0, 0}, {});

physics::acc(1.0, b, 0.0);                 // into b.ax, b.ay, b.az
b.acceleration()[0];                       // (1, 0, 0)
b.acceleration()[1];                       // (-1, 0, 0)

soa::Vec3Array out(b.size());
physics::acc(1.0, b, 0.0, out.view());     // into a separate buffer; b is not modified
```

Since `acc` is overloaded, `physics::acc` cannot be passed by name where a callable is expected; a lambda selects the overload (see Section 12):

```cpp
auto direct = [](double G, Object& b, double e) { physics::acc(G, b, e); };
```

---

## 11. Force evaluation: the Barnes-Hut octree

Declared in `include/physics/barnes_hut.h++`. This module is the complete C++ counterpart of the notebook's Barnes-Hut cell and offers an $\mathcal{O}(N\log N)$ replacement for `physics::acc`.

### 11.1 The method in brief

Direct summation (Eq. (4)) evaluates every pair of bodies and costs $\mathcal{O}(N^2)$ per force evaluation. The Barnes-Hut method (Barnes & Hut, 1986) replaces a compact group of distant bodies by a single pseudo-body carrying their total mass at their centre of mass. The bodies are organised in an *octree*: a hierarchy of cubic cells, each of which can be split into eight equal octants.

**Construction.** The root cell is the smallest cube that encloses all bodies. A body descends from the root into the octant that contains it (Eq. (16)),

$$k=[x\geq c_x]+2\,[y\geq c_y]+4\,[z\geq c_z],\qquad k\in\{0,\dots,7\}$$

until it reaches an empty leaf. If the leaf is occupied, it is subdivided (the centre of child $k$ is shifted by a quarter of the parent side along each axis, Eq. (17)), the resident body is pushed down, and the descent continues. Cells at the depth `Nodes::max_depth` (32) are never subdivided, so that bodies at (almost) identical positions cannot cause unbounded recursion; they share a leaf.

**Mass moments.** Every cell stores its total mass $M_c$ and its centre of mass $\hat{R}_c$ (Eq. (18)),

$$M_c=\sum_{j\in c}m_j,\qquad \hat{R}_c=\frac{1}{M_c}\sum_{j\in c}m_j\,\hat{r}_j$$

computed from the leaves upwards.

**Tree walk.** The acceleration of body $i$ is obtained by walking the tree from the root. A cell of side $s=2\,\mathrm{half}$ at the distance $d=\lvert\hat{r}_i-\hat{R}_c\rvert$ is accepted as one softened point mass if (Eq. (19))

$$\frac{s}{d}<\theta,$$

and is otherwise opened and its children are examined. A cell that contains body $i$ is always opened, regardless of the criterion (Salmon & Warren, 1994). Leaves are summed exactly, skipping body $i$ itself. The opening angle $\theta$ controls the trade-off between accuracy and cost: $\theta=0$ opens every cell and reproduces direct summation, and $\theta\approx0.5$ to $1$ is the usual range. The softening is the Plummer softening of Eq. (4), identical to that of `physics::acc`.

The tree is rebuilt at every force evaluation, because the bodies move between evaluations.

### 11.2 Data layout

The whole tree is held by one `physics::Nodes` object as a Structure of Arrays, in the same way as `Object` holds all bodies. Node `n` is the `n`-th entry of every array, and nodes refer to each other by index, never by pointer, so the arrays may grow during construction. The root is node 0.

| Array | Type | Meaning |
|---|---|---|
| `centre` | `soa::Vec3Array` | geometric centre of the cell |
| `com` | `soa::Vec3Array` | centre of mass of the bodies in the cell, $\hat{R}_c$ |
| `half` | `std::vector<double>` | half of the side length of the cell |
| `mass` | `std::vector<double>` | total mass in the cell, $M_c$ |
| `first_child` | `std::vector<std::size_t>` | index of child 0, or `NONE` for a leaf |
| `first_body` | `std::vector<std::size_t>` | first body held by a leaf, or `NONE` if empty |
| `next_body` | `std::vector<std::size_t>` | next body in the same leaf; one entry per **body**, not per node |

Three design decisions follow from this layout.

- **Contiguous children.** The eight children of a node are allocated together, so child `o` of node `n` is node `first_child[n] + o`. Only one index per node is needed, and the octant number of Eq. (16) maps directly to a child index.
- **Parents precede children.** A child is always appended after its parent, so every child has a larger index than its parent. The bottom-up pass that computes the mass moments, and the depth computation of `stats()`, therefore reduce to a plain sweep over the arrays, with no recursion.
- **Leaf lists.** A leaf holds one body, except at `max_depth`, where the bodies of a leaf are chained through `next_body`. The tree walk always sums a leaf body by body and exactly, so sharing a leaf costs no accuracy; only the cost of that leaf grows with the number of its bodies.

The memory cost is 80 bytes per node (ten `double`-sized entries) and 8 bytes per body. For a Plummer sphere of 4000 bodies the tree has 15 529 nodes, about 3.9 nodes per body.

### 11.3 Interface

**Force routines**

| Function | Description |
|---|---|
| `acc_bh(G, const Object& bodies, e, theta, soa::Vec3Span out)` | Writes the Barnes-Hut accelerations into `out` (size $N$). The counterpart of `acc`; builds a fresh tree. Throws `std::invalid_argument` if the size of `out` is wrong or `theta` is negative. |
| `acc_bh(G, Object& bodies, e, theta = 0.5)` | Writes the accelerations into the system's `ax`, `ay`, `az`. |
| `acc_bh_theta(theta)` | Returns a callable `void(double G, Object& bodies, double e)` with `theta` fixed, for code that expects the three-argument form of a force routine (as the integrators will). |

**Diagnostics**

| Function | Description |
|---|---|
| `force_error(a_tree, a_direct)` | Relative force error $\lvert\hat{a}_i^{\mathrm{BH}}-\hat{a}_i^{\mathrm{direct}}\rvert/\lvert\hat{a}_i^{\mathrm{direct}}\rvert$ of every body. Throws if the two views differ in size. |
| `Nodes::stats()` | `TreeStats{nodes, leaves, depth}`: number of nodes, number of non-empty leaves and maximum depth. |

**The class `Nodes`**

| Member | Description |
|---|---|
| `build(bodies)` | Builds the tree for the current positions: bounding cube, insertion of all bodies, mass moments. Reuses the memory of the previous tree. |
| `acceleration(i, bodies, G, e, theta)` | Acceleration of body `i` from the tree walk. `const`, so concurrent calls for different bodies are safe. |
| `size()`, `is_leaf(n)`, `child(n, o)` | Number of nodes; whether node `n` is a leaf; index of child `o` of node `n`. |
| `octant(n, r)`, `contains(n, r)` | The octant of Eq. (16), and the containment test used by the opening rule. |
| `NONE`, `max_depth` | The "no child / no body" marker, and the depth limit (32). |

Construction helpers (`add`, `subdivide`, `append`, `insert`, `moments`) are private. All entities carry Doxygen comments in the header.

### 11.4 Usage

**Accelerations of a system, directly into its acceleration arrays:**

```cpp
#include "physics/barnes_hut.h++"

physics::acc_bh(G, bodies, e, 0.5);          // fills bodies.ax, ay, az
physics::acc_bh(G, bodies, e);               // the same: theta defaults to 0.5
```

**Comparing with direct summation:**

```cpp
soa::Vec3Array a_direct(bodies.size()), a_tree(bodies.size());
physics::acc(G, bodies, e, a_direct.view());
physics::acc_bh(G, bodies, e, 0.5, a_tree.view());

std::vector<double> error = physics::force_error(a_tree.view(), a_direct.view());
```

**Inspecting the tree:**

```cpp
physics::Nodes tree;
tree.build(bodies);
const physics::TreeStats s = tree.stats();   // s.nodes, s.leaves, s.depth
const double M = tree.mass[0];               // total mass, stored at the root
const Vec3 R = tree.com.view()[0];           // centre of mass of the whole system

// Descend to the leaf that holds body 0 and list the bodies in it
const Vec3 r0 = bodies.displacement()[0];
std::size_t n = 0;
while (!tree.is_leaf(n)) {
    n = tree.child(n, tree.octant(n, r0));
}
for (std::size_t j = tree.first_body[n]; j != physics::Nodes::NONE; j = tree.next_body[j]) {
    std::cout << "leaf body " << j << '\n';  // normally only body 0 itself
}
```

**Acceleration of a single body.** `Nodes::acceleration` returns the same value that `acc_bh` writes for that body:

```cpp
const Vec3 a7 = tree.acceleration(7, bodies, G, e, 0.5);
```

**Reusing one tree across steps.** `acc_bh` constructs a new `Nodes` object, and therefore reallocates its arrays, at every call. A program that evaluates forces many times can keep one tree and call `build` repeatedly, which reuses the memory:

```cpp
physics::Nodes tree;
soa::Vec3Span a = bodies.acceleration();
for (std::size_t step = 0; step < steps; ++step) {
    /* ... move the bodies ... */
    tree.build(bodies);
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        a.set(i, tree.acceleration(i, bodies, G, e, theta));
    }
}
```

**Fixing the opening angle for an integrator:**

```cpp
for (const double theta : {0.25, 0.5, 1.0}) {
    const auto force = physics::acc_bh_theta(theta);   // theta captured by value
    force(G, bodies, e);                               // same call form as acc(G, bodies, e)
}
```

---

## 12. Combining the modules: a time integrator

The integrators of the notebook are not yet ported, but the modules above already contain everything an integrator needs. The following sketch implements one Velocity Verlet step,

$$r_{n+1}=r_n+h\,v_n+\tfrac{h^2}{2}\,a_n,\qquad v_{n+1}=v_n+\tfrac{h}{2}\,(a_n+a_{n+1}),$$

for any force routine with the call form `force(G, bodies, e)`. It is an illustration of the interfaces, not part of the framework.

```cpp
#include "helpers/object.h++"
#include "helpers/minimal.h++"
#include "helpers/verificators.h++"
#include "procgen/plummer.h++"
#include "procgen/virialise.h++"
#include "physics/accelerations.h++"
#include "physics/barnes_hut.h++"

/// One Velocity Verlet step; a_old is scratch storage of size N, allocated once by the caller
template <class Force>
void vv_step(Object& b, double G, double e, double h, Force&& force, soa::Vec3Array& a_old) {
    soa::copy(b.acceleration(), a_old.view());                 // a_n
    soa::axpy(b.displacement(), h, b.velocity());              // r += h v
    soa::axpy(b.displacement(), 0.5 * h * h, a_old.view());    //    + h^2 / 2 a_n
    force(G, b, e);                                            // a_{n+1} into b.ax, ay, az
    soa::axpy(b.velocity(), 0.5 * h, a_old.view());            // v += h / 2 a_n
    soa::axpy(b.velocity(), 0.5 * h, b.acceleration());        //    + h / 2 a_{n+1}
}

/// Runs a copy of the system and returns the relative energy error at the end
template <class Force>
double run(Object b, double G, double e, double h, std::size_t steps, Force&& force) {
    soa::Vec3Array a_old(b.size());
    force(G, b, e);                                            // a_0
    const auto s0 = verify::conservation_state(b, G, e);
    for (std::size_t n = 0; n < steps; ++n) {
        vv_step(b, G, e, h, force, a_old);
    }
    return verify::h_rel(s0.H, verify::conservation_state(b, G, e).H);
}

int main() {
    const double G = 1.0, M = 1.0, e = 0.05, h = 1e-3;
    Object start = procgen::init_plummer(500, G, M, procgen::plummer_scale_length(G, M));
    procgen::virialise(start, G, e);

    auto direct = [](double G, Object& b, double e) { physics::acc(G, b, e); };
    run(minimal::copy_bodies(start), G, e, h, 1000, direct);                       // 1.5e-7
    run(minimal::copy_bodies(start), G, e, h, 1000, physics::acc_bh_theta(1.0));   // 4.8e-4
}
```

The relative energy errors after 1000 steps are given in the comments. Because both force routines share one call form, the comparison between direct summation and the tree is a matter of passing a different callable, and `copy_bodies` guarantees that both runs start from the same state. The larger error of the tree code reflects the force errors of Section 15 and the lack of momentum conservation noted in Section 16.

---

## 13. Correspondence with the notebook

**Helpers, initial conditions and conserved quantities**

| Notebook | C++ | Remark |
|---|---|---|
| `class Body` (list of objects) | `Object` (one object for all bodies) | SoA instead of a list of bodies |
| NumPy vector operations | `Vec3`, `soa::` bulk operations | |
| `rnd.seed(seed)`, `rnd.random()` | `minimal::set_seed`, `minimal::random_uniform` | same default seed, 256; a different generator, so the numbers differ |
| `copy_bodies` | `minimal::copy_bodies` | |
| `random_direction` | `minimal::random_direction` | |
| `to_com_frame` | `minimal::to_com_frame` | |
| `euclidian_distance` | `minimal::euclidean_distance` | spelling corrected |
| `vel_init`, `get_perpendicular`, `vel_direction` | `minimal::` functions of the same names | |
| `init_plummer` | `procgen::init_plummer` | random numbers drawn in the same order |
| `virialise` | `procgen::virialise` | |
| `init_random_sphere` | | not yet ported |
| `t_abs` | `verify::t_abs` | $m\lvert v\rvert^2/2$ instead of $\lvert p\rvert^2/(2m)$; no division |
| `v_abs`, `h_abs`, `h_rel` | `verify::` functions of the same names | |
| `p_abs`, `p_scale`, `l_abs`, `l_scale`, `vec_rel` | `verify::` functions of the same names | plus the aliases `p_rel`, `l_rel` |
| `conservation_state` | `verify::conservation_state` | returns a struct `{H, P, L}` |
| `acc` | `physics::acc` | Newton's third law halves the work |

**Barnes-Hut**

| Notebook | C++ | Remark |
|---|---|---|
| `max_depth_bh` | `Nodes::max_depth` | same value, 32 |
| `class Node` (one object per cell) | `Nodes` (one object for all cells) | SoA instead of a linked tree |
| `octant_bh` | `Nodes::octant` | Eq. (16) |
| `subdivide_bh` | `Nodes::subdivide` | Eq. (17) |
| `insert_bh` | `Nodes::insert` | recursive, as in the notebook |
| `moments_bh` | `Nodes::moments` | Eq. (18); reverse sweep instead of recursion |
| `build_tree_bh` | `Nodes::build` | same bounding cube and margin |
| `contains_bh` | `Nodes::contains` | |
| `node_acc_bh` | `Nodes::acceleration` | Eq. (19); explicit stack instead of recursion |
| `acc_bh` | `physics::acc_bh` | |
| `acc_bh_theta` | `physics::acc_bh_theta` | |
| `force_error_bh` | `physics::force_error` | |
| `tree_stats_bh` | `Nodes::stats` | |

The following differences in the Barnes-Hut port are intentional and do not change the algorithm:

- The tree walk visits the children in the order 0 to 7, as the notebook does, but uses an explicit stack of fixed size, since the depth is bounded by `max_depth`.
- The opening test is evaluated as $s^2<\theta^2d^2$, which avoids a square root and is equivalent to $s<\theta d$ for $\theta\geq0$; a negative $\theta$ is therefore rejected.
- Softened distances use `s * sqrt(s)` instead of `pow(s, 1.5)`, as in `physics::acc`.
- The notebook's leaf stores a Python list of bodies; here the list is chained through `next_body`.

---

## 14. Validation of the Barnes-Hut port

The port was checked against the notebook itself by running the notebook's `build_tree_bh`, `acc_bh` and `tree_stats_bh` and the C++ code on the same input: 300 bodies with unequal masses, distributed as a Plummer sphere, which contains one pair at exactly the same position and one pair at a separation of $10^{-9}$, so that the depth limit and the leaf lists are exercised.

| Check | Result |
|---|---|
| C++ against notebook, $\theta=0$, largest relative difference in $\hat{a}_i$ | $2.0\times10^{-15}$ |
| C++ against notebook, $\theta=0.5$ | $1.5\times10^{-15}$ |
| C++ against notebook, $\theta=1$ | $1.5\times10^{-15}$ |
| Tree statistics (nodes, non-empty leaves, depth) | 1657, 298, 32 in both |
| $\theta=0$ against `physics::acc`, largest relative difference | $3.1\times10^{-15}$ |
| Total mass and centre of mass at the root against direct sums | agree to rounding |

The differences are at the level of floating-point rounding, as expected from a different order of additions. These checks were run with temporary test programs; no automated test suite exists in the repository yet, so a unit test built from the table above is a natural addition.

---

## 15. Accuracy and cost of the Barnes-Hut method

Measured on a Plummer sphere in Hénon units ($G=M=1$, $a=3\pi/16$, softening $0.05$), generated with `procgen::init_plummer`.

**Force error** ($N=4000$; tree of 15 529 nodes, depth 12):

| $\theta$ | median | 95th percentile | maximum |
|---|---|---|---|
| 0 | $1.8\times10^{-15}$ | $3.8\times10^{-15}$ | $1.3\times10^{-14}$ |
| 0.25 | $3.8\times10^{-4}$ | $9.1\times10^{-4}$ | $9.7\times10^{-3}$ |
| 0.5 | $2.1\times10^{-3}$ | $6.4\times10^{-3}$ | $6.5\times10^{-2}$ |
| 0.75 | $6.3\times10^{-3}$ | $1.8\times10^{-2}$ | $2.5\times10^{-1}$ |
| 1.0 | $1.5\times10^{-2}$ | $4.0\times10^{-2}$ | $2.5\times10^{-1}$ |

The error grows by a factor of 5 to 7 per doubling of $\theta$, which is somewhat steeper than the $\theta^2$ scaling of a single accepted cell (Eq. (20)), because the errors of many accepted cells accumulate. The maxima are attained by bodies near the centre of the cluster, where the contributions of the surrounding bodies nearly cancel, so that the small net force has a large relative error.

**Cost** of one force evaluation, including the construction of the tree (single thread, `-O2`, GCC 11.4 on a 4-core Linux virtual machine; one run per entry):

| $N$ | direct [s] | Barnes-Hut, $\theta=0.5$ [s] | Barnes-Hut, $\theta=1$ [s] | direct / BH(0.5) |
|---|---|---|---|---|
| 1 000 | 0.0010 | 0.0044 | 0.0014 | 0.2 |
| 4 000 | 0.0160 | 0.0341 | 0.0074 | 0.5 |
| 16 000 | 0.2475 | 0.1750 | 0.0335 | 1.4 |
| 64 000 | 4.0024 | 0.8282 | 0.1575 | 4.8 |

Direct summation grows by a factor of about 16 per quadrupling of $N$, whereas the tree code grows by a factor of 7.7, 5.1 and 4.7 in the three successive quadruplings of the full measurement (1 000 to 64 000), approaching the factor of about 4 of $\mathcal{O}(N\log N)$. The crossover at $\theta=0.5$ lies near $N\approx10^4$: below it the direct sum is faster, because the pair loop of `physics::acc` is very light and uses Newton's third law to halve its work. The absolute times depend on the machine and the build type, so only the scaling and the crossover should be taken from this table.

---

## 16. Known limitations

**General**

- **`verify::v_abs` is $\mathcal{O}(N^2)$.** Monitoring the energy of a large system is as expensive as one direct force evaluation, also when the forces come from the tree, so it should be done at intervals rather than at every step.
- **One shared random generator.** Reproducibility depends on the order of all draws in a program; `minimal::set_seed` should be called before each generator whose output is to be compared between programs. The C++ generator (`std::mt19937_64`) is not the one used by Python's `random`, so the same seed gives different, though statistically equivalent, initial conditions in the notebook and here.
- **Views are not owning.** A `soa::Vec3Span` taken from an `Object` or a `Vec3Array` dangles after `add`, `resize` or `reserve` on that container.

**Barnes-Hut**

- **Momentum is not conserved exactly.** The tree forces do not obey Newton's third law, so the total linear momentum and angular momentum drift, and the energy error is larger than for direct summation (see the notebook's Barnes-Hut explanation). The conserved quantities of `verify` always use the exact pair sums.
- **The tree is rebuilt, and its memory reallocated, at every `acc_bh` call.** A persistent `Nodes` object whose `build` is called repeatedly reuses its memory (Section 11.4); this is not yet exposed through `acc_bh`.
- **Monopole approximation only.** Each accepted cell contributes a point mass; higher multipole moments are not used.
- **Single-threaded.** `Nodes::acceleration` is `const` and independent for different bodies, so the force loop can be parallelised, but this has not been done.
- **Relative force errors are ill-defined where the exact force vanishes.** `force_error` returns infinity or NaN for a body whose direct acceleration is exactly zero.
- **Bodies closer together than about $2^{-32}$ of the root side share a leaf** and are never separated by further subdivision. They are summed exactly, but a dense clump of such bodies would make that leaf expensive. This is the behaviour of the notebook.

---

## References

- Aarseth, S. J., Hénon, M. and Wielen, R. (1974) A comparison of numerical methods for the study of star cluster dynamics. *Astronomy and Astrophysics*, 37, 183–187.
- Barnes, J. and Hut, P. (1986) A hierarchical $\mathcal{O}(N\log N)$ force-calculation algorithm. *Nature*, 324, 446–449.
- Plummer, H. C. (1911) On the problem of distribution in globular star clusters. *Monthly Notices of the Royal Astronomical Society*, 71, 460–470.
- Salmon, J. K. and Warren, M. S. (1994) Skeletons from the treecode closet. *Journal of Computational Physics*, 111(1), 136–155.
