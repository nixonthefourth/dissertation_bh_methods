//
// Created by Mykyta Khomiakov on 19/11/2026.
//

#include "physics/barnes_hut.h++"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace physics {
    std::size_t Nodes::octant(const std::size_t n, const Vec3& r) const {
        return (r.x >= centre.x[n] ? 1u : 0u)
             | (r.y >= centre.y[n] ? 2u : 0u)
             | (r.z >= centre.z[n] ? 4u : 0u);
    }

    bool Nodes::contains(const std::size_t n, const Vec3& r) const {
        return std::abs(r.x - centre.x[n]) <= half[n]
            && std::abs(r.y - centre.y[n]) <= half[n]
            && std::abs(r.z - centre.z[n]) <= half[n];
    }

    std::size_t Nodes::add(const Vec3& c, const double h) {
        const std::size_t n = size();

        // Vec3Array: grow the three axes together with resize(), then write the whole vector with set()
        centre.resize(n + 1);
        centre.view().set(n, c);
        com.resize(n + 1);      // zero-initialised: the sums of moments() start from zero

        half.push_back(h);
        mass.push_back(0.0);
        first_child.push_back(NONE);
        first_body.push_back(NONE);
        return n;
    }

    void Nodes::subdivide(const std::size_t n) {
        // Copy what is needed BEFORE add(), which may reallocate the arrays
        const Vec3 c = centre.view()[n];
        const double quarter = 0.5 * half[n];
        const std::size_t first = size();

        for (std::size_t k = 0; k < 8; ++k) {
            add({c.x + ((k & 1u) ? quarter : -quarter),
                 c.y + ((k & 2u) ? quarter : -quarter),
                 c.z + ((k & 4u) ? quarter : -quarter)}, quarter);
        }
        first_child[n] = first;
    }

    void Nodes::append(const std::size_t n, const std::size_t b) {
        next_body[b] = NONE;
        if (first_body[n] == NONE) {
            first_body[n] = b;
            return;
        }

        std::size_t last = first_body[n];
        while (next_body[last] != NONE) {
            last = next_body[last];
        }
        next_body[last] = b;
    }

    void Nodes::insert(const std::size_t n, const std::size_t b, const Object& bodies, const int depth) {
        const soa::ConstVec3Span r = bodies.displacement();

        if (is_leaf(n)) {
            // Empty leaf (or the depth limit is reached): simply store the body here
            if (first_body[n] == NONE || depth >= max_depth) {
                append(n, b);
                return;
            }

            // Occupied leaf: subdivide and push the resident bodies one level down
            std::size_t resident = first_body[n];
            first_body[n] = NONE;
            subdivide(n);
            while (resident != NONE) {
                const std::size_t following = next_body[resident];     // insert() rewrites next_body
                insert(child(n, octant(n, r[resident])), resident, bodies, depth + 1);
                resident = following;
            }
        }

        // Internal node: descend into the octant that contains body b
        insert(child(n, octant(n, r[b])), b, bodies, depth + 1);
    }

    void Nodes::moments(const Object& bodies) {
        const soa::ConstVec3Span r = bodies.displacement();
        const soa::Vec3Span R = com.view();

        for (std::size_t n = size(); n-- > 0;) {
            if (is_leaf(n)) {
                for (std::size_t j = first_body[n]; j != NONE; j = next_body[j]) {
                    mass[n] += bodies.mass[j];
                    R.add(n, bodies.mass[j] * r[j]);
                }
            } else {
                for (std::size_t k = 0; k < 8; ++k) {
                    const std::size_t c = child(n, k);
                    mass[n] += mass[c];
                    R.add(n, mass[c] * R[c]);           // R[c] is already a centre of mass, not a sum
                }
            }

            if (mass[n] > 0.0) {
                R.set(n, R[n] / mass[n]);               // the numerator above is sum m r
            }
        }
    }

    void Nodes::build(const Object& bodies) {
        // Discard the previous tree (the allocated memory is kept for the next one)
        centre.resize(0);
        com.resize(0);
        half.clear();
        mass.clear();
        first_child.clear();
        first_body.clear();

        const std::size_t N = bodies.size();
        next_body.assign(N, NONE);
        if (N == 0) {
            return;
        }

        // Smallest axis-aligned box around all bodies, made into a cube
        const soa::ConstVec3Span r = bodies.displacement();
        Vec3 lower = r[0], upper = r[0];
        for (std::size_t i = 1; i < N; ++i) {
            const Vec3 v = r[i];
            lower = {std::min(lower.x, v.x), std::min(lower.y, v.y), std::min(lower.z, v.z)};
            upper = {std::max(upper.x, v.x), std::max(upper.y, v.y), std::max(upper.z, v.z)};
        }
        const Vec3 extent = upper - lower;
        const double h = 0.5 * std::max({extent.x, extent.y, extent.z}) * (1.0 + 1e-6) + 1e-12;  // +1e-12 guards a single body

        add(0.5 * (lower + upper), h);
        for (std::size_t i = 0; i < N; ++i) {
            insert(0, i, bodies, 0);
        }
        moments(bodies);
    }

    Vec3 Nodes::acceleration(const std::size_t i, const Object& bodies, const double G, const double e,
                             const double theta) const {
        const soa::ConstVec3Span r = bodies.displacement();
        const soa::ConstVec3Span R = com.view();
        const Vec3 r_i = r[i];
        const double e2 = e * e;

        // Depth-first walk with an explicit stack. Every opened node replaces itself by its 8 children,
        // so at most 7 nodes per level stay waiting, plus the 8 children of the node being opened.
        // The children are pushed in reverse order, so they are visited in the order 0..7.
        std::array<std::size_t, 7 * max_depth + 8> stack;
        std::size_t top = 0;
        stack[top++] = 0;

        Vec3 a{};
        while (top > 0) {
            const std::size_t n = stack[--top];
            if (mass[n] == 0.0) {
                continue;                               // empty cube: no contribution
            }

            // Leaf: exact summation over its bodies, skipping the self-interaction
            if (is_leaf(n)) {
                for (std::size_t j = first_body[n]; j != NONE; j = next_body[j]) {
                    if (j == i) {
                        continue;
                    }
                    const Vec3 r_ij = r_i - r[j];
                    const double s = norm_sq(r_ij) + e2;
                    a += (-G * bodies.mass[j] / (s * std::sqrt(s))) * r_ij;
                }
                continue;
            }

            // Internal node: accept it as one pseudo-body or open it. s < theta d is the criterion
            // s / d < theta written without a division, and the containment test is only evaluated for
            // nodes that already pass it.
            const Vec3 r_iN = r_i - R[n];               // from the centre of mass of the node to body i
            const double d2 = norm_sq(r_iN);
            const double side = 2.0 * half[n];
            if (side * side < theta * theta * d2 && !contains(n, r_i)) {
                const double s = d2 + e2;
                a += (-G * mass[n] / (s * std::sqrt(s))) * r_iN;
                continue;
            }

            for (std::size_t k = 8; k-- > 0;) {
                stack[top++] = child(n, k);
            }
        }
        return a;
    }

    TreeStats Nodes::stats() const {
        TreeStats s{size(), 0, 0};
        std::vector<int> depth(size(), 0);

        // A parent always comes before its children, so depth[n] is final when node n is reached
        for (std::size_t n = 0; n < size(); ++n) {
            s.depth = std::max(s.depth, depth[n]);
            if (is_leaf(n)) {
                s.leaves += (first_body[n] != NONE);
            } else {
                for (std::size_t k = 0; k < 8; ++k) {
                    depth[child(n, k)] = depth[n] + 1;
                }
            }
        }
        return s;
    }

    void acc_bh(const double G, const Object& bodies, const double e, const double theta, const soa::Vec3Span out) {
        const std::size_t N = bodies.size();
        if (out.size() != N) {
            throw std::invalid_argument("acc_bh: output view must have one entry per body");
        }
        if (!(theta >= 0.0)) {
            throw std::invalid_argument("acc_bh: the opening angle theta must not be negative");
        }

        Nodes tree;
        tree.build(bodies);
        for (std::size_t i = 0; i < N; ++i) {
            out.set(i, tree.acceleration(i, bodies, G, e, theta));
        }
    }

    void acc_bh(const double G, Object& bodies, const double e, const double theta) {
        acc_bh(G, std::as_const(bodies), e, theta, bodies.acceleration());
    }

    std::vector<double> force_error(const soa::ConstVec3Span a_tree, const soa::ConstVec3Span a_direct) {
        if (a_tree.size() != a_direct.size()) {
            throw std::invalid_argument("force_error: both views must have one entry per body");
        }

        std::vector<double> error(a_tree.size());
        for (std::size_t i = 0; i < error.size(); ++i) {
            error[i] = norm(a_tree[i] - a_direct[i]) / norm(a_direct[i]);
        }
        return error;
    }
}
