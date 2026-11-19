//
// Created by Mykyta Khomiakov on 19/11/2026.
//

#ifndef MAIN_BARNES_HUT_H
#define MAIN_BARNES_HUT_H

#include <cstddef>
#include <limits>
#include <vector>

#include "helpers/Vec3.h++"
#include "helpers/object.h++"
#include "helpers/soa.h++"

/// @brief Force evaluation
namespace physics {
    /// @brief Size and shape of an octree (port of the notebook's tree_stats_bh)
    /// @details The depth shows how finely the tree has to resolve the closest pair of bodies: a tight
    ///          binary forces many subdivisions before its two members end up in separate leaves.
    struct TreeStats {
        std::size_t nodes;      // total number of nodes, empty ones included
        std::size_t leaves;     // number of leaves that hold at least one body
        int depth;              // depth of the deepest node (the root has depth 0)
    };

    /// @brief Barnes-Hut octree (Barnes & Hut 1986), with all nodes stored as a Structure of Arrays (SoA)
    /// @details Node n is the n-th entry of every array, exactly like body i is the i-th entry of every
    ///          array of Object. Nodes refer to each other by index, never by pointer, so the arrays can
    ///          grow while the tree is being built. The root is node 0.
    ///              - centre[n]       geometric centre of the cubic cell
    ///              - half[n]         half of the side length of the cell, so that the side is s = 2 half
    ///              - mass[n]         total mass M_c of all bodies inside the cell, Eq. (18)
    ///              - com[n]          centre of mass R_c of those bodies, Eq. (18)
    ///              - first_child[n]  index of child 0 of an internal node, NONE for a leaf
    ///              - first_body[n]   index of the first body held by a leaf, NONE for an empty leaf
    ///              - next_body[b]    index of the body that follows body b in the same leaf, NONE at the end
    ///                                (one entry per BODY, not per node)
    ///
    ///          The 8 children of a node are allocated together, so child o (0..7) of node n is the node
    ///          first_child[n] + o. A child therefore always has a larger index than its parent, which
    ///          lets the bottom-up passes run as a simple reverse sweep over the arrays.
    ///
    ///          A leaf normally holds a single body. It holds several only at max_depth, where a cell can
    ///          no longer be subdivided; the bodies of such a leaf are chained through next_body.
    ///
    ///          The tree is a snapshot of one configuration: it has to be rebuilt with build() whenever
    ///          the bodies have moved (acc_bh does this at every call).
    class Nodes {
        public:
            /// @brief Marks "no child" / "no body"
            static constexpr std::size_t NONE = std::numeric_limits<std::size_t>::max();

            /// @brief Maximum depth of the tree
            /// @details Two bodies at (almost) the same position can never be separated by subdividing
            ///          the cube, so without a limit the insertion would recurse forever. At depth 32 a
            ///          cell is 2^32 times smaller than the root, far below any physical separation here,
            ///          so the limit only acts as a safety net.
            static constexpr int max_depth = 32;

            soa::Vec3Array centre;                  // geometric centres of the cells
            soa::Vec3Array com;                     // centres of mass of the cells
            std::vector<double> half;               // half side lengths
            std::vector<double> mass;               // total masses
            std::vector<std::size_t> first_child;   // index of child 0, NONE for a leaf
            std::vector<std::size_t> first_body;    // first body of a leaf, NONE if empty
            std::vector<std::size_t> next_body;     // next body in the same leaf, one entry per body

            /// @brief Number of nodes in the tree
            [[nodiscard]] std::size_t size() const { return mass.size(); }

            /// @brief True if node n has no children
            [[nodiscard]] bool is_leaf(std::size_t n) const { return first_child[n] == NONE; }

            /// @brief Index of child o (0..7) of the internal node n
            [[nodiscard]] std::size_t child(std::size_t n, std::size_t o) const { return first_child[n] + o; }

            /// @brief Octant of node n that contains the displacement r (Eq. (16))
            /// @param n Index of the node
            /// @param r Displacement vector
            /// @details Each bit of the result answers one question: bit 0 is set if x >= c_x, bit 1 if
            ///          y >= c_y and bit 2 if z >= c_z, so that k = [x >= c_x] + 2 [y >= c_y] + 4 [z >= c_z].
            /// @return Returns k in 0..7
            [[nodiscard]] std::size_t octant(std::size_t n, const Vec3& r) const;

            /// @brief True if the displacement r lies inside the cube of node n (boundary included)
            [[nodiscard]] bool contains(std::size_t n, const Vec3& r) const;

            /// @brief Builds the complete octree for the current positions of all bodies
            /// @param bodies The system of stellar bodies (displacements and masses are read)
            /// @details Any previous tree is discarded (its memory is reused). The root is the smallest
            ///          cube that encloses every body, enlarged by a tiny margin so that no body lies
            ///          exactly on its outer boundary. The bodies are then inserted one by one, and the
            ///          mass and centre of mass of every node are computed (Eq. (18)).
            void build(const Object& bodies);

            /// @brief Gravitational acceleration of body i, summed over the tree (Eq. (19))
            /// @param i Index of the body that feels the force
            /// @param bodies The system of stellar bodies the tree was built from
            /// @param G Constant of gravitational attraction
            /// @param e Plummer softening factor (the same as in acc)
            /// @param theta Opening angle: theta = 0 opens every node and reproduces direct summation,
            ///        larger values are faster but less accurate (0.5 to 1.0 is the usual range)
            /// @details A node of side s = 2 half seen from the distance d = |r_i - R_c| is accepted as a
            ///          single softened point mass at its centre of mass if s < theta d, and is otherwise
            ///          opened and its children are visited. A node that contains body i itself is always
            ///          opened, whatever the criterion says (Salmon & Warren 1994). Leaves are summed
            ///          exactly, which is Eq. (4) restricted to the leaf, skipping the body itself.
            /// @return Returns a_i
            [[nodiscard]] Vec3 acceleration(std::size_t i, const Object& bodies, double G, double e, double theta) const;

            /// @brief Number of nodes, number of non-empty leaves and maximum depth of the tree
            [[nodiscard]] TreeStats stats() const;

        private:
            /// @brief Appends one empty leaf cell
            /// @param c Centre of the cell
            /// @param h Half side length of the cell
            /// @return Index of the new node
            std::size_t add(const Vec3& c, double h);

            /// @brief Splits the leaf n into 8 empty children (Eq. (17))
            /// @details Child k has half the side length of its parent and is shifted from the parent
            ///          centre by a quarter of the parent side along every axis, with the sign of each
            ///          shift taken from the corresponding bit of k (the same rule as octant()).
            void subdivide(std::size_t n);

            /// @brief Appends body b to the end of the list of bodies of the leaf n
            void append(std::size_t n, std::size_t b);

            /// @brief Inserts body b into the subtree that starts at node n
            /// @param n Current node of the descent
            /// @param b Index of the body
            /// @param bodies The system of stellar bodies
            /// @param depth Depth of node n (the root has depth 0)
            /// @details A leaf may hold a single body. When a second body arrives, the leaf is
            ///          subdivided, the body already living there is pushed down into its octant, and the
            ///          insertion continues into the octant of the new body. This repeats until every body
            ///          sits alone in a leaf, or until max_depth is reached, where bodies simply share a leaf.
            void insert(std::size_t n, std::size_t b, const Object& bodies, int depth);

            /// @brief Computes the total mass and the centre of mass of every node (Eq. (18))
            /// @param bodies The system of stellar bodies
            /// @details Leaf: M = sum of the masses of its bodies, R = their mass-weighted mean position.
            ///          Internal node: M = sum_k M_k, R = sum_k M_k R_k / M over its 8 children. Because a
            ///          child has a larger index than its parent, a sweep from the last node to the first
            ///          is a post-order traversal.
            void moments(const Object& bodies);
    };

    /// @brief Computes the gravitational acceleration of every body with the Barnes-Hut approximation
    /// @param G Constant of gravitational attraction
    /// @param bodies The system of stellar bodies (displacements and masses are read)
    /// @param e Softening factor to prevent division by zero if bodies get too close
    /// @param theta Opening angle (not negative), see Nodes::acceleration
    /// @param out Output view of size N, the counterpart of acc: a drop-in replacement with the same
    ///        output, O(N log N) instead of O(N^2) operations
    /// @details The tree is rebuilt at every call, since the bodies move between force evaluations.
    ///          Unlike acc, the pair forces are not applied equally and oppositely, so the total momentum
    ///          is no longer conserved exactly.
    void acc_bh(double G, const Object& bodies, double e, double theta, soa::Vec3Span out);

    /// @brief Computes the Barnes-Hut accelerations and stores them in the system's ax, ay, az
    /// @param G Constant of gravitational attraction
    /// @param bodies The system of stellar bodies
    /// @param e Softening factor
    /// @param theta Opening angle (default 0.5, the usual accuracy setting)
    void acc_bh(double G, Object& bodies, double e, double theta = 0.5);

    /// @brief Barnes-Hut force routine with a fixed opening angle
    /// @param theta Opening angle
    /// @details acc_bh takes theta as an additional argument, whereas an integrator calls the force routine
    ///          as f(G, bodies, e). This returns such a callable with theta captured by value, so it is
    ///          safe to create inside a loop over several opening angles.
    /// @return Returns a callable with the signature void(double G, Object& bodies, double e)
    inline auto acc_bh_theta(double theta) {
        return [theta](double G, Object& bodies, double e) { acc_bh(G, bodies, e, theta); };
    }

    /// @brief Relative force error of every body, |a_i(BH) - a_i(direct)| / |a_i(direct)|
    /// @param a_tree Accelerations from acc_bh
    /// @param a_direct Accelerations from acc for the same configuration
    /// @details Used to validate the tree (theta = 0 must give zero up to rounding) and to measure its
    ///          accuracy. A body whose direct acceleration is exactly zero gives inf or NaN.
    /// @return Returns one relative error per body
    std::vector<double> force_error(soa::ConstVec3Span a_tree, soa::ConstVec3Span a_direct);
}

#endif //MAIN_BARNES_HUT_H
