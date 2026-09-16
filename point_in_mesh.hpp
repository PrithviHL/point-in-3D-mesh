#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <random>
#include <vector>

struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    double operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
};

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

struct AABB {
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    void expand(const Vec3& p) {
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
        hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    }
    // Slab test for a ray starting at o with precomputed 1/dir
    bool hitBy(const Vec3& o, const Vec3& inv) const {
        double t0 = 0, t1 = 1e300;
        for (int a = 0; a < 3; ++a) {
            double tn = (lo[a] - o[a]) * inv[a], tf = (hi[a] - o[a]) * inv[a];
            if (tn > tf) std::swap(tn, tf);
            t0 = std::max(t0, tn);
            t1 = std::min(t1, tf);
            if (t0 > t1) return false;
        }
        return true;
    }
};

class TriangleMesh {
public:
    using Tri = std::array<int, 3>;

    TriangleMesh(std::vector<Vec3> verts, std::vector<Tri> tris)
        : v_(std::move(verts)), t_(std::move(tris)) {
        order_.resize(t_.size());
        for (size_t i = 0; i < t_.size(); ++i) order_[i] = static_cast<int>(i);
        if (!t_.empty()) root_ = build(0, static_cast<int>(t_.size()));
    }

    // Generalized winding number (Jacobson et al. 2013): sum of signed solid angles / 4pi.
    // ~1 inside, ~0 outside for a closed, consistently oriented mesh; degrades gracefully
    // when the mesh has holes, which is where plain ray parity breaks.
    double windingNumber(const Vec3& p) const {
        double total = 0;
        for (const Tri& t : t_) {
            Vec3 a = v_[t[0]] - p, b = v_[t[1]] - p, c = v_[t[2]] - p;
            double la = norm(a), lb = norm(b), lc = norm(c);
            // Van Oosterom-Strackee solid angle formula
            double num = dot(a, cross(b, c));
            double den = la * lb * lc + dot(a, b) * lc + dot(b, c) * la + dot(c, a) * lb;
            total += 2.0 * std::atan2(num, den);
        }
        return total / (4.0 * M_PI);
    }

    bool containsWinding(const Vec3& p) const { return std::abs(windingNumber(p)) > 0.5; }

    // Parity ray casting over a BVH. A ray can graze an edge or vertex and be counted
    // twice or not at all, so cast three random directions and take the majority vote.
    bool containsRayCast(const Vec3& p, uint32_t seed = 12345) const {
        std::mt19937 rng(seed);
        std::normal_distribution<double> g(0.0, 1.0);
        int insideVotes = 0;
        for (int r = 0; r < 3; ++r) {
            Vec3 d{g(rng), g(rng), g(rng)};
            if (countHits(p, d) % 2 == 1) ++insideVotes;
        }
        return insideVotes >= 2;
    }

    size_t triangleCount() const { return t_.size(); }

private:
    struct Node {
        AABB box;
        int start = 0, end = 0;  // leaf range into order_
        std::unique_ptr<Node> left, right;
    };

    std::vector<Vec3> v_;
    std::vector<Tri> t_;
    std::vector<int> order_;
    std::unique_ptr<Node> root_;

    Vec3 centroid(int i) const {
        const Tri& t = t_[i];
        return (v_[t[0]] + v_[t[1]] + v_[t[2]]) * (1.0 / 3.0);
    }

    std::unique_ptr<Node> build(int start, int end) {
        auto node = std::make_unique<Node>();
        node->start = start;
        node->end = end;
        AABB cbox;
        for (int i = start; i < end; ++i) {
            for (int k : t_[order_[i]]) node->box.expand(v_[k]);
            cbox.expand(centroid(order_[i]));
        }
        if (end - start <= 4) return node;

        // Median split along the widest centroid axis
        Vec3 ext = cbox.hi - cbox.lo;
        int axis = (ext.x > ext.y && ext.x > ext.z) ? 0 : (ext.y > ext.z ? 1 : 2);
        int mid = (start + end) / 2;
        std::nth_element(order_.begin() + start, order_.begin() + mid, order_.begin() + end,
                         [&](int a, int b) { return centroid(a)[axis] < centroid(b)[axis]; });
        node->left = build(start, mid);
        node->right = build(mid, end);
        return node;
    }

    // Moller-Trumbore, counting only hits strictly in front of the origin
    bool rayHitsTriangle(const Vec3& o, const Vec3& d, const Tri& t) const {
        const double eps = 1e-12;
        Vec3 e1 = v_[t[1]] - v_[t[0]], e2 = v_[t[2]] - v_[t[0]];
        Vec3 pv = cross(d, e2);
        double det = dot(e1, pv);
        if (std::abs(det) < eps) return false;
        double inv = 1.0 / det;
        Vec3 tv = o - v_[t[0]];
        double u = dot(tv, pv) * inv;
        if (u < 0.0 || u > 1.0) return false;
        Vec3 qv = cross(tv, e1);
        double w = dot(d, qv) * inv;
        if (w < 0.0 || u + w > 1.0) return false;
        return dot(e2, qv) * inv > eps;
    }

    int countHits(const Vec3& o, const Vec3& d) const {
        Vec3 inv{1.0 / d.x, 1.0 / d.y, 1.0 / d.z};
        int hits = 0;
        std::vector<const Node*> stack{root_.get()};
        while (!stack.empty()) {
            const Node* n = stack.back();
            stack.pop_back();
            if (!n || !n->box.hitBy(o, inv)) continue;
            if (!n->left) {
                for (int i = n->start; i < n->end; ++i)
                    if (rayHitsTriangle(o, d, t_[order_[i]])) ++hits;
            } else {
                stack.push_back(n->left.get());
                stack.push_back(n->right.get());
            }
        }
        return hits;
    }
};
