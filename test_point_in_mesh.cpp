// Build: g++ -std=c++17 -O2 -Wall -Wextra test_point_in_mesh.cpp -o test && ./test
#include "point_in_mesh.hpp"

#include <cassert>
#include <chrono>
#include <cstdio>
#include <functional>
#include <map>

using Tri = TriangleMesh::Tri;

TriangleMesh makeCube(bool dropTopFace = false) {
    std::vector<Vec3> v = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                           {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    std::vector<Tri> t = {{0, 2, 1}, {0, 3, 2}, {0, 1, 5}, {0, 5, 4}, {1, 2, 6}, {1, 6, 5},
                          {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7}};
    if (!dropTopFace) {
        t.push_back({4, 5, 6});
        t.push_back({4, 6, 7});
    }
    return TriangleMesh(v, t);
}

// Unit icosphere, outward-facing
TriangleMesh makeSphere(int subdivisions) {
    const double p = (1.0 + std::sqrt(5.0)) / 2.0;
    std::vector<Vec3> v = {{-1, p, 0}, {1, p, 0}, {-1, -p, 0}, {1, -p, 0}, {0, -1, p}, {0, 1, p},
                           {0, -1, -p}, {0, 1, -p}, {p, 0, -1}, {p, 0, 1}, {-p, 0, -1}, {-p, 0, 1}};
    for (auto& x : v) x = x * (1.0 / norm(x));
    std::vector<Tri> t = {{0, 11, 5}, {0, 5, 1},  {0, 1, 7},   {0, 7, 10}, {0, 10, 11},
                          {1, 5, 9},  {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
                          {3, 9, 4},  {3, 4, 2},  {3, 2, 6},   {3, 6, 8},  {3, 8, 9},
                          {4, 9, 5},  {2, 4, 11}, {6, 2, 10},  {8, 6, 7},  {9, 8, 1}};
    for (int s = 0; s < subdivisions; ++s) {
        std::map<std::pair<int, int>, int> mids;
        auto mid = [&](int a, int b) {
            auto key = std::minmax(a, b);
            auto it = mids.find(key);
            if (it != mids.end()) return it->second;
            Vec3 m = (v[a] + v[b]) * 0.5;
            v.push_back(m * (1.0 / norm(m)));
            return mids[key] = static_cast<int>(v.size()) - 1;
        };
        std::vector<Tri> next;
        for (auto& f : t) {
            int ab = mid(f[0], f[1]), bc = mid(f[1], f[2]), ca = mid(f[2], f[0]);
            next.push_back({f[0], ab, ca});
            next.push_back({f[1], bc, ab});
            next.push_back({f[2], ca, bc});
            next.push_back({ab, bc, ca});
        }
        t.swap(next);
    }
    return TriangleMesh(v, t);
}

// Torus in the xy-plane: non-convex, genus 1
TriangleMesh makeTorus(double R, double r, int nu, int nv) {
    std::vector<Vec3> v;
    std::vector<Tri> t;
    for (int i = 0; i < nu; ++i)
        for (int j = 0; j < nv; ++j) {
            double u = 2 * M_PI * i / nu, w = 2 * M_PI * j / nv;
            v.push_back({(R + r * std::cos(w)) * std::cos(u), (R + r * std::cos(w)) * std::sin(u),
                         r * std::sin(w)});
        }
    auto id = [&](int i, int j) { return ((i % nu) * nv) + (j % nv); };
    for (int i = 0; i < nu; ++i)
        for (int j = 0; j < nv; ++j) {
            t.push_back({id(i, j), id(i + 1, j), id(i + 1, j + 1)});
            t.push_back({id(i, j), id(i + 1, j + 1), id(i, j + 1)});
        }
    return TriangleMesh(v, t);
}

int main() {
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> U(-1.3, 1.3);

    // 1. Cube: hand-picked points, including ones aligned with edges and vertices
    TriangleMesh cube = makeCube();
    for (Vec3 p : {Vec3{0.5, 0.5, 0.5}, Vec3{0.1, 0.9, 0.2}, Vec3{0.5, 0.5, 0.0001}}) {
        assert(cube.containsWinding(p) && cube.containsRayCast(p));
    }
    for (Vec3 p : {Vec3{1.5, 0.5, 0.5}, Vec3{-0.2, -0.2, -0.2}, Vec3{2, 2, 2}, Vec3{0.5, 0.5, 1.2}}) {
        assert(!cube.containsWinding(p) && !cube.containsRayCast(p));
    }
    std::printf("cube: hand-picked cases passed\n");

    // 2. Sphere and torus vs analytic ground truth (skip a thin band where the
    //    faceted mesh and the true surface legitimately disagree)
    struct Case {
        const char* name;
        TriangleMesh mesh;
        std::function<double(const Vec3&)> sdf;  // < 0 inside
        double band;
    };
    std::vector<Case> cases;
    cases.push_back({"sphere", makeSphere(5), [](const Vec3& p) { return norm(p) - 1.0; }, 0.01});
    cases.push_back({"torus", makeTorus(0.8, 0.3, 96, 48),
                     [](const Vec3& p) {
                         double q = std::sqrt(p.x * p.x + p.y * p.y) - 0.8;
                         return std::sqrt(q * q + p.z * p.z) - 0.3;
                     },
                     0.01});

    const int N = 20000;
    for (auto& c : cases) {
        std::vector<Vec3> pts;
        while (static_cast<int>(pts.size()) < N) {
            Vec3 p{U(rng), U(rng), U(rng)};
            if (std::abs(c.sdf(p)) > c.band) pts.push_back(p);
        }

        int wOk = 0, rOk = 0;
        auto t0 = std::chrono::steady_clock::now();
        for (auto& p : pts) wOk += c.mesh.containsWinding(p) == (c.sdf(p) < 0);
        auto t1 = std::chrono::steady_clock::now();
        for (auto& p : pts) rOk += c.mesh.containsRayCast(p) == (c.sdf(p) < 0);
        auto t2 = std::chrono::steady_clock::now();

        double wUs = std::chrono::duration<double, std::micro>(t1 - t0).count() / N;
        double rUs = std::chrono::duration<double, std::micro>(t2 - t1).count() / N;
        std::printf("%-7s tris=%-6zu winding %d/%d (%.1f us/query) | BVH ray cast %d/%d (%.2f us/query) | speedup %.0fx\n",
                    c.name, c.mesh.triangleCount(), wOk, N, wUs, rOk, N, rUs, wUs / rUs);
        assert(wOk == N && rOk == N);
    }

    // 3. Open mesh (cube missing its top): winding number still classifies the
    //    interior; single-ray parity breaks whenever the ray escapes through the hole
    TriangleMesh open = makeCube(true);
    Vec3 center{0.5, 0.5, 0.5};
    std::printf("open cube: winding(center) = %.4f -> inside=%d\n", open.windingNumber(center),
                open.containsWinding(center));
    assert(open.containsWinding(center));
    assert(!open.containsWinding({0.5, 0.5, 1.5}));

    std::printf("all tests passed\n");
}
