#pragma once
#include "vec2.hpp"
#include <vector>
#include <algorithm>

namespace geom {

inline double signed_area(const std::vector<Vec2>& poly) {
    double a = 0.0;
    const size_t n = poly.size();
    if (n < 3) return 0.0;
    for (size_t i = 0; i < n; ++i) {
        const Vec2& p = poly[i];
        const Vec2& q = poly[(i + 1) % n];
        a += p.x * q.y - q.x * p.y;
    }
    return 0.5 * a;
}

inline double area_abs(const std::vector<Vec2>& poly) {
    return std::abs(signed_area(poly));
}

inline void make_ccw(std::vector<Vec2>& poly) {
    if (signed_area(poly) < 0.0) {
        std::reverse(poly.begin(), poly.end());
    }
}

// Inside test for half-plane defined by directed edge A->B of a CCW clip polygon.
// P is inside if on the left side of the edge.
inline bool inside_halfplane(const Vec2& A, const Vec2& B, const Vec2& P, double eps = 1e-12) {
    return cross_z(B - A, P - A) >= -eps;
}

// Intersection between segment P1->P2 and infinite line A->B.
inline Vec2 line_intersection(const Vec2& P1, const Vec2& P2, const Vec2& A, const Vec2& B, double eps = 1e-12) {
    const Vec2 r = P2 - P1;
    const Vec2 s = B - A;

    const double denom = cross_z(r, s);
    if (std::abs(denom) < eps) {
        return (P1 + P2) * 0.5; // fallback for near-parallel
    }

    const double t = cross_z(A - P1, s) / denom;
    return P1 + r * t;
}

// Sutherland–Hodgman clipping of subject polygon by a convex clip polygon.
// For our use: triangles -> convex.
inline std::vector<Vec2> convex_clip(const std::vector<Vec2>& subject_in,
                                     const std::vector<Vec2>& clip_in,
                                     double eps = 1e-12) {
    if (subject_in.size() < 3 || clip_in.size() < 3) return {};

    std::vector<Vec2> subject = subject_in;
    std::vector<Vec2> clip = clip_in;

    make_ccw(subject);
    make_ccw(clip);

    std::vector<Vec2> output = subject;

    for (size_t i = 0; i < clip.size(); ++i) {
        const Vec2 A = clip[i];
        const Vec2 B = clip[(i + 1) % clip.size()];

        std::vector<Vec2> input = output;
        output.clear();
        if (input.empty()) break;

        Vec2 S = input.back();
        for (const Vec2& E : input) {
            const bool Ein = inside_halfplane(A, B, E, eps);
            const bool Sin = inside_halfplane(A, B, S, eps);

            if (Ein) {
                if (!Sin) output.push_back(line_intersection(S, E, A, B, eps));
                output.push_back(E);
            } else if (Sin) {
                output.push_back(line_intersection(S, E, A, B, eps));
            }
            S = E;
        }
    }

    auto near_eq = [&](const Vec2& p, const Vec2& q) {
        const double dx = p.x - q.x;
        const double dy = p.y - q.y;
        return (dx*dx + dy*dy) < (eps*eps);
    };

    if (!output.empty()) {
        std::vector<Vec2> cleaned;
        cleaned.reserve(output.size());
        for (const auto& p : output) {
            if (cleaned.empty() || !near_eq(cleaned.back(), p)) cleaned.push_back(p);
        }
        if (cleaned.size() >= 2 && near_eq(cleaned.front(), cleaned.back())) cleaned.pop_back();
        output.swap(cleaned);
    }

    if (output.size() < 3) return {};
    make_ccw(output);
    return output;
}

} // namespace geom
