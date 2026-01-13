#pragma once
#include <stdexcept>
#include <cmath>
#include "../geom/triangle.hpp"
#include "../geom/plane.hpp"
#include "../geom/vec3.hpp"
#include "overlap_analyzer.hpp"

namespace em {

// Build a bottom plane at distance h from the triangle centroid along the given normal (unit).
// This corresponds to defining thickness along bottom_normal.
inline geom::Plane bottom_plane_from_triangle_centroid(const geom::Triangle& tri,
                                                      const geom::Vec3& bottom_normal_unit,
                                                      double h) {
    const geom::Vec3 c = geom::centroid(tri);
    return geom::Plane::fromPointNormal(c + bottom_normal_unit * h, bottom_normal_unit);
}

// Find the maximum thickness h in [h_min, h_max] such that overlap metric >= threshold.
// - If use_iou=true: metric = IoU = A_int / A_union
// - else: metric = A_int / min(Aa, Ab)
// Assumes (empirically) metric decreases with h for fixed directions. If not monotonic,
// bisection may find a local boundary. For a thesis prototype it's typically adequate.
inline double find_max_thickness_for_overlap(
    const geom::Triangle& tri_interface,
    const geom::Vec3& bottom_normal_unit,
    const geom::Vec3& dir_a,
    const geom::Vec3& dir_b,
    double threshold,
    double h_min,
    double h_max,
    bool use_iou = true,
    int max_iter = 60,
    double eps = 1e-12
) {
    if (threshold < 0.0 || threshold > 1.0) throw std::runtime_error("threshold must be in [0,1]");
    if (h_min < 0.0 || h_max <= h_min) throw std::runtime_error("invalid h range");

    auto metric_at = [&](double h) {
        const geom::Plane bottom = bottom_plane_from_triangle_centroid(tri_interface, bottom_normal_unit, h);
        const OverlapMetrics m = triangle_overlap_on_plane(tri_interface, bottom, dir_a, dir_b, eps);
        return use_iou ? m.iou : m.overlap_min;
    };

    const double m_lo = metric_at(h_min);
    if (m_lo < threshold) {
        // Already failing at minimal thickness -> no feasible thickness in the interval.
        return h_min;
    }

    const double m_hi = metric_at(h_max);
    if (m_hi >= threshold) {
        // Still ok at the maximum bound -> feasible up to h_max.
        return h_max;
    }

    double lo = h_min;
    double hi = h_max;
    for (int i = 0; i < max_iter; ++i) {
        const double mid = 0.5 * (lo + hi);
        const double m_mid = metric_at(mid);
        if (m_mid >= threshold) lo = mid;
        else hi = mid;
    }
    return lo;
}

} // namespace em
