#pragma once
#include <algorithm>
#include <stdexcept>
#include <vector>

#include "../geom/plane.hpp"
#include "../geom/triangle.hpp"
#include "../geom/project.hpp"
#include "../geom/basis.hpp"
#include "../geom/plane2d.hpp"
#include "../geom/clip2d.hpp"

namespace em {

struct OverlapMetrics {
    double area_a{0.0};
    double area_b{0.0};
    double area_intersection{0.0};
    double iou{0.0};          // intersection / union
    double overlap_min{0.0};  // intersection / min(area_a, area_b)
};

inline OverlapMetrics triangle_overlap_on_plane(
    const geom::Triangle& tri_on_interface,
    const geom::Plane& bottom_plane,
    const geom::Vec3& dir_a,
    const geom::Vec3& dir_b,
    double eps = 1e-12
) {
    geom::Triangle A3, B3;
    if (!geom::project_triangle_along(tri_on_interface, dir_a, bottom_plane, A3, eps)) {
        throw std::runtime_error("Projection failed for dir_a (ray parallel to bottom plane)");
    }
    if (!geom::project_triangle_along(tri_on_interface, dir_b, bottom_plane, B3, eps)) {
        throw std::runtime_error("Projection failed for dir_b (ray parallel to bottom plane)");
    }

    // 2D coordinate system on the bottom plane
    const geom::OrthoBasis basis = geom::make_plane_basis(bottom_plane.n, eps);
    const geom::Vec3 origin = bottom_plane.p0;

    std::vector<geom::Vec2> A2 = geom::triangle_to2D(A3, origin, basis);
    std::vector<geom::Vec2> B2 = geom::triangle_to2D(B3, origin, basis);

    const double areaA = geom::area_abs(A2);
    const double areaB = geom::area_abs(B2);

    // intersection polygon via Sutherland–Hodgman
    const std::vector<geom::Vec2> I = geom::convex_clip(A2, B2, eps);
    const double areaI = geom::area_abs(I);

    OverlapMetrics m;
    m.area_a = areaA;
    m.area_b = areaB;
    m.area_intersection = areaI;

    const double unionA = areaA + areaB - areaI;
    if (unionA > eps) m.iou = areaI / unionA;

    const double minA = std::min(areaA, areaB);
    if (minA > eps) m.overlap_min = areaI / minA;

    return m;
}

} // namespace em
