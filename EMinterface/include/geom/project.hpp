#pragma once
#include "triangle.hpp"
#include "plane.hpp"

namespace geom {

// Projects triangle vertices along direction d onto plane pl.
// Returns false if projection fails (ray parallel to plane) for any vertex.
inline bool project_triangle_along(const Triangle& tri, const Vec3& d, const Plane& pl, Triangle& out, double eps = 1e-12) {
    double t0, t1, t2;
    if (!intersect_ray_plane(tri.v0, d, pl, t0, eps)) return false;
    if (!intersect_ray_plane(tri.v1, d, pl, t1, eps)) return false;
    if (!intersect_ray_plane(tri.v2, d, pl, t2, eps)) return false;

    out.v0 = tri.v0 + d * t0;
    out.v1 = tri.v1 + d * t1;
    out.v2 = tri.v2 + d * t2;
    return true;
}

} // namespace geom
