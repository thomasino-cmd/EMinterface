#pragma once
#include "vec3.hpp"

namespace geom {

// Plane represented as: n · (p - p0) = 0
// Convention: n points to the "positive" side.
struct Plane {
    Vec3 n;   // preferably unit
    Vec3 p0;  // a point on the plane

    static Plane fromPointNormal(const Vec3& p0_in, const Vec3& n_unit) {
        return Plane{n_unit, p0_in};
    }
};

inline double signed_distance(const Plane& pl, const Vec3& p) {
    return dot(pl.n, p - pl.p0);
}

// Ray-plane intersection: p(t) = p + t*d
// Returns true if not parallel. t_out is written.
inline bool intersect_ray_plane(const Vec3& p, const Vec3& d, const Plane& pl, double& t_out, double eps = 1e-12) {
    const double denom = dot(pl.n, d);
    if (std::abs(denom) < eps) return false;
    t_out = dot(pl.n, (pl.p0 - p)) / denom;
    return true;
}

} // namespace geom
