#pragma once
#include "vec3.hpp"

namespace geom {

struct OrthoBasis {
    Vec3 e1; // in-plane
    Vec3 e2; // in-plane
    Vec3 n;  // unit normal
};

// Build a stable orthonormal basis given a (possibly non-unit) normal.
// e1,e2 span the plane orthogonal to n.
inline OrthoBasis make_plane_basis(const Vec3& normal, double eps = 1e-12) {
    OrthoBasis b;
    b.n = normalize(normal, eps);

    // Pick a helper vector not parallel to n.
    Vec3 a = (std::abs(b.n.z) < 0.9) ? Vec3{0.0, 0.0, 1.0} : Vec3{0.0, 1.0, 0.0};

    b.e1 = normalize(cross(a, b.n), eps);
    b.e2 = cross(b.n, b.e1); // already unit
    return b;
}

inline Vec3 to_basis_coords(const Vec3& p, const Vec3& origin, const OrthoBasis& b) {
    const Vec3 d = p - origin;
    return {dot(d, b.e1), dot(d, b.e2), dot(d, b.n)};
}

inline Vec3 from_basis_coords(const Vec3& c, const Vec3& origin, const OrthoBasis& b) {
    return origin + b.e1*c.x + b.e2*c.y + b.n*c.z;
}

} // namespace geom
