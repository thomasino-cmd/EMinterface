#pragma once
#include "vec3.hpp"

namespace geom {

struct Triangle {
    Vec3 v0, v1, v2;
};

inline Vec3 normal_unnormalized(const Triangle& t) {
    return cross(t.v1 - t.v0, t.v2 - t.v0);
}

inline Vec3 normal_unit(const Triangle& t, double eps = 1e-12) {
    return normalize(normal_unnormalized(t), eps);
}

inline double area(const Triangle& t) {
    return 0.5 * norm(normal_unnormalized(t));
}

inline Vec3 centroid(const Triangle& t) {
    return (t.v0 + t.v1 + t.v2) / 3.0;
}

} // namespace geom
