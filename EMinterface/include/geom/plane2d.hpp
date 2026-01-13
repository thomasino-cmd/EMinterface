#pragma once
#include "vec2.hpp"
#include "vec3.hpp"
#include "basis.hpp"
#include "triangle.hpp"
#include <vector>

namespace geom {

inline Vec2 to2D(const Vec3& p, const Vec3& origin, const OrthoBasis& b) {
    const Vec3 d = p - origin;
    return {dot(d, b.e1), dot(d, b.e2)};
}

inline std::vector<Vec2> triangle_to2D(const Triangle& t, const Vec3& origin, const OrthoBasis& b) {
    return {to2D(t.v0, origin, b), to2D(t.v1, origin, b), to2D(t.v2, origin, b)};
}

} // namespace geom
