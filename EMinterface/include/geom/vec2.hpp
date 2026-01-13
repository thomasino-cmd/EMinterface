#pragma once
#include <cmath>

namespace geom {

struct Vec2 {
    double x{0.0}, y{0.0};

    Vec2() = default;
    Vec2(double x_, double y_) : x(x_), y(y_) {}

    Vec2 operator+(const Vec2& o) const { return {x + o.x, y + o.y}; }
    Vec2 operator-(const Vec2& o) const { return {x - o.x, y - o.y}; }
    Vec2 operator*(double s) const { return {x * s, y * s}; }
    Vec2 operator/(double s) const { return {x / s, y / s}; }
};

inline double dot(const Vec2& a, const Vec2& b) { return a.x * b.x + a.y * b.y; }
inline double cross_z(const Vec2& a, const Vec2& b) { return a.x * b.y - a.y * b.x; }
inline double norm2(const Vec2& v) { return dot(v, v); }
inline double norm(const Vec2& v) { return std::sqrt(norm2(v)); }

} // namespace geom
