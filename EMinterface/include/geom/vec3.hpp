#pragma once
#include <cmath>
#include <stdexcept>

namespace geom {

struct Vec3 {
    double x{0.0}, y{0.0}, z{0.0};

    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }

    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
};

inline double dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
}

inline double norm2(const Vec3& v) { return dot(v, v); }
inline double norm(const Vec3& v) { return std::sqrt(norm2(v)); }

inline Vec3 normalize(const Vec3& v, double eps = 1e-12) {
    const double n = norm(v);
    if (n < eps) throw std::runtime_error("normalize(): zero-length vector");
    return v / n;
}

inline Vec3 safe_normalize(const Vec3& v, double eps = 1e-12) {
    const double n = norm(v);
    if (n < eps) return {0.0, 0.0, 0.0};
    return v / n;
}

inline bool is_zero(const Vec3& v, double eps = 1e-12) {
    return norm2(v) < eps*eps;
}

inline double clamp(double x, double lo, double hi) {
    return (x < lo) ? lo : (x > hi) ? hi : x;
}

} // namespace geom
