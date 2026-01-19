#pragma once
#include <cmath>
#include <stdexcept>
#include "../geom/vec3.hpp"

namespace em {

    // Convenzioni coerenti con il testo:
    // - (x1,y1,z1) globale
    // - O2 = (x_af, 0, f2(x_af,0)) in coordinate x1
    // - x2 = x1 - x_af, y2 = y1, z2 = z1 - z0
    // - x3 ottenuto da x2 con: Rz(phi) poi Ry(theta)  (come nel documento)

    struct FrameParams {
        double x_af{ 0.0 };   // [m]
        double theta{ 0.0 };  // [rad]
        double phi{ 0.0 };    // [rad]
        double z0{ 0.0 };     // z0 = f2(x_af,0)  [m]
    };

    // Piccola matrice 3x3 solo per rotazioni
    struct Mat3 {
        double m[3][3]{};

        static Mat3 Rz(double a) {
            Mat3 R{};
            double c = std::cos(a), s = std::sin(a);
            R.m[0][0] = c;  R.m[0][1] = s;  R.m[0][2] = 0;
            R.m[1][0] = -s; R.m[1][1] = c;  R.m[1][2] = 0;
            R.m[2][0] = 0;  R.m[2][1] = 0;  R.m[2][2] = 1;
            return R;
        }
        static Mat3 Ry(double a) {
            Mat3 R{};
            double c = std::cos(a), s = std::sin(a);
            R.m[0][0] = c;  R.m[0][1] = 0; R.m[0][2] = -s;
            R.m[1][0] = 0;  R.m[1][1] = 1; R.m[1][2] = 0;
            R.m[2][0] = s;  R.m[2][1] = 0; R.m[2][2] = c;
            return R;
        }

        geom::Vec3 apply(const geom::Vec3& v) const {
            return {
                m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
                m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
                m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z
            };
        }

        static Mat3 mul(const Mat3& A, const Mat3& B) {
            Mat3 C{};
            for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
                double s = 0;
                for (int k = 0; k < 3; ++k) s += A.m[i][k] * B.m[k][j];
                C.m[i][j] = s;
            }
            return C;
        }
    };

    // Gestore trasformazioni
    class Frames {
    public:
        explicit Frames(const FrameParams& p) : p_(p) {
            // Forward: x3 = Ry(theta)*Rz(phi)*x2
            R_x2_to_x3_ = Mat3::mul(Mat3::Ry(p_.theta), Mat3::Rz(p_.phi));

            // Inverse: x2 = Rz(-phi)*Ry(-theta)*x3
            R_x3_to_x2_ = Mat3::mul(Mat3::Rz(-p_.phi), Mat3::Ry(-p_.theta));
        }

        // x1 -> x2
        geom::Vec3 x2_from_x1(const geom::Vec3& p1) const {
            return { p1.x - p_.x_af, p1.y, p1.z - p_.z0 };
        }
        // x2 -> x1
        geom::Vec3 x1_from_x2(const geom::Vec3& p2) const {
            return { p2.x + p_.x_af, p2.y, p2.z + p_.z0 };
        }

        // x2 -> x3
        geom::Vec3 x3_from_x2(const geom::Vec3& p2) const {
            return R_x2_to_x3_.apply(p2);
        }
        // x3 -> x2
        geom::Vec3 x2_from_x3(const geom::Vec3& p3) const {
            return R_x3_to_x2_.apply(p3);
        }

        // x1 -> x3
        geom::Vec3 x3_from_x1(const geom::Vec3& p1) const {
            return x3_from_x2(x2_from_x1(p1));
        }
        // x3 -> x1
        geom::Vec3 x1_from_x3(const geom::Vec3& p3) const {
            return x1_from_x2(x2_from_x3(p3));
        }

        // Versore asse z3 espresso in coordinate x1 (solo rotazione, no traslazione)
        geom::Vec3 z3_axis_in_x1() const {
            // direzione in x3
            geom::Vec3 ez3{ 0,0,1 };
            // ez3 -> x2 -> x1 (x2 e x1 paralleli, quindi la rotazione è la stessa)
            geom::Vec3 ez2 = R_x3_to_x2_.apply(ez3);
            return geom::normalize(ez2);
        }

        const FrameParams& params() const { return p_; }

    private:
        FrameParams p_;
        Mat3 R_x2_to_x3_{};
        Mat3 R_x3_to_x2_{};
    };

} // namespace em
