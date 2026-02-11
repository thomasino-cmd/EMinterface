#pragma once
#include <cmath>
#include <limits>
#include <string>
#include <stdexcept>
#include <algorithm>


#include "../geom/vec2.hpp"
#include "../geom/vec3.hpp"
#include "frames.hpp"
#include "meniscus_model.hpp"

namespace em {

    struct ProjectionResult {
        bool ok{ false };
        double t_z3{ 0.0 };       // valore di z3 trovato
        geom::Vec3 p1;          // punto su menisco in coordinate x1

        int reason_code{ 0 };
        const char* reason_label{ "ok" };

        double disc{ std::numeric_limits<double>::quiet_NaN() };
        double root1{ std::numeric_limits<double>::quiet_NaN() };
        double root2{ std::numeric_limits<double>::quiet_NaN() };
        double chosen_root{ std::numeric_limits<double>::quiet_NaN() };

        std::string checks;
    };

    // Risolve a*t^2 + b*t + c = 0 e sceglie la radice più piccola >=0
    inline bool smallest_nonneg_root(double a, double b, double c, double& t_out, double eps = 1e-14) {
        if (std::abs(a) < eps) {
            // lineare
            if (std::abs(b) < eps) return false;
            double t = -c / b;
            if (t < 0.0) return false;
            t_out = t;
            return true;
        }

        double disc = b * b - 4 * a * c;
        if (disc < 0.0) return false;

        double s = std::sqrt(std::max(0.0, disc));
        double t1 = (-b - s) / (2 * a);
        double t2 = (-b + s) / (2 * a);

        double best = std::numeric_limits<double>::infinity();
        if (t1 >= 0.0) best = std::min(best, t1);
        if (t2 >= 0.0) best = std::min(best, t2);
        if (!std::isfinite(best)) return false;

        t_out = best;
        return true;
    }

    // Proietta un nodo della mesh (x3,y3,0) lungo +z3 fino a colpire il menisco.
    inline ProjectionResult project_node_to_meniscus(
        const Frames& fr,
        const ParaboloidMeniscus& men,
        double x3, double y3,
        double eps = 1e-12
    ) {
        ProjectionResult out;

        auto mark_fail = [&](int code, const char* label, const std::string& checks = std::string{}) {
            out.ok = false;
            out.reason_code = code;
            out.reason_label = label;
            out.checks = checks;
            return out;
            };

        const auto& P = fr.params();
        const double th = P.theta;
        const double ph = P.phi;

        const double cth = std::cos(th), sth = std::sin(th);
        const double cph = std::cos(ph), sph = std::sin(ph);

        // Parametrizzazione punto in x1 come funzione di t=z3:
        // (x3,y3,t) -> (x1(t),y1(t),z1(t))
        //
        // Derivata direttamente dalle formule inverse che hai scritto:
        // x1 = x_af + (x3 cosθ + t sinθ) cosφ - y3 sinφ
        // y1 =        (x3 cosθ + t sinθ) sinφ + y3 cosφ
        // z1 = z0 + (-x3 sinθ + t cosθ)
        //
        const double x1_0 = P.x_af + (x3 * cth) * cph - y3 * sph;
        const double x1_1 = (sth)*cph;

        const double y1_0 = (x3 * cth) * sph + y3 * cph;
        const double y1_1 = (sth)*sph;

        const double z1_0 = P.z0 + (-x3 * sth);
        const double z1_1 = cth;

        // f2 = k*(x1^2 + y1^2) + apm
        // imponi z1(t)=f2(x1(t),y1(t)) -> quadratica in t
        const double k = men.k();
        const double apm = men.apm();

        // r^2(t)= (x1_0+x1_1 t)^2 + (y1_0+y1_1 t)^2 = c0 + c1 t + c2 t^2
        const double c2 = x1_1 * x1_1 + y1_1 * y1_1; // = sin^2(theta)
        const double c1 = 2.0 * (x1_0 * x1_1 + y1_0 * y1_1);
        const double c0 = x1_0 * x1_0 + y1_0 * y1_0;

        // k*c2 t^2 + (k*c1 - z1_1)t + (k*c0 + apm - z1_0)=0
        const double A = k * c2;
        const double B = k * c1 - z1_1;
        const double C = k * c0 + apm - z1_0;

        const double disc = B * B - 4.0 * A * C;
        out.disc = disc;
        if (!std::isfinite(disc)) return mark_fail(11, "disc_nan");
        if (disc < 0.0) return mark_fail(12, "disc_negative");

        if (std::abs(A) < 1e-14) {
            if (std::abs(B) < 1e-14) return mark_fail(13, "degenerate_equation");
            const double t = -C / B;
            out.root1 = t;
            out.root2 = std::numeric_limits<double>::quiet_NaN();
            out.chosen_root = t;
            out.t_z3 = t;
            if (t < 0.0) {
                out.checks = "t_negative_used";
            }
                 
        }
        else {
            const double sqrt_disc = std::sqrt(std::max(0.0, disc));
            const double t1 = (-B - sqrt_disc) / (2.0 * A);
            const double t2 = (-B + sqrt_disc) / (2.0 * A);
            out.root1 = t1;
            out.root2 = t2;

            const bool t1_ok = t1 >= 0.0;
            const bool t2_ok = t2 >= 0.0;

            double t;
            if (t1_ok && t2_ok) {
                t = std::min(t1, t2);
            }
            else if (t1_ok) {
                t = t1;
            }
            else if (t2_ok) {
                t = t2;
            }
            else {
                // Both roots are negative: the meniscus lies in the -z3 direction from the z3=0 plane.
                // Accept the intersection anyway by picking the root closest to zero (largest t).
                t = std::max(t1, t2);
                out.checks = "t_negative_used";
            }
            out.chosen_root = t;
            out.t_z3 = t;

        }

        const double t = out.t_z3;

        // Punto su menisco
        const double x1 = x1_0 + x1_1 * t;
        const double y1 = y1_0 + y1_1 * t;
        const double z1 = z1_0 + z1_1 * t;

        // Check r <= Rin (opzionale ma consigliato)
        if (!men.inside_cuvette(x1, y1, eps)) {
            return mark_fail(10, "outside_Rin", "r>Rin");
        }

        const double z_expected = men.f2(x1, y1);
        if (std::abs(z1 - z_expected) > 1e-9) {
            return mark_fail(15, "surface_mismatch", "abs(z1-f2)>1e-9");
        }

        out.ok = true;
        out.t_z3 = t;
        out.p1 = { x1, y1, z1 };
        out.reason_code = 0;
        out.reason_label = "ok";
        return out;
    }

} // namespace em
