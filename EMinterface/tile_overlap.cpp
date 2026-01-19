#define EMSOLVER_NO_MAIN
#include "em_solver.cpp"

#include <iostream>
#include <vector>
#include <array>

#include <algorithm>

#include "include/geom/triangle.hpp"
#include "include/geom/vec3.hpp"

#include "include/utils/constants.hpp"

#include "include/em/meniscus_model.hpp"
#include "include/em/frames.hpp"
#include "include/em/beam_mesh.hpp"
#include "include/em/meniscus_projection.hpp"
#include "include/em/postprocess_grid.hpp"
#include "include/em/tile_wrapper.hpp"

static geom::Vec3 tri_normal_pointing_with_incident(const geom::Triangle& t, const geom::Vec3& s_in) {
    geom::Vec3 n = geom::normal_unit(t);
    if (geom::dot(n, s_in) < 0.0) n = n * -1.0;
    return n;
}

int main() {
    try {
        // =========================
        // PARAMETRI ( demostrativi ) 
        // =========================
        // Cuvette (interno)
        const double Rin = 6.3e-3;      // [m]
        const double V = 1.0e-6;      // [m^3] 1 mL = 1e-6 m^3 (esempio)
        const double Dm = 1.0e-3;      // [m] profondità menisco (esempio)

        // Fascio 
        const double Rf = 2.0e-3;     // [m] raggio fascio
        const double x_af = 1.0e-3;     // [m] offset asse fascio
        const double theta = 10.0 * utils::pi / 180.0;
        const double phi = 30.0 * utils::pi / 180.0;

        // Mesh
        const int n_pr = 10;          // punti per raggio (>=2)
        const int n_rays = 6 * (n_pr - 1);  

        // Post-processing plane
        // z1 = Cpp con 0 <= Cpp < f2(0,0)=a_pm 
        const double Cpp = 0.6e-3;      // [m] esempio

        // Griglia output su (x1,y1) (copri tutto il disco interno)
        const int Nx = 250, Ny = 250;

        // Wave / media (metti i tuoi valori veri qui)
        const double f = 100e9; // [Hz]
        const double omega = 2.0 * utils::pi * f;

        em::MediaParams media;
        media.omega = omega;

        // Mezzo 1 (aria)
        media.eps1_r = { 1.0, 0.0 };
        media.mu1_r = { 1.0, 0.0 };
        media.sigma1 = 0.0;

        // Mezzo 2 (soluzione salina - esempi)
        media.eps2_r = { 80.0, -5.0 };  // esempio complesso
        media.mu2_r = { 1.0, 0.0 };
        media.sigma2 = 1.0;           // S/m (esempio)

        // Polarizzazione / ampiezza
        em::BeamParams beam;
        beam.A_inc = { 1.0, 0.0 };
        beam.isTE = true;


        const int beam_sign = -1; // direzione fascio -1 perchè diretto verso basso 
        // =========================

        // =========================
        // 1) Menisco
        // =========================
        em::ParaboloidMeniscus men(em::MeniscusParams{ V, Dm, Rin });
        const double apm = men.apm();

        if (!(Cpp >= 0.0 && Cpp < apm)) {
            std::cerr << "WARNING: Cpp fuori range suggerito [0, a_pm). a_pm=" << apm << "\n";
        }

        // z0 = f2(x_af,0)
        const double z0 = men.f2(x_af, 0.0);

        // =========================
        // 2) Frames x1<->x3
        // =========================
        em::Frames fr(em::FrameParams{ x_af, theta, phi, z0 });

        // Direzione asse z3 in coordinate x1
        geom::Vec3 z3_dir_x1 = fr.z3_axis_in_x1();
        geom::Vec3 s_in = z3_dir_x1 * (double)beam_sign;
        s_in = geom::normalize(s_in);

        beam.s_in_global = s_in;

        // =========================
        // 3) Mesh nel piano (x3,y3), z3=0
        // =========================
        em::Mesh2D mesh = em::generate_disk_mesh(Rf, n_pr, n_rays);

        // =========================
        // 4) Proiezione mesh sul menisco (vera)
        // =========================
        std::vector<geom::Vec3> nodes_on_meniscus(mesh.nodes.size());
        std::vector<char> ok(mesh.nodes.size(), 0);

        for (size_t i = 0; i < mesh.nodes.size(); ++i) {
            const double x3 = mesh.nodes[i].x;
            const double y3 = mesh.nodes[i].y;

            // Proietta lungo +z3 (poi la direzione fisica la gestiamo con beam_sign)
            em::ProjectionResult pr = em::project_node_to_meniscus(fr, men, x3, y3);
            if (pr.ok) {
                nodes_on_meniscus[i] = pr.p1;
                ok[i] = 1;
            }
        }

        // =========================
        // 5) Postprocess grid su z1=Cpp
        // =========================
        em::Grid2D grid(-Rin, Rin, -Rin, Rin, Nx, Ny);

        double P_total_in = 0.0;
        double P_total_on_plane = 0.0;
        int tri_used = 0;

        // =========================
        // 6) Loop triangoli: solve locale + deposita su griglia
        // =========================
        for (const auto& tri_idx : mesh.tris) {
            int i0 = tri_idx[0], i1 = tri_idx[1], i2 = tri_idx[2];
            if (!ok[i0] || !ok[i1] || !ok[i2]) continue;

            geom::Triangle tri;
            tri.v0 = nodes_on_meniscus[i0];
            tri.v1 = nodes_on_meniscus[i1];
            tri.v2 = nodes_on_meniscus[i2];

            // Solve locale su faccia triangolare
            em::TileSolveResult res = em::solve_on_triangle(tri, media, beam, /*flip_normal=*/true);

            // Potenza che entra nel mezzo 2 attraverso quella faccia:
            // P_in ≈ (S · n_interface) * Area_face
            geom::Vec3 n_int = tri_normal_pointing_with_incident(tri, s_in);
            const double A_face = geom::area(tri);
            const double Pin = std::max(0.0, geom::dot(res.S_avg, n_int)) * A_face;
            if (Pin <= 0.0) continue;

            P_total_in += Pin;

            // Proietta triangolo sul piano z1=Cpp lungo direzione di trasporto energia (Poynting)
            geom::Vec3 dir = res.d_poynting;
            if (std::abs(dir.z) < 1e-12) continue; // quasi parallelo al piano

            auto proj_to_plane = [&](const geom::Vec3& p)->std::pair<bool, geom::Vec3> {
                double t = (Cpp - p.z) / dir.z;
                if (t < 0.0) return { false, {0,0,0} };
                return { true, p + dir * t };
                };

            auto p0 = proj_to_plane(tri.v0);
            auto p1 = proj_to_plane(tri.v1);
            auto p2 = proj_to_plane(tri.v2);
            if (!p0.first || !p1.first || !p2.first) continue;

            std::array<geom::Vec2, 3> tri2d = {
                geom::Vec2{p0.second.x, p0.second.y},
                geom::Vec2{p1.second.x, p1.second.y},
                geom::Vec2{p2.second.x, p2.second.y}
            };

            // Area footprint sul piano (x1,y1)
            auto area2 = [&](const std::array<geom::Vec2, 3>& t2)->double {
                double a = (t2[0].x * (t2[1].y - t2[2].y)
                    + t2[1].x * (t2[2].y - t2[0].y)
                    + t2[2].x * (t2[0].y - t2[1].y)) * 0.5;
                return std::abs(a);
                };
            const double A_fp = area2(tri2d);
            if (A_fp < 1e-18) continue;

            // Attenuazione (prima versione): usa Im(k) lungo lo spostamento del baricentro
            geom::Vec3 c_face = (tri.v0 + tri.v1 + tri.v2) / 3.0;
            auto pc = proj_to_plane(c_face);
            if (!pc.first) continue;
            geom::Vec3 disp = pc.second - c_face;

            // alpha_vec = -Im(kp2)  (coerente con la tua scelta Im(kz)<=0 nel solver)
            geom::Vec3 alpha = {
                -std::imag(res.kp2_global.x),
                -std::imag(res.kp2_global.y),
                -std::imag(res.kp2_global.z)
            };
            double att = std::exp(-2.0 * geom::dot(alpha, disp));
            att = std::clamp(att, 0.0, 1.0);

            const double P_on_plane = Pin * att;
            P_total_on_plane += P_on_plane;

            // Intensità “uniforme” sul footprint: I = P / Area
            const double I = P_on_plane / A_fp;

            grid.add_triangle_uniform(tri2d, I);
            tri_used++;
        }

        // =========================
        // 7) Output
        // =========================
        grid.save_csv("field_map_cpp.csv");

        std::cout << "=== MENISCUS SIM ===\n";
        std::cout << "Rin=" << Rin << "  V=" << V << "  Dm=" << Dm << "\n";
        std::cout << "a_pm=" << apm << "  z0=f2(x_af,0)=" << z0 << "  Cpp=" << Cpp << "\n";
        std::cout << "Mesh: nodes=" << mesh.nodes.size() << "  tris=" << mesh.tris.size()
            << "  used=" << tri_used << "\n";
        std::cout << "Total power into medium2 (approx) = " << P_total_in << " W\n";
        std::cout << "Total power on plane z1=Cpp (approx) = " << P_total_on_plane << " W\n";
        std::cout << "Saved: field_map_cpp.csv  (x,y,value) where value ~ intensity [W/m^2]\n";

    }
    catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
