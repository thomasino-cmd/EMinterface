#define EMSOLVER_NO_MAIN
#include "em_solver.cpp"

#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <iomanip>
#include <limits>
#include <cmath>

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
        const double V = 3.0e-7;      // [m^3] 300 mm^3
        const double Dm = 1.0e-3;      // [m] profondità menisco (esempio)

        // Fascio 
        const double Rf = 5.64e-3;     // [m] raggio fascio USUAL VALUE: Rf = 5.64e-3;  

        //DISALLINEAMENTO LATERALE
        const double x_af = 0.0e-3;     // [m] offset asse fascio oscilla tra 0.0mm E 0.5mm
        const double theta = 0 * utils::pi / 180.0;      //10.0 gradi di inclinazione rispetto alla verticale (esempio)
        const double phi = 0 * utils::pi / 180.0;        //30 gradi di rotazione attorno all'asse z (esempio) questo così è inutile perchè il fascio è circolare quindi simmetrico 

        // Mesh
        const int n_pr = 10;          // punti per raggio (>=2)
        const int n_rays = 6 * (n_pr - 1);  

        // Post-processing plane
        // z1 = Cpp con 0 <= Cpp < f2(0,0)=a_pm 
        const double Cpp = 0.0 ;      // [m] esempio 0.6e-3

        // Griglia output su (x1,y1) (copri tutto il disco interno)
        const int Nx = 250, Ny = 250;

        // Wave / media (metti i tuoi valori veri qui)
        const double f = 3.701e14; // [Hz] esempio 100 GHz ((100e9)) (onda millimetrica), lambda ~ 3 mm
        const double omega = 2.0 * utils::pi * f;

        em::MediaParams media;
        media.omega = omega;

        // Mezzo 1 (aria)
        media.eps1_r = { 1.0, 0.0 };
        media.mu1_r = { 1.0, 0.0 };
        media.sigma1 = 0.0;

        // Mezzo 2 (soluzione salina - esempi)
        media.eps2_r = { 1.7716, -3.346e-7 };  
        media.mu2_r = { 1.0, 0.0 };
        media.sigma2 = 1.0;           // S/m (esempio 1.0 )

        // Polarizzazione / ampiezza
        em::BeamParams beam;
        beam.A_inc = { 2744.9 , 0.0 };  // ampiezza complessa (esempio 1.0)
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
        std::ofstream node_status_csv("node_status.csv");
        if (!node_status_csv) throw std::runtime_error("Cannot open node_status.csv");
        node_status_csv << "node_id,x3_m,y3_m,x1_m,y1_m,ok,reason_code,reason_label,"
                           "disc,root1,root2,chosen_root,checks\n";

        std::vector<geom::Vec3> nodes_on_meniscus(mesh.nodes.size());
        std::vector<char> ok(mesh.nodes.size(), 0);

        for (size_t i = 0; i < mesh.nodes.size(); ++i) {
            const double x3 = mesh.nodes[i].x;
            const double y3 = mesh.nodes[i].y;

            // Coordinate del nodo nel sistema globale sul piano z3=0
            const geom::Vec3 p1 = fr.x1_from_x3({ x3, y3, 0.0 });

            // Proietta lungo +z3 (poi la direzione fisica la gestiamo con beam_sign)
            em::ProjectionResult pr = em::project_node_to_meniscus(fr, men, x3, y3);
            if (pr.ok) {
                nodes_on_meniscus[i] = pr.p1;
                ok[i] = 1;
            }

            node_status_csv << i << "," << std::setprecision(16)
                << x3 << "," << y3 << ","
                << p1.x << "," << p1.y << ","
                << (pr.ok ? 1 : 0) << ","
                << pr.reason_code << "," << pr.reason_label << ","
                << pr.disc << ","
                << pr.root1 << ","
                << pr.root2 << ","
                << pr.chosen_root << ","
                << std::quoted(pr.checks) << "\n";
        }

        // =========================
        // 5) Postprocess grid su z1=Cpp
        // =========================
        em::Grid2D grid(-Rin, Rin, -Rin, Rin, Nx, Ny);

        // Output dettagliato per confronto menisco vs piano fondo (per triangolo)
        std::ofstream meniscus_csv("triangle_power_meniscus.csv");
        if (!meniscus_csv) throw std::runtime_error("Cannot open triangle_power_meniscus.csv");
        meniscus_csv << "tri_id,cx_m,cy_m,cz_m,area_face_m2,Sdotn_Wm2,Pin_W\n";

        std::ofstream bottom_csv("triangle_power_bottom.csv");
        if (!bottom_csv) throw std::runtime_error("Cannot open triangle_power_bottom.csv");
        bottom_csv << "tri_id,foot_cx_m,foot_cy_m,foot_area_m2,path_len_m,alpha_dot_disp,att,Ponplane_W,I_Wm2\n";

        // Geometria triangoli per plotting continuo (menisco 3D e fondo 2D)
        std::ofstream meniscus_mesh_csv("triangle_mesh_meniscus.csv");
        if (!meniscus_mesh_csv) throw std::runtime_error("Cannot open triangle_mesh_meniscus.csv");
        meniscus_mesh_csv << "tri_id,v0x_m,v0y_m,v0z_m,v1x_m,v1y_m,v1z_m,v2x_m,v2y_m,v2z_m,area_face_m2,Sdotn_Wm2,Pin_W,R,T,A\n";

        std::ofstream bottom_mesh_csv("triangle_mesh_bottom.csv");
        if (!bottom_mesh_csv) throw std::runtime_error("Cannot open triangle_mesh_bottom.csv");
        bottom_mesh_csv << "tri_id,valid,p0x_m,p0y_m,p1x_m,p1y_m,p2x_m,p2y_m,foot_area_m2,att,Ponplane_W,I_Wm2,status_code,status_label\n";

        std::ofstream tri_status_csv("triangle_status.csv");
        if (!tri_status_csv) throw std::runtime_error("Cannot open triangle_status.csv");
        tri_status_csv << "tri_id,status_code,status_label,has_meniscus_geom,has_bottom_geom,Pin_W,Ponplane_W,att\n";

        std::ofstream tri_fail_nodes_csv("tri_fail_nodes.csv");
        if (!tri_fail_nodes_csv) throw std::runtime_error("Cannot open tri_fail_nodes.csv");
        tri_fail_nodes_csv << "tri_id,i0_ok,i1_ok,i2_ok,i0_id,i1_id,i2_id\n";

        const double nanv = std::numeric_limits<double>::quiet_NaN();

        double P_total_in = 0.0;
        double P_total_on_plane = 0.0;
        int tri_used = 0;
        int tri_id_counter = 0;

        // =========================
        // 6) Loop triangoli: solve locale + deposita su griglia
        // =========================
        for (const auto& tri_idx : mesh.tris) {
            const int tri_id = tri_id_counter++;
            int i0 = tri_idx[0], i1 = tri_idx[1], i2 = tri_idx[2];

            if (!ok[i0] || !ok[i1] || !ok[i2]) {
                tri_fail_nodes_csv << tri_id << ","
                    << int(ok[i0]) << "," << int(ok[i1]) << "," << int(ok[i2]) << ","
                    << i0 << "," << i1 << "," << i2 << "\n";
                tri_status_csv << tri_id << ",1,invalid_meniscus_nodes,0,0,"
                    << nanv << "," << nanv << "," << nanv << "\n";
                bottom_mesh_csv << tri_id << ",0,"
                    << nanv << "," << nanv << "," << nanv << "," << nanv << ","
                    << nanv << "," << nanv << "," << nanv << "," << nanv << ","
                    << nanv << "," << nanv << ",1,invalid_meniscus_nodes\n";
                continue;
            }

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
            const double Sdotn = geom::dot(res.S_avg, n_int);
            const double Pin = std::max(0.0, Sdotn) * A_face;

            const geom::Vec3 c_face = (tri.v0 + tri.v1 + tri.v2) / 3.0;
            meniscus_csv << tri_id << ","
                << std::setprecision(16)
                << c_face.x << "," << c_face.y << "," << c_face.z << ","
                << A_face << "," << Sdotn << "," << Pin << "\n";

            meniscus_mesh_csv << tri_id << ","
                << std::setprecision(16)
                << tri.v0.x << "," << tri.v0.y << "," << tri.v0.z << ","
                << tri.v1.x << "," << tri.v1.y << "," << tri.v1.z << ","
                << tri.v2.x << "," << tri.v2.y << "," << tri.v2.z << ","
                << A_face << "," << Sdotn << "," << Pin << ","
                << res.R << "," << res.T << "," << res.A << "\n";

            if (Pin <= 0.0) {
                tri_status_csv << tri_id << ",2,pin_nonpositive,1,0,"
                    << Pin << "," << nanv << "," << nanv << "\n";
                bottom_mesh_csv << tri_id << ",0,"
                    << nanv << "," << nanv << "," << nanv << "," << nanv << ","
                    << nanv << "," << nanv << "," << nanv << "," << nanv << ","
                    << nanv << "," << nanv << ",2,pin_nonpositive\n";
                continue;
            }

            P_total_in += Pin;

            // Proietta triangolo sul piano z1=Cpp lungo direzione di trasporto energia (Poynting)
            geom::Vec3 dir = res.d_poynting;
            if (std::abs(dir.z) < 1e-12) {
                tri_status_csv << tri_id << ",3,dir_parallel_to_plane,1,0,"
                    << Pin << "," << nanv << "," << nanv << "\n";
                bottom_mesh_csv << tri_id << ",0,"
                    << nanv << "," << nanv << "," << nanv << "," << nanv << ","
                    << nanv << "," << nanv << "," << nanv << "," << nanv << ","
                    << nanv << "," << nanv << ",3,dir_parallel_to_plane\n";
                continue;
            } // quasi parallelo al piano

            auto proj_to_plane = [&](const geom::Vec3& p)->std::pair<bool, geom::Vec3> {
                double t = (Cpp - p.z) / dir.z;
                if (t < 0.0) return { false, {0,0,0} };
                return { true, p + dir * t };
                };

            auto p0 = proj_to_plane(tri.v0);
            auto p1 = proj_to_plane(tri.v1);
            auto p2 = proj_to_plane(tri.v2);
            if (!p0.first || !p1.first || !p2.first) {
                tri_status_csv << tri_id << ",4,vertex_projection_fail,1,0,"
                    << Pin << "," << nanv << "," << nanv << "\n";
                bottom_mesh_csv << tri_id << ",0,"
                    << nanv << "," << nanv << "," << nanv << "," << nanv << ","
                    << nanv << "," << nanv << "," << nanv << "," << nanv << ","
                    << nanv << "," << nanv << ",4,vertex_projection_fail\n";
                continue;
            }

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
            if (A_fp < 1e-18) {
                tri_status_csv << tri_id << ",5,footprint_too_small,1,0,"
                    << Pin << "," << nanv << "," << nanv << "\n";
                bottom_mesh_csv << tri_id << ",0,"
                    << nanv << "," << nanv << "," << nanv << "," << nanv << ","
                    << nanv << "," << nanv << "," << nanv << "," << nanv << ","
                    << nanv << "," << nanv << ",5,footprint_too_small\n";
                continue;
            }

            // Attenuazione (prima versione): usa Im(k) lungo lo spostamento del baricentro
            auto pc = proj_to_plane(c_face);
            if (!pc.first) {
                tri_status_csv << tri_id << ",6,centroid_projection_fail,1,0,"
                    << Pin << "," << nanv << "," << nanv << "\n";
                bottom_mesh_csv << tri_id << ",0,"
                    << nanv << "," << nanv << "," << nanv << "," << nanv << ","
                    << nanv << "," << nanv << "," << nanv << "," << nanv << ","
                    << nanv << "," << nanv << ",6,centroid_projection_fail\n";
                continue;
            }
            geom::Vec3 disp = pc.second - c_face;
            const double path_len = geom::norm(disp);

            // alpha_vec = -Im(kp2)  (coerente con la tua scelta Im(kz)<=0 nel solver)
            geom::Vec3 alpha = {
                -std::imag(res.kp2_global.x),
                -std::imag(res.kp2_global.y),
                -std::imag(res.kp2_global.z)
            };
            const double alpha_dot_disp = geom::dot(alpha, disp);
            double att = std::exp(-2.0 * alpha_dot_disp);
            att = std::clamp(att, 0.0, 1.0);

            const double P_on_plane = Pin * att;
            P_total_on_plane += P_on_plane;

            // Intensità “uniforme” sul footprint: I = P / Area
            const double I = P_on_plane / A_fp;

            const geom::Vec2 c_fp{
                (tri2d[0].x + tri2d[1].x + tri2d[2].x) / 3.0,
                (tri2d[0].y + tri2d[1].y + tri2d[2].y) / 3.0
            };

            bottom_csv << tri_id << ","
                << std::setprecision(16)
                << c_fp.x << "," << c_fp.y << ","
                << A_fp << "," << path_len << ","
                << alpha_dot_disp << "," << att << ","
                << P_on_plane << "," << I << "\n";

            bottom_mesh_csv << tri_id << ",1,"
                << std::setprecision(16)
                << tri2d[0].x << "," << tri2d[0].y << ","
                << tri2d[1].x << "," << tri2d[1].y << ","
                << tri2d[2].x << "," << tri2d[2].y << ","
                << A_fp << "," << att << "," << P_on_plane << "," << I
                << ",0,ok" << "\n";

            tri_status_csv << tri_id << ",0,ok,1,1,"
                << Pin << "," << P_on_plane << "," << att << "\n";

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
        std::cout << "Saved: triangle_power_meniscus.csv  (per-triangle interface power diagnostics)\n";
        std::cout << "Saved: triangle_power_bottom.csv  (per-triangle delivered power at bottom plane)\n";
        std::cout << "Saved: triangle_mesh_meniscus.csv  (3D meniscus triangles with power diagnostics)\n";
        std::cout << "Saved: triangle_mesh_bottom.csv  (2D bottom footprints with delivered power)\n";
        std::cout << "Saved: triangle_status.csv  (per-triangle validity/exclusion reason)\n";
        std::cout << "Saved: node_status.csv  (per-node meniscus projection diagnostics)\n";
        std::cout << "Saved: tri_fail_nodes.csv  (triangle -> failing node mapping)\n";

    }
    catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
