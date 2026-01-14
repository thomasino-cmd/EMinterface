// Demo: compute overlap between two projected triangles on the bottom plane.

#define EMSOLVER_NO_MAIN
#include "em_solver.cpp"         
#include "include/em/meniscus_overlap.hpp"
#include <iostream>
#include <iomanip>

int main() {
    try {
        //stai prendendo un tile  1mm * 1mm
        //approssimato con un triangolo
        //posto su z = 0
        // 1) Define one interface triangle (units: meters)
        geom::Triangle tri;
        tri.v0 = {0.0, 0.0, 0.0};
        tri.v1 = {1e-3, 0.0, 0.0};
        tri.v2 = {0.0, 1e-3, 0.0};


        // 2) parametri dei mezzi settati in advance così evito  input a runtime (SOLO PER DEMO)

        em::MediaParams mp;
        mp.omega  = 2.0 * PI * 3e14; // ~ 1 micron light (rough), just as a demo
        mp.eps1_r  = {1.0, 0.0};
        mp.mu1_r   = {1.0, 0.0};
        mp.sigma1  = 0.0;
        mp.eps2_r  = {1.77, -0.02};   // lossy-ish saline-ish placeholder RICORDA CHE EPS'' DEVE ESSERE NEGATIVO!3
        mp.mu2_r   = {1.0, 0.0};
        mp.sigma2  = 0.0;

        // 3) Beam parameters
        em::BeamParams bp;
        bp.isTE = true;
        bp.A_inc = {1.0, 0.0};
        bp.s_in_global = geom::normalize({0.2, 0.0, 1.0}); // oblique incidence

        // 4) risolve la fisica sul triangolo in tile_wrapper.hpp che chiama il em_solver.cpp (solo la classe)
        const em::TileSolveResult r = em::solve_on_triangle(tri, mp, bp);

        // Choose which two directions you want to compare
        const geom::Vec3 dir_a = r.d_phase;     // phase direction
        const geom::Vec3 dir_b = r.d_atten;     // attenuation direction


        // 5) Define a bottom plane: normal along +z, at thickness h (meters)
        const geom::Vec3 nb = geom::normalize({0.0, 0.0, 1.0});

        const double h = 1.0e-3;

        const geom::Plane bottom = em::bottom_plane_from_triangle_centroid(tri, nb, h);

        // 6) Compute overlap
        const em::OverlapMetrics m = em::triangle_overlap_on_plane(tri, bottom, dir_a, dir_b);

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "theta_i (local) = " << r.theta_i << " rad\n";
        std::cout << "R=" << r.R << "  T=" << r.T << "  A=" << r.A << "\n";
        std::cout << "AreaA=" << m.area_a << "  AreaB=" << m.area_b << "  AreaI=" << m.area_intersection << "\n";
        std::cout << "IoU=" << m.iou << "  overlap_min=" << m.overlap_min << "\n";


        // 7) find maximum thickness satisfying IoU >= 0.95
        const double hmax = em::find_max_thickness_for_overlap(tri, nb, dir_a, dir_b, 0.95, 0.0, 5e-3);
        std::cout << "h_max(IoU>=0.95) ~= " << hmax << " m\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
