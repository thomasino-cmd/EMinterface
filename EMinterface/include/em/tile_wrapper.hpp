#pragma once

#include <complex>
#include <stdexcept>
#include <algorithm>

#include "../geom/triangle.hpp"
#include "../geom/vec3.hpp"

namespace em {

struct MediaParams {
    std::complex<double> eps1_r{1.0, 0.0};
    std::complex<double> mu1_r{1.0, 0.0};
    double sigma1{0.0};

    std::complex<double> eps2_r{1.0, 0.0};
    std::complex<double> mu2_r{1.0, 0.0};
    double sigma2{0.0};

    double omega{0.0};
};

struct BeamParams {
    geom::Vec3 s_in_global;                // unit direction in medium 1 towards the interface
    std::complex<double> A_inc{1.0, 0.0};  // complex amplitude (scalar)
    bool isTE{true};
};

struct LocalFrame {
    geom::Vec3 ex;   // tangential axis in the plane of incidence
    geom::Vec3 ey;   // orthogonal tangential axis
    geom::Vec3 ez;   // normal axis (points from medium 1 into medium 2)
    double theta_i{0.0}; // local incidence angle [rad]
};

struct TileSolveResult {
    // Global-space complex vectors (phasors)
    Vector3D Ap2_global;
    Vector3D kp2_global;

    // Useful directions (real unit vectors)
    geom::Vec3 d_phase;     // ~ Re(k)
    geom::Vec3 d_atten;     // ~ -Im(k)
    geom::Vec3 d_poynting;  // ~ <S>
    geom::Vec3 S_avg;  // [W/m^2] time-average Poynting (real)


    // Diagnostics
    double theta_i{0.0};
    double R{0.0}, T{0.0}, A{0.0};
};

inline double clamp01(double v) {
    return std::max(-1.0, std::min(1.0, v));
}

inline geom::Vec3 vec_re(const Vector3D& v) {
    return {std::real(v.x), std::real(v.y), std::real(v.z)};
}

inline geom::Vec3 vec_im(const Vector3D& v) {
    return {std::imag(v.x), std::imag(v.y), std::imag(v.z)};
}

inline Vector3D conj_vec(const Vector3D& v) {
    return Vector3D(std::conj(v.x), std::conj(v.y), std::conj(v.z));
}

// Rotate a complex vector from local components (x,y,z) into global coordinates using the basis (ex,ey,ez).
inline Vector3D rotate_local_to_global(const Vector3D& v_local,
                                       const geom::Vec3& ex,
                                       const geom::Vec3& ey,
                                       const geom::Vec3& ez) {
    const std::complex<double> vx = v_local.x;
    const std::complex<double> vy = v_local.y;
    const std::complex<double> vz = v_local.z;

    return Vector3D(
        vx * ex.x + vy * ey.x + vz * ez.x,
        vx * ex.y + vy * ey.y + vz * ez.y,
        vx * ex.z + vy * ey.z + vz * ez.z
    );
}

// Build local frame on the triangle.
inline LocalFrame make_local_frame(const geom::Triangle& tri,
                                   geom::Vec3 s_in_global_unit,
                                   bool flip_normal_to_match_incident = true,
                                   double eps = 1e-12) 
{
    s_in_global_unit = geom::normalize(s_in_global_unit, eps);

    geom::Vec3 ez = geom::normal_unit(tri, eps);

    // Make ez point such that dot(s_in, ez) >= 0 for conventional incidence.
    if (flip_normal_to_match_incident && geom::dot(s_in_global_unit, ez) < 0.0) {
        ez = ez * -1.0;
    }

    const double cos_th = clamp01(geom::dot(s_in_global_unit, ez));
    const double theta = std::acos(cos_th);

    // Tangential component of s_in defines the plane of incidence.
    geom::Vec3 ex = s_in_global_unit - ez * cos_th;
    if (geom::is_zero(ex, eps)) {
        // Normal incidence: pick any stable tangent.
        geom::Vec3 a = (std::abs(ez.z) < 0.9) ? geom::Vec3{0.0, 0.0, 1.0} : geom::Vec3{0.0, 1.0, 0.0};
        ex = geom::cross(a, ez);
    }
    ex = geom::normalize(ex, eps);

    geom::Vec3 ey = geom::cross(ez, ex);
    ey = geom::normalize(ey, eps);

    return LocalFrame{ex, ey, ez, theta};
}


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////


// Build incident field vector in the LOCAL solver coordinate system.

inline Vector3D make_incident_E_local(double theta_i, std::complex<double> A_inc, bool isTE) {
    if (isTE) {
        return Vector3D(0.0, A_inc, 0.0);
    }
    const double c = std::cos(theta_i);
    const double s = std::sin(theta_i);
    return Vector3D(A_inc * c, 0.0, A_inc * (-s));
}

// Compute time-average Poynting vector <S> (real) from E and k.
inline geom::Vec3 poynting_avg_real(const Vector3D& E, const Vector3D& k_vec,
                                   std::complex<double> mu_r, double omega) {
    if (omega <= 0.0) throw std::runtime_error("poynting_avg_real: omega must be > 0");

    const std::complex<double> mu_abs = mu_r * Medium::MU0;
    Vector3D H = (1.0 / (mu_abs * omega)) * cross(k_vec, E);
    Vector3D S = 0.5 * cross(E, conj_vec(H));
    return {std::real(S.x), std::real(S.y), std::real(S.z)};
}



/// <summary>
/// main function of the file
/// </summary>
/// <param name="tri"></param>
/// <param name="media"></param>
/// <param name="beam"></param>
/// <param name="flip_normal_to_match_incident"></param>
/// <param name="eps"></param>
/// <returns></returns>
inline TileSolveResult solve_on_triangle(const geom::Triangle& tri,
                                        const MediaParams& media,
                                        const BeamParams& beam,
                                        bool flip_normal_to_match_incident = true,
                                        double eps = 1e-12) {
    if (media.omega <= 0.0) throw std::runtime_error("solve_on_triangle: omega must be > 0");

    // 1) Local frame
    LocalFrame fr = make_local_frame(tri, beam.s_in_global, flip_normal_to_match_incident, eps);

    // 2) Local incident field
    Vector3D Ap1_local = make_incident_E_local(fr.theta_i, beam.A_inc, beam.isTE);

    // 3) lancia il solver originale in coordinate locali
    EMInterfaceSolver solver(
        media.eps1_r, media.mu1_r, media.sigma1,
        media.eps2_r, media.mu2_r, media.sigma2,
        media.omega, fr.theta_i,
        Ap1_local, beam.isTE
    );

    Solution sol_local = solver.solve();
    const Vector3D kp2_local = solver.get_kp2(); 

    // 4) Rotate to global
    TileSolveResult out;
    out.Ap2_global = rotate_local_to_global(sol_local.Ap2, fr.ex, fr.ey, fr.ez);
    out.kp2_global = rotate_local_to_global(kp2_local, fr.ex, fr.ey, fr.ez);

    out.theta_i = fr.theta_i;
    out.R = sol_local.R;
    out.T = sol_local.T;
    out.A = sol_local.A;

    // 5) Directions
    const geom::Vec3 k_re = vec_re(out.kp2_global);
    const geom::Vec3 k_im = vec_im(out.kp2_global);

    out.d_phase = geom::safe_normalize(k_re, eps);
    out.d_atten = geom::safe_normalize(k_im * (-1.0), eps); // direction of decay

    const geom::Vec3 S = poynting_avg_real(out.Ap2_global, out.kp2_global, media.mu2_r, media.omega);
    out.d_poynting = geom::safe_normalize(S, eps);
    out.S_avg = S;

    return out;
}

} // namespace em
