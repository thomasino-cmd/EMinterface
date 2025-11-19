#include <iostream>     
#include <sstream> 
#include <iomanip>      
#include <string>      
#include <vector>       
#include <complex>      
#include <cmath>        
#include <stdexcept>    
#include <algorithm>    
#include <exception>    
#include <Eigen/Dense> 





using namespace std;
using namespace Eigen; 
using complexd = complex<double>;
const double EPS = 1e-12;
const double PI = acos(-1.0);

struct Vector3D {
    complexd x = 0.0, y = 0.0, z = 0.0;
    Vector3D() = default;
    Vector3D(complexd _x, complexd _y, complexd _z) : x(_x), y(_y), z(_z) {}
    Vector3D operator+(const Vector3D& o) const { return Vector3D(x + o.x, y + o.y, z + o.z); }
    Vector3D operator-(const Vector3D& o) const { return Vector3D(x - o.x, y - o.y, z - o.z); }
    Vector3D operator*(const complexd& s) const { return Vector3D(x * s, y * s, z * s); }
    friend Vector3D operator*(const complexd& s, const Vector3D& v) { return v * s; }
};

complexd dot(const Vector3D& a, const Vector3D& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vector3D cross(const Vector3D& a, const Vector3D& b) {
    return Vector3D(a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x);
}
string cstr(const complexd& c) {
    ostringstream ss;
    ss << fixed << setprecision(6);
    ss << real(c);
    if (imag(c) >= 0) ss << " + j" << imag(c);
    else ss << " - j" << abs(imag(c));
    return ss.str();
}

struct Medium {
    complexd eps; // absolute permittivity
    complexd mu;  // absolute permeability
    complexd k;   // scalar wave number (complex)
    complexd Z;   // wave impedance sqrt(mu/eps)
    Medium() = default;
    // eps_r and mu_r are relative complex values; omega in rad/s
    Medium(complexd eps_r, complexd mu_r, double omega) {
        const double EPS0 = 8.854187817e-12;
        const double MU0 = 1.2566370614359173e-6; // correct mu0
        eps = eps_r * EPS0;
        mu = mu_r * MU0;
        // k = omega * sqrt(mu * eps)
        complexd prod = mu * eps;
        complexd s = sqrt(prod);
        k = complexd(omega) * s;
        // choose sign of k such that Im(k) <= 0 if possible (attenuating)
        if (imag(k) > 0) k = -k;
        Z = sqrt(mu / eps);
    }
};

struct Solution {
    Vector3D Ar1;
    Vector3D Ap2;
    double R = 0.0, T = 0.0, A = 0.0;
    double theta_t_deg = 0.0;
};

class EMInterfaceSolver {
public:
    EMInterfaceSolver(complexd eps1_r, complexd mu1_r,
        complexd eps2_r, complexd mu2_r,
        double omega, double theta_i_deg,
        const Vector3D& Ap1_in, bool isTE)
        : m1(eps1_r, mu1_r, omega), m2(eps2_r, mu2_r, omega),
        omega(omega), theta_i(theta_i_deg* PI / 180.0), Ap1(Ap1_in), isTE(isTE)
    {
        calculate_wave_vectors();
        // TEM check: Ap1 · kp1 == 0 (within tolerance)
        complexd d = dot(Ap1, kp1);
        if (abs(d) > 1e-8) {
            throw runtime_error("L'onda incidente non è TEM: Ap1 · kp1 != 0 (valore = " + cstr(d) + "). Abort.");
        }
    }

    Solution solve() {

        /*
        cerr << "DEBUG k-vectors:\n";
        cerr << "k1 = " << m1.k << "\n";
        cerr << "k2 = " << m2.k << "\n";
        cerr << "kp1 = (" << kp1.x << ", " << kp1.y << ", " << kp1.z << ")\n";
        cerr << "kp2 = (" << kp2.x << ", " << kp2.y << ", " << kp2.z << ")\n";
        cerr << "kx = " << kp1.x << "\n";
        */
        if (isTE) return solve_TE();
        else return solve_TM();
    }

private:
    Medium m1, m2;
    double omega;
    double theta_i;
    Vector3D Ap1;
    bool isTE;

    // wave vectors (only x and z are used)
    Vector3D kp1, kr1, kp2;

    void calculate_wave_vectors() {
        complexd k1 = m1.k;
        complexd k2 = m2.k;
        // tangential component kx
        complexd kx = k1 * sin(theta_i); // conserved
        // kp1z = sqrt(k1^2 - kx^2), choose root with Re(kp1z) > 0 if possible (incident towards +z)
        complexd arg1 = k1 * k1 - kx * kx;
        complexd kz1 = sqrt(arg1);
        if (real(kz1) < 0) kz1 = -kz1;
        // If kz1 has positive imag we might flip sign only if that makes sense; keep real>0 for incident
        if (imag(kz1) > 0 && real(kz1) >= 0) kz1 = -kz1;

        kp1 = Vector3D(kx, complexd(0.0), kz1);
        // reflected has same kx and opposite kz
        kr1 = Vector3D(kx, complexd(0.0), -kz1);

        // kp2z = sqrt(k2^2 - kx^2) with Im(kp2z) <= 0
        complexd arg2 = k2 * k2 - kx * kx;
        complexd kz2 = sqrt(arg2);
        if (imag(kz2) > 0) kz2 = -kz2; // enforce Im <= 0
        if (abs(imag(kz2)) < 1e-15 && real(kz2) < 0) kz2 = -kz2; // prefer positive real part when purely real
        kp2 = Vector3D(kx, complexd(0.0), kz2);
    }

    


    //
    //// Solve TE case using 2x2 linear system 
    //Solution solve_TE() {
    //    // unknowns: Ar1y, Ap2y
    //    complexd Ap1y = Ap1.y;
    //    complexd kp1z = kp1.z;
    //    complexd kr1z = kr1.z;
    //    complexd kp2z = kp2.z;
    //    complexd mu1 = m1.mu;
    //    complexd mu2 = m2.mu;

    //    // Equations:
    //    // 1) Ap1y + Ar1y - Ap2y = 0
    //    // 2) (kp1z*Ap1y + kr1z*Ar1y)/mu1 - (kp2z*Ap2y)/mu2 = 0

    //    // matrix form M * [Ar1y, Ap2y]^T = rhs
    //    // from 1) Ar1y - Ap2y = -Ap1y
    //    // from 2) (kr1z/mu1) * Ar1y + (-kp2z/mu2) * Ap2y = -kp1z*Ap1y / mu1

    //    complexd a11 = complexd(1.0);               // Ar1y coef from eq1
    //    complexd a12 = complexd(-1.0);              // Ap2y coef from eq1
    //    complexd b1 = -Ap1y;

    //    complexd a21 = kr1z / mu1;
    //    complexd a22 = -kp2z / mu2;
    //    complexd b2 = -(kp1z * Ap1y) / mu1;

    //    // Solve 2x2
    //    complexd det = a11 * a22 - a12 * a21;
    //    if (abs(det) < EPS) throw runtime_error("Sistema singolare (TE).");
    //    complexd Ar1y = (b1 * a22 - a12 * b2) / det;
    //    complexd Ap2y = (a11 * b2 - b1 * a21) / det;

    //    Solution sol;
    //    sol.Ar1 = Vector3D(complexd(0.0), Ar1y, complexd(0.0));
    //    sol.Ap2 = Vector3D(complexd(0.0), Ap2y, complexd(0.0));

    //    // Power calculations
    //    return calculate_power(sol);
    //}
    //





    
    //per ora preferisco risolvere direttamente con le formule chiuse
    Solution solve_TE() {
        // Dati noti
        complexd Ap1y = Ap1.y;
        complexd K1 = m1.k;
        complexd K2 = m2.k;
        complexd mu1 = m1.mu;
        complexd mu2 = m2.mu;

        complexd Kp1z = kp1.z;
        complexd Kp2z = kp2.z;

      
        complexd denom = Kp1z * mu2 + Kp2z * mu1;

       
        // metto Tolleranza relativa per evitare divisioni per numeri troppo piccoli 
        double tol = 1e-14;
        double scale = (double)abs(denom);
        if (scale < 1e-30) scale = 1.0;

        Solution sol;

        if (abs(denom) > tol * scale) {
           
            complexd Ap2y = (complexd(2.0, 0.0) * Kp1z * mu2 / denom) * Ap1y;                    // (3.11)
            complexd Ar1y = ((Kp1z * mu2 - Kp2z * mu1) / denom) * Ap1y;                        // (3.12)

            sol.Ap2 = Vector3D(complexd(0.0), Ap2y, complexd(0.0));
            sol.Ar1 = Vector3D(complexd(0.0), Ar1y, complexd(0.0));
            return calculate_power(sol);
        }
        else {
            //////////////////// PROBABILMENTE NON NECESSARIO //////////////////////
            
            // Fallback numerico stabile: calcola i coefficienti di Fresnel tramite impedenze (
            complexd cos_i = Kp1z / K1;   // cos(theta_i) in termini di k
            complexd cos_t = Kp2z / K2;   // cos(theta_t) in termini di k
            complexd Z1 = m1.Z;
            complexd Z2 = m2.Z;

            complexd den_r = Z2 * cos_i + Z1 * cos_t;
            if (abs(den_r) < 1e-30) {
                // ultima risorsa: evita crash e segnala ritorno (riflessione totale simulata)
                cerr << "solve_TE fallback: denominatore Fresnel quasi zero. Restituisco riflessione totale.\n";
                sol.Ap2 = Vector3D(complexd(0.0), complexd(0.0), complexd(0.0));
                sol.Ar1 = Vector3D(complexd(0.0), Ap1y, complexd(0.0));
                return calculate_power(sol);
            }

            complexd rTE = (Z2 * cos_i - Z1 * cos_t) / den_r;
            complexd tTE = (complexd(2.0) * Z2 * cos_i) / den_r;

            sol.Ar1 = Vector3D(complexd(0.0), rTE * Ap1y, complexd(0.0));
            sol.Ap2 = Vector3D(complexd(0.0), tTE * Ap1y, complexd(0.0));

            // Debug (commenta se non desideri output)
            cerr << "[solve_TE fallback] denom troppo piccolo: |denom|=" << abs(denom)
                << ", using Fresnel rTE=" << rTE << ", tTE=" << tTE << "\n";

            return calculate_power(sol);
        }
    }


 






                    //// FUNZIONANTE ////
                      /// !!!! !!! ///

    // ========================================================================
    // NEW solve_TM IMPLEMENTATION (using Eigen)
    // This function implements the physically correct 4x4 linear system
    // derived from Maxwell's boundary conditions for the TM case.
    //
    // System: M * x = rhs
    // Unknowns: x = [Ar1x, Ar1z, Ap2x, Ap2z]^T
    //
    // Equations:
    // 1. Continuity of E_x:     Ap1x + Ar1x = Ap2x
    // 2. Continuity of H_y:     (1/mu1)(k_p1 x Ap1)_y + (1/mu1)(k_r1 x Ar1)_y = (1/mu2)(k_p2 x Ap2)_y
    // 3. Reflected TEM:         Ar1 . k_r1 = 0
    // 4. Transmitted TEM:       Ap2 . k_p2 = 0
    // ========================================================================
Solution solve_TM() {
    // Known: Ap1.x and Ap1.z 
    complexd Ap1x = Ap1.x;
    complexd Ap1z = Ap1.z;

    // Unknowns vector x = [Ar1x, Ar1z, Ap2x, Ap2z]
    Matrix4cd M = Matrix4cd::Zero();
    Vector4cd rhs = Vector4cd::Zero();

    // Get constants and vectors
    Vector3D kp1v = kp1;
    Vector3D kp2v = kp2;
    complexd mu1 = m1.mu;
    complexd mu2 = m2.mu;

    // --- Build Matrix M and Vector rhs ---

    // Eq 1: Continuity of E_x
    // (1)Ar1x + (0)Ar1z + (-1)Ap2x + (0)Ap2z = -Ap1x
    M(0, 0) = 1.0;
    M(0, 2) = -1.0;
    rhs(0) = -Ap1x;

    // Eq 2: Continuity of H_y
    // H_y = (-j/mu) * (k_z*A_x - k_x*A_z)
    // (1/mu1)(k_p1z*Ap1x - k_p1x*Ap1z) + (1/mu1)(-k_p1z*Ar1x - k_p1x*Ar1z) = (1/mu2)(k_p2z*Ap2x - k_p1x*Ap2z)
    // Rearranged:
    // (-k_p1z/mu1)Ar1x + (-k_p1x/mu1)Ar1z + (-k_p2z/mu2)Ap2x + (k_p1x/mu2)Ap2z = -(1/mu1)(k_p1z*Ap1x - k_p1x*Ap1z)
    M(1, 0) = -kp1v.z / mu1;
    M(1, 1) = -kp1v.x / mu1;
    M(1, 2) = -kp2v.z / mu2;
    M(1, 3) = kp1v.x / mu2;
    rhs(1) = -(kp1v.z * Ap1x - kp1v.x * Ap1z) / mu1;

    // Eq 3: Reflected TEM (Ar1 . k_r1 = 0)
    // k_r1 = (k_p1x, 0, -k_p1z)
    // (k_p1x)Ar1x + (-k_p1z)Ar1z = 0
    M(2, 0) = kp1v.x;
    M(2, 1) = -kp1v.z;
    rhs(2) = 0.0;

    // Eq 4: Transmitted TEM (Ap2 . k_p2 = 0)
    // k_p2 = (k_p1x, 0, k_p2z)
    // (k_p1x)Ap2x + (k_p2z)Ap2z = 0
    M(3, 2) = kp1v.x;
    M(3, 3) = kp2v.z;
    rhs(3) = 0.0;


    // === Solve 4x4 Linear System using Eigen ===
    Eigen::FullPivLU<Matrix4cd> lu(M);

    if (!lu.isInvertible()) {
        throw runtime_error("Sistema lineare TM singolare o mal condizionato (Eigen).");
    }

    Vector4cd solu = lu.solve(rhs);

    // Extract solution
    Solution sol;
    sol.Ar1 = Vector3D(solu(0), complexd(0.0), solu(1));
    sol.Ap2 = Vector3D(solu(2), complexd(0.0), solu(3));

    return calculate_power(sol);
}





/*

///  VERSIONE CON FORMULE CHIUSE ////
Solution solve_TM() {
    // === Costanti e dati noti ===
    complexd Ap1x = Ap1.x;
    complexd Ap1z = Ap1.z;

    complexd K1 = m1.k;
    complexd mu1 = m1.mu;
    complexd mu2 = m2.mu;

    complexd Kp1vx = kp1.x;
    complexd Kp1vz = kp1.z;
    complexd Kp2vz = kp2.z;

    // === Denominatore comune ===
    complexd denom =
        (Kp1vx * Kp1vx) * (K1 * K1) * (Kp1vz * Kp1vz) * mu1 +
        (Kp1vx * Kp1vx) * (Kp2vz * Kp2vz) * K1 * mu2 +
        (Kp1vz * Kp1vz) * (Kp2vz * Kp2vz) * K1 * mu2 +
        (Kp1vz * Kp1vz) * (Kp2vz * Kp2vz) * mu1;

    // === Numeratori  ===

    // A_r1x
    complexd num_Ar1x =
        Kp1vz * (
            -std::pow(Kp1vx, 2.0) * mu1 * std::pow(K1, 2.0) * Ap1x
            - Kp1vx * Kp2vz * K1 * mu2 * Ap1z
            + Kp1vz * Kp2vz * K1 * mu2 * Ap1x
            - std::pow(Kp2vz, 2.0) * mu1 * Ap1x
            );

    complexd Ar1x = num_Ar1x / denom;

    // A_r1z
    complexd num_Ar1z =
        Kp1vx * (
            -std::pow(Kp1vx, 2.0) * std::pow(K1, 2.0) * mu1 * Ap1z
            - Kp1vx * Kp2vz * K1 * mu2 * Ap1z
            + Kp1vz * Kp2vz * K1 * mu2 * Ap1x
            - std::pow(Kp2vz, 2.0) * mu1 * Ap1x
            );

    complexd Ar1z = num_Ar1z / denom;

    // A_p2x
    complexd num_Ap2x =
        Kp2vz * K1 * mu2 * (
            std::pow(Kp1vx, 2.0) * Ap1x
            + 2.0 * std::pow(Kp1vz, 2.0) * Ap1x
            - Kp1vx * Kp1vz * Ap1z
            );

    complexd Ap2x = num_Ap2x / denom;

    // A_p2z
    complexd num_Ap2z =
        Kp1vx * std::pow(K1, 2.0) * mu2 * (
            -std::pow(Kp1vx, 2.0) * Ap1z
            - 2.0 * std::pow(Kp1vz, 2.0) * Ap1z
            + Kp1vx * Kp1vz * Ap1x
            );

    complexd Ap2z = num_Ap2z / denom;

    // === il risultato ===
    Solution sol;
    sol.Ar1 = Vector3D(Ar1x, complexd(0.0), Ar1z);
    sol.Ap2 = Vector3D(Ap2x, complexd(0.0), Ap2z);

    return calculate_power(sol);
}




*/






    Solution calculate_power(const Solution& partial) {
        Solution sol = partial;
        // compute Poynting z component: S_z = 0.5 * Re{ E × H* }_z
        // For amplitude A, H_amp = (1/(omega*mu)) (k × A)
        auto poynting_z = [&](const Vector3D& A, const Vector3D& kvec, const Medium& m)->complexd {
            Vector3D H = (1.0 / complexd(omega)) * ((1.0 / m.mu) * cross(kvec, A));
            Vector3D ExHc = cross(A, Vector3D(conj(H.x), conj(H.y), conj(H.z)));
            return complexd(0.5) * ExHc.z; // may be complex; physical S_z is Re(...)
            };
        complexd Sp = poynting_z(Ap1, kp1, m1);
        complexd Sr = poynting_z(sol.Ar1, kr1, m1);
        complexd St = poynting_z(sol.Ap2, kp2, m2);

        double Sp_z = real(Sp);
        double Sr_z = real(Sr);
        double St_z = real(St);

        // Incident power should be positive in +z direction for the incident wave (depending on sign conv)
        // In our geometry incident is coming toward +z (from z<0), so Sp_z should be positive.
        // Reflected power should be negative (propagating away -> negative z Poynting).
        double R = 0.0, T = 0.0;
        if (abs(Sp_z) > 1e-18) {
            R = fabs(Sr_z) / Sp_z;
            T = St_z / Sp_z;
        }
        else {
            R = NAN; T = NAN;
        }
        double Aabs = 1.0 - R - T;

        sol.R = R;
        sol.T = T;
        sol.A = Aabs;

        // angle of refraction from real parts of kp2
        double kpxr = real(kp2.x);
        double kpzr = real(kp2.z);
        double theta_t = atan2(kpxr, kpzr); // note: atan2(y,x) but we want atan2(kx, kz)
        sol.theta_t_deg = theta_t * 180.0 / PI;

        return sol;
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    try {
        cout << "=== EM Interface Solver (CLI) ===\n";
        cout << "Inserire frequenza f [Hz]: "; double f; cin >> f;
        double omega = 2.0 * PI * f;

        auto read_complex_rel = [&](const string& name)->complexd {
            double re, im;
            cout << "Inserire " << name << " parte reale (valore relativo): "; cin >> re;
            cout << "Inserire " << name << " parte immaginaria (valore relativo): "; cin >> im;
            return complexd(re, im);
            };

        cout << "--- Mezzo 1 ---\n";
        complexd eps1_r = read_complex_rel("eps_r1");
        complexd mu1_r = read_complex_rel("mu_r1");

        cout << "--- Mezzo 2 ---\n";
        complexd eps2_r = read_complex_rel("eps_r2");
        complexd mu2_r = read_complex_rel("mu_r2");

        cout << "Angolo di incidenza theta_i [gradi]: "; double theta_i; cin >> theta_i;
        cout << "Polarizzazione (TE inserire 1, TM inserire 0): "; int pol; cin >> pol;
        bool isTE = (pol == 1);

        Vector3D Ap1;
        if (isTE) {
            cout << "Inserire ampiezza complessa Ap1_y (parte reale): "; double ar, ai; cin >> ar; cout << " (parte immaginaria): "; cin >> ai;
            Ap1 = Vector3D(0.0, complexd(ar, ai), 0.0);
        }
        else {
            
            cout << "Inserire ampiezza totale TM (A_tot) (parte reale): "; double ar, ai; cin >> ar; cout << " (parte immaginaria): "; cin >> ai;
            complexd A_total = complexd(ar, ai);

            // Calcola automaticamente le componenti Ap1x e Ap1z per garantire la condizione TEM
            // A_p1 deve essere ortogonale a k_p1 = (k1*sin_i, 0, k1*cos_i)
            // Il vettore ortogonale nel piano x-z è A_p1_dir = (cos_i, 0, -sin_i)
            double theta_i_rad = theta_i * PI / 180.0;
            complexd Ap1x = A_total * cos(theta_i_rad);
            complexd Ap1z = A_total * (-sin(theta_i_rad));
            Ap1 = Vector3D(Ap1x, complexd(0.0), Ap1z);

            cout << "  -> Calcolata Ap1.x = " << cstr(Ap1x) << "\n";
            cout << "  -> Calcolata Ap1.z = " << cstr(Ap1z) << "\n";
        }





        EMInterfaceSolver solver(eps1_r, mu1_r, eps2_r, mu2_r, omega, theta_i, Ap1, isTE);
        Solution sol = solver.solve();

        cout << "\n=== Risultati ===\n";
        cout << "Ampiezza riflessa Ar1: \n";
        cout << "  Ar1.x = " << cstr(sol.Ar1.x) << "\n";
        cout << "  Ar1.y = " << cstr(sol.Ar1.y) << "\n";
        cout << "  Ar1.z = " << cstr(sol.Ar1.z) << "\n";
        cout << "Ampiezza rifratta Ap2: \n";
        cout << "  Ap2.x = " << cstr(sol.Ap2.x) << "\n";
        cout << "  Ap2.y = " << cstr(sol.Ap2.y) << "\n";
        cout << "  Ap2.z = " << cstr(sol.Ap2.z) << "\n";

        cout << fixed << setprecision(4);
        cout << "\nRiflettanza R = " << sol.R * 100.0 << " %\n";
        cout << "Trasmittanza T = " << sol.T * 100.0 << " %\n";
        cout << "Assorbanza (dissipata) A = " << sol.A * 100.0 << " %\n";
        cout << "Angolo di penetrazione (theta_t) = " << sol.theta_t_deg << " deg\n";

    }
    catch (exception& e) {
        cerr << "Errore: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
