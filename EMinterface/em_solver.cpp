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

using namespace std;
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

    


    /*
    // Solve TE case using 2x2 linear system 
    Solution solve_TE() {
        // unknowns: Ar1y, Ap2y
        complexd Ap1y = Ap1.y;
        complexd kp1z = kp1.z;
        complexd kr1z = kr1.z;
        complexd kp2z = kp2.z;
        complexd mu1 = m1.mu;
        complexd mu2 = m2.mu;

        // Equations:
        // 1) Ap1y + Ar1y - Ap2y = 0
        // 2) (kp1z*Ap1y + kr1z*Ar1y)/mu1 - (kp2z*Ap2y)/mu2 = 0

        // matrix form M * [Ar1y, Ap2y]^T = rhs
        // from 1) Ar1y - Ap2y = -Ap1y
        // from 2) (kr1z/mu1) * Ar1y + (-kp2z/mu2) * Ap2y = -kp1z*Ap1y / mu1

        complexd a11 = complexd(1.0);               // Ar1y coef from eq1
        complexd a12 = complexd(-1.0);              // Ap2y coef from eq1
        complexd b1 = -Ap1y;

        complexd a21 = kr1z / mu1;
        complexd a22 = -kp2z / mu2;
        complexd b2 = -(kp1z * Ap1y) / mu1;

        // Solve 2x2
        complexd det = a11 * a22 - a12 * a21;
        if (abs(det) < EPS) throw runtime_error("Sistema singolare (TE).");
        complexd Ar1y = (b1 * a22 - a12 * b2) / det;
        complexd Ap2y = (a11 * b2 - b1 * a21) / det;

        Solution sol;
        sol.Ar1 = Vector3D(complexd(0.0), Ar1y, complexd(0.0));
        sol.Ap2 = Vector3D(complexd(0.0), Ap2y, complexd(0.0));

        // Power calculations
        return calculate_power(sol);
    }
    */






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






    /*
    
    // posso risolvere il sistema lineare 4x4 direttamente con eliminazione gaussiana e trovare così subito le incognite (Ar1x, Ar1z, Ap2x, Ap2z) 
    Solution solve_TM() {
        // Known: Ap1.x and Ap1.z 

        complexd Ap1x = Ap1.x;
        complexd Ap1z = Ap1.z;

        // vettore di incognite x = [Ar1x, Ar1z, Ap2x, Ap2z]
        const int N = 4;
        vector<vector<complexd>> M(N, vector<complexd>(N, complexd(0.0)));
        vector<complexd> rhs(N, complexd(0.0));

        // equazioni che uso per costruire sistema:
        // (I)  ẑ × (Ap1 + Ar1 - Ap2) = 0  -> two scalar eqs (x and y components)
        // (II) (1/mu1) ẑ × (kp1×Ap1 + kr1×Ar1) - (1/mu2) ẑ × (kp2×Ap2) = 0 -> two scalar eqs (x and y)
        // Additionally Ar1 must satisfy Ar1·kr1 = 0 (orthogonality) -> 1 eq
        
        // a questo punto assumo che l'onda sia già TEM, quindi Ap1·kp1 = 0

        // But the vector eqs produce the necessary independent equations. 
        // Eq0: x-component of ẑ×(...) = 0
        // Eq1: y-component of ẑ×(...) = 0
        // Eq2: x-component of (1/mu1) ẑ×(kp1×Ap1 + kr1×Ar1) - (1/mu2) ẑ×(kp2×Ap2) = 0
        // Eq3: y-component of same as Eq2

        // Unknown vectors:
        // Ar1 = (u0, 0, u1)
        // Ap2 = (u2, 0, u3)

        // Precompute some constant contributions
        Vector3D Ap1v = Ap1;
        Vector3D Ar1v; // symbolic
        Vector3D Ap2v;

        // Helper lambda: compute zhat x V for symbolic coefficients on unknowns
        // For each equation, we find coefficients for u0..u3 and known term from Ap1.

        auto zcross = [](const Vector3D& v)->Vector3D {
            // zhat = (0,0,1) cross v = (-v.y, v.x, 0)
            return Vector3D(-v.y, v.x, complexd(0.0));
            };

        // Define symbolic Ar1 = [u0, 0, u1], Ap2 = [u2, 0, u3]
        // We compute each lhs as linear combination: lhs = C0*u0 + C1*u1 + C2*u2 + C3*u3 + known_term

        // Precompute known part for Eq0 and Eq1 from z×Ap1:
        Vector3D zAp1 = zcross(Ap1v); // known
        // z×Ar1 contributes: -Ar1.y, Ar1.x, 0 -> Ar1.y = 0 => (0, Ar1x,0)
        // z×Ar1 = (0, u0, 0)
        // z×Ap2 = (0, Ap2x, 0) but with minus sign since Ap1+Ar1-Ap2 = 0

        // Eq0: x-component of z×(...) = 0 -> zAp1.x + zAr1.x - zAp2.x = 0
        // the only nontrivial component from z×E continuity is the y-component:
        // EqY1: y-comp: zAp1.y + zAr1.y - zAp2.y = 0 => Ap1x + Ar1x - Ap2x = 0

   
        
        //Riga 0: Continuità di E_y
        //cioè ho Ap1x+Ar1x - Ap2x = 0 che riarrangio per il sistema Mx = b e diventa 
        // (1)Ar1x + (0)Ar1z + (-1)Ap2x + (0)Ap2z = -Ap1x

        M[0][0] = complexd(1.0); // u0 cioè il coefficente di Ar1x
        M[0][1] = complexd(0.0); // u1    // Coeff. per u1 (Ar1z)
        M[0][2] = complexd(-1.0); // u2     cioè il coefficente di Ap2x
        M[0][3] = complexd(0.0); // u3     // Coeff. per u3 (Ap2z)

        rhs[0] = -Ap1x;  // questo è termine noto

        // Next build Eq1 and Eq2 from H continuity vector (cioe devo calolare i termini z × (k×A) )
        // Compute known term: z×(kp1×Ap1) (this is known)
        Vector3D kp1v = kp1;
        Vector3D kr1v = kr1;
        Vector3D kp2v = kp2;

        Vector3D kpxAp1 = cross(kp1v, Ap1v);
        Vector3D z_kpxAp1 = zcross(kpxAp1);

        // For Ar1: term z×(kr1 × Ar1) is linear in Ar1 components u0,u1
        // Let Ar1 = [u0, 0, u1]
        // compute kr1 × Ar1 = cross(kr1, Ar1) symbolically: cross(k, [u0,0,u1])
        // cross = ( ky*u1 - kz*0, kz*u0 - kx*u1, kx*0 - ky*u0 ) but ky=0
        // So kr1×Ar1 = (0 - kr1.z*0, kr1.z*u0 - kr1.x*u1, kr1.x*0 - 0*u0) => (0, kr1.z*u0 - kr1.x*u1, 0)
        // Then z×(kr1×Ar1) = (- (kr1.z*u0 - kr1.x*u1).y , (kr1.z*u0 - kr1.x*u1).x, 0) but it's simpler to compute directly with symbolic coefs:
        // Actually kr1×Ar1 has only y-component = (kr1.z*u0 - kr1.x*u1)
        // z×(vector with only y-component Y) = (-Y, 0, 0)

        // So z×(kr1×Ar1) = ( -(kr1.z*u0 - kr1.x*u1), 0, 0 )

        // Similarly for kp2×Ap2 with Ap2 = [u2,0,u3]:
        // kp2×Ap2 has y-component = kp2.z*u2 - kp2.x*u3
        // z×(kp2×Ap2) = ( -(kp2.z*u2 - kp2.x*u3), 0, 0 )

        // H continuity: (1/mu1)[ z×(kp1×Ap1) + z×(kr1×Ar1) ] - (1/mu2)[ z×(kp2×Ap2) ] = 0
        // This gives components: x and y (but as above y are zeros); main non-trivial component is x
        // We'll form two scalar equations: x-component and y-component of that vector equation.
        // From the algebra above, z×(kr1×Ar1) contributes only to x-component.
        // z×(kp1×Ap1) may have both x and y components.

        // Compute known parts numeric:
        Vector3D z_kp1xAp1 = z_kpxAp1;
        // z×(kr1×Ar1) contributions coefficients:
        complexd kr1_z = kr1v.z;
        complexd kr1_x = kr1v.x;
        complexd kp2_z = kp2v.z;
        complexd kp2_x = kp2v.x;

        // Eq1: x-component of H-continuity:
        // (1/mu1) * ( z_kp1xAp1.x + ( - (kr1.z*u0 - kr1.x*u1) ) ) - (1/mu2) * ( - (kp2.z*u2 - kp2.x*u3) ) = 0
        // Rearranged as coefficients for u0..u3 and rhs = - (1/mu1) * z_kp1xAp1.x

        complexd mu1 = m1.mu;
        complexd mu2 = m2.mu;

        // coefficient u0:
        M[1][0] = (-kr1_z) / mu1; // because inside there's -(kr1.z*u0 - kr1.x*u1) => term for u0 is -kr1.z, divided by mu1
        // coefficient u1:
        M[1][1] = (kr1_x) / mu1;  // because -(kr1.z*u0 - kr1.x*u1) gives +kr1.x * u1
        // coefficient u2:
        M[1][2] = (-(-kp2_z)) / mu2; // from - (1/mu2)* ( - (kp2.z*u2 - kp2.x*u3) ) => +kp2.z/mu2
        // coefficient u3:
        M[1][3] = (-(-(-kp2_x))) / mu2; // be careful: simplifies to (-(-(-kp2_x)))/mu2 => -( -kp2_x)/mu2 = kp2_x/mu2?
        // Let's compute directly: - (1/mu2) * ( - (kp2.z*u2 - kp2.x*u3) ) = (1/mu2)*(kp2.z*u2 - kp2.x*u3)
        // So coef u2 = kp2.z / mu2, coef u3 = -kp2.x / mu2
        M[1][2] = kp2_z / mu2;
        M[1][3] = -kp2_x / mu2;

        rhs[1] = -(z_kp1xAp1.x) / mu1;

        // Eq2: y-component of H-continuity:
        // compute z×(kp1×Ap1).y + z×(kr1×Ar1).y - (mu1/mu2) * z×(kp2×Ap2).y = 0
        // But from symbolic earlier, z×(kr1×Ar1).y = 0; z×(kp2×Ap2).y = 0
        // Thus Eq2 reduces to (1/mu1) * z_kp1xAp1.y = 0 -> gives value maybe nonzero depending on Ap1
        // So Eq2 becomes simply that known term must be canceled by nothing -> this is an equation for unknowns only if z×(kr1×Ar1).y or z×(kp2×Ap2).y nonzero
        // Let's compute contributions more generally using full vector algebra to be safe.

        // We'll compute coefficients for Eq2 by explicitly computing symbolic expressions using the generic cross formulas.

        // Build linear coefficients by constructing matrices of influence:
        // For a given unknown u = u0..u3, compute z×(k×A) contribution when A has that single unknown =1 and others 0.

        auto compute_z_k_cross_A_coeffs = [&](const Vector3D& k, int unknownIndex) -> Vector3D {
            // unknownIndex: 0->Ar1x, 1->Ar1z, 2->Ap2x, 3->Ap2z
            Vector3D A(0.0, 0.0, 0.0);
            if (unknownIndex == 0) A.x = complexd(1.0), A.y = complexd(0.0), A.z = complexd(0.0); // Ar1x
            if (unknownIndex == 1) A.x = complexd(0.0), A.y = complexd(0.0), A.z = complexd(1.0); // Ar1z
            if (unknownIndex == 2) A.x = complexd(1.0), A.y = complexd(0.0), A.z = complexd(0.0); // Ap2x
            if (unknownIndex == 3) A.x = complexd(0.0), A.y = complexd(0.0), A.z = complexd(1.0); // Ap2z
            Vector3D kxA = cross(k, A);
            return zcross(kxA);
            };

        // For Eq2 and Eq3 we'll build general coefficients using above lambda:
        // Eq index 2: x-component of (1/mu1)[ z×(kp1×Ap1) + z×(kr1×Ar1) ] - (1/mu2)[ z×(kp2×Ap2) ] = 0
        // Eq index 3: y-component

        // Known part:
        Vector3D knownH = z_kp1xAp1; // z×(kp1×Ap1)

        // For Ar1 unknowns (0,1) using kr1:
        for (int ui = 0; ui < 2; ++ui) {
            Vector3D coeff = compute_z_k_cross_A_coeffs(kr1v, ui);
            // add (1/mu1) * coeff to matrix row 2 (x) and row3 (y), for u0,u1 respectively
            M[2][ui] += coeff.x / mu1;
            M[3][ui] += coeff.y / mu1;
        }
        // For Ap2 unknowns (2,3) using kp2 with minus (1/mu2) factor
        for (int ui = 2; ui < 4; ++ui) {
            Vector3D coeff = compute_z_k_cross_A_coeffs(kp2v, ui);
            M[2][ui] += -coeff.x / mu2;
            M[3][ui] += -coeff.y / mu2;
        }
        // RHS are negative of known part divided by mu1
        rhs[2] = -knownH.x / mu1;
        rhs[3] = -knownH.y / mu1;

        // NOTE: earlier we already set row 1 (index 0) as continuity y-component of E
        // Now we must ensure rows 2..3 include contributions from Ar1 and Ap2. But we also must ensure unknown Ar1 orthogonality Ar1·kr1 = 0 is satisfied.
        // However the system we built (with H-continuity x & y) plus E-continuity y and the orthogonality Ar1·kr1 maybe redundant/consistent.
        // To ensure a full-rank 4x4 system we will replace row 3 by orthogonality constraint Ar1·kr1 = 0:
        // Ar1x*kr1x + Ar1z*kr1z = 0  => coefficients: [kr1x, kr1z, 0,0] , rhs = 0

        // la riga 3 è la condizione di ortogonalità di Ar1 con kr1 (Ar1·kr1 = 0) che si traduce in  Ar1x*kr1x + Ar1z*kr1z = 0
        M[3][0] = kr1v.x;   // coefficient for Ar1x
        M[3][1] = kr1v.z;
        M[3][2] = complexd(0.0);
        M[3][3] = complexd(0.0);
        rhs[3] = complexd(0.0);

        // Now solve 4x4 linear system using Gaussian elimination (partial pivot)
        // Convert to augmented matrix
        vector<vector<complexd>> A(N, vector<complexd>(N + 1, complexd(0.0)));
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) A[i][j] = M[i][j];
            A[i][N] = rhs[i];
        }

        // Gaussian elimination
        for (int i = 0; i < N; i++) {
            // pivot
            int piv = i;
            for (int r = i + 1; r < N; r++) {
                if (abs(A[r][i]) > abs(A[piv][i])) piv = r;
            }
            if (abs(A[piv][i]) < 1e-15) throw runtime_error("Sistema lineare TM singolare o mal condizionato.");
            if (piv != i) swap(A[piv], A[i]);
            complexd diag = A[i][i];
            for (int c = i; c <= N; c++) A[i][c] /= diag;
            for (int r = 0; r < N; r++) {
                if (r == i) continue;
                complexd fac = A[r][i];
                if (abs(fac) < 1e-18) continue;
                for (int c = i; c <= N; c++) A[r][c] -= fac * A[i][c];
            }
        }
        vector<complexd> solu(N);

        // adesso a questo punto la matrice è in forma ridotta e posso leggere le soluzioni
        // Le prime N colonne (indici 0 a N-1) formano la matrice identità I
        //L'ultima colonna (la N-esima) contiene i valori della soluzione del sistema Ax=b
        //L'ultimo passaggio solu[i] = A[i][N] semplicemente copia questa colonna finale nel vettore solu.

        for (int i = 0; i < N; i++) solu[i] = A[i][N];

        Solution sol;
        sol.Ar1 = Vector3D(solu[0], complexd(0.0), solu[1]);
        sol.Ap2 = Vector3D(solu[2], complexd(0.0), solu[3]);

        return calculate_power(sol);
    }
    */


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

    // === Numeratori (formule dal PDF) ===

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
            cout << "Inserire ampiezza complessa Ap1_x (parte reale): "; double arx, aix; cin >> arx; cout << " (parte immaginaria): "; cin >> aix;
            cout << "Inserire ampiezza complessa Ap1_z (parte reale): "; double arz, aiz; cin >> arz; cout << " (parte immaginaria): "; cin >> aiz;
            Ap1 = Vector3D(complexd(arx, aix), complexd(0.0), complexd(arz, aiz));
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
