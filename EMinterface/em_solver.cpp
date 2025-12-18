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
const double PI = acos(-1.0);

// --- Strutture di Base ---

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

    // scientific notation with precision to 8
    ss << scientific << setprecision(8);
    ss << real(c);
    if (imag(c) >= 0) ss << "+j" << imag(c);
    else ss << "-j" << abs(imag(c));
    return ss.str();
}

// --- Struct Medium  ---

struct Medium {
    static constexpr double EPS0 = 8.854187817e-12;
    static constexpr double MU0 = 1.2566370614359173e-6;

    complexd eps_r;
    complexd mu_r;
    double sigma;

    complexd eps_eff;
    complexd mu_abs;
    complexd k;
    complexd eta;

    Medium(complexd _eps_r, complexd _mu_r, double _sigma, double omega)
        : eps_r(_eps_r), mu_r(_mu_r), sigma(_sigma) {

        mu_abs = mu_r * MU0;

        // Calcolo permittività efficace: eps_eff = eps_real - j(sigma/omega)
        complexd eps_absolute = eps_r * EPS0;
        complexd conductivity_term = complexd(0, -1.0) * (sigma / omega);
        eps_eff = eps_absolute + conductivity_term;

        // k = omega * sqrt(mu * eps_eff)
        // La radice principale ha parte reale positiva. Per mezzi passivi, Im(k) deve essere <= 0
        // (Assumendo convenzione fasori exp(jwt - kr)) o >=0 se exp(jwt - jkr).
        // Qui usiamo standard physics: k = beta - j*alpha
        k = omega * std::sqrt(mu_abs * eps_eff);

        // eta = sqrt(mu / eps_eff) questa è l'impedenza intrinseca del mezzo
        eta = std::sqrt(mu_abs / eps_eff);
    }
};


constexpr double Medium::EPS0;
constexpr double Medium::MU0;



// --- Struct Soluzione ---

struct Solution {
    Vector3D Ar1;
    Vector3D Ap2;
    double R = 0.0, T = 0.0, A = 0.0;
    double theta_t_rad = 0.0;
};





// --- Solver ---

class EMInterfaceSolver {
public:
    EMInterfaceSolver(complexd eps1_r, complexd mu1_r, double sigma1,
        complexd eps2_r, complexd mu2_r, double sigma2,
        double omega, double theta_i_rad,
        const Vector3D& Ap1_in, bool isTE)
        : m1(eps1_r, mu1_r, sigma1, omega), m2(eps2_r, mu2_r, sigma2, omega),
        omega(omega), theta_i(theta_i_rad), Ap1(Ap1_in), isTE(isTE)
    {
        // Inizializza immediatamente i vettori d'onda per tutti i mezzi
        initialize_wave_vectors();
    }

    Solution solve() {
        cerr << "--- Debug Info ---\n";
        cerr << "Medium 1: k=" << cstr(m1.k) << ", eta=" << cstr(m1.eta) << "\n";
        cerr << "Medium 2: k=" << cstr(m2.k) << ", eta=" << cstr(m2.eta) << "\n";
        cerr << "wavevector kp1 (inc): " << cstr(kp1.x) << ", " << cstr(kp1.y) << ", " << cstr(kp1.z) << "\n";
        cerr << "wavevectorkp2 (tra): " << cstr(kp2.x) << ", " << cstr(kp2.y) << ", " << cstr(kp2.z) << "\n";
        cerr << "versori di Kp1v:" << cstr(kp1vx) << "," <<cstr(kp1vy) <<","<< cstr(kp1vz)<<"\n";
        cerr << "------------------\n";

        if (isTE) return solve_TE();
        else return solve_TM();
    }

private:
    Medium m1, m2;
    double omega;
    double theta_i;
    Vector3D Ap1;
    bool isTE;


    // wavevectors
    Vector3D kp1, kr1, kp2;


    //// versori dei wavevectors
    double kp1vx,kp1vy,kp1vz, kr1vx,kr1vy,kr1vz;              //tanto questi sono compoenti di un versore quindi sono reali li metto double e sposto kp2 sopra










    //// Calcolo componenti k basato su angolo (helper)
    //WaveVector get_k_components(const Medium& m, complexd theta) {
    //    WaveVector kv;
    //    kv.kx = m.k * std::sin(theta);
    //    kv.ky = 0.0;
    //    kv.kz = m.k * std::cos(theta);
    //    return kv;
    //}

 // Calcola i vettori (o versori) d'onda iniziali basandosi sulla geometria descritta
    void initialize_wave_vectors() {
        // Conversione angolo in radianti
        // Usiamo complexd per coerenza, anche se theta_i input è reale         //INVECE NO LO MODIFICO
        double theta_rad = theta_i;
        double sin_th = std::sin(theta_rad);
        double cos_th = std::cos(theta_rad);

        // -----------------------------------------------------------
        // 1. ONDA INCIDENTE (Mezzo 1) -> VERSORE (Unit Vector)
        // Riferimento Immagine: K_p1vx = sin(theta), K_p1vz = cos(theta)
        // -----------------------------------------------------------
        kp1vx = sin_th;      // Componente x del versore
        kp1vy = 0.0;
        kp1vz = cos_th;      // Componente z del versore

        // -----------------------------------------------------------
        // 2. ONDA RIFLESSA (Mezzo 1) -> VERSORE (Unit Vector)
        // Riferimento Immagine: K_r1vx = K_p1vx, K_r1vz = -K_p1vz
        // -----------------------------------------------------------
        kr1vx = sin_th;       // Speculare rispetto a z: x uguale
        kr1vy = 0.0;
        kr1vz = -cos_th;      // Speculare rispetto a z: z opposto

        // -----------------------------------------------------------
        // 3. ONDA TRASMESSA (Mezzo 2) -> VETTORE D'ONDA COMPLETO
        // Riferimento Immagine: K_p2x = K1 * K_p1vx
        // -----------------------------------------------------------

        // Recuperiamo il numero d'onda scalare completo del mezzo 1 e 2
        complexd K1 = m1.k;
        complexd K2 = m2.k;

        // Conservazione della componente tangenziale (Legge di Snell generalizzata)
        // K_p2x deve essere un vettore d'onda completo (dimensione 1/m)
        complexd Kp2x = K1 * kp1vx;


        // Calcolo componente normale K_p2z
        complexd arg_sq = (K2 * K2) - (Kp2x * Kp2x);
        complexd Kp2z = std::sqrt(arg_sq);

        // FIX APPLIED: apllica la condizione fisica Im(kz) <= 0
        // questo assicura che l'onda decada (o che quantomeno non cresca ) nella direzione +z.
        if (imag(Kp2z) > 1e-15) { //small tolerance for floating point safety
            Kp2z = -Kp2z;
        }

        kp2.x = Kp2x;
        kp2.y = 0.0;
        kp2.z = Kp2z;

        // -----------------------------------------------------------
        // 4. ONDA INCIDENTE (Mezzo 1) -> VETTORE D'ONDA COMPLETO
        kp1.x = m1.k * kp1vx;
        kp1.y = 0.0;
        kp1.z = m1.k * kp1vz;
        // -----------------------------------------------------------
        // 5. ONDA RIFLESSA (Mezzo 1) -> VETTORE D'ONDA COMPLETO
        kr1.x = m1.k * kr1vx;
        kr1.y = 0.0;
        kr1.z = m1.k * kr1vz;
        // -----------------------------------------------------------


    }

    // --- Calcolo Potenza (Poynting) ---
    // Generalizzato per vettori complessi
    Solution calculate_power(const Solution& partial) {
        Solution sol = partial;

        // Funzione locale per Sz = 0.5 * Re(E x H*)_z
        auto calc_Sz = [&](const Vector3D& E, const Vector3D& k_vec, const Medium& m) -> double {
            // H = (1 / omega*mu) * (k x E)
            Vector3D H = (1.0 / (m.mu_abs * omega)) * cross(k_vec, E);

            // S = 0.5 * E x conj(H)
            Vector3D H_conj(std::conj(H.x), std::conj(H.y), std::conj(H.z));
            Vector3D S = 0.5 * cross(E, H_conj);
            return std::real(S.z);
            };

        double Sz_inc = calc_Sz(Ap1, kp1, m1);  /// qua invece  vuole i Kp1 completo che devo ancora definire
        double Sz_ref = calc_Sz(sol.Ar1, kr1, m1);
        double Sz_tra = calc_Sz(sol.Ap2, kp2, m2);

        // Sz_ref sarà negativo perché viaggia in -z. Prendiamo il valore assoluto per il rapporto.
        if (std::abs(Sz_inc) > 1e-20) {
            sol.R = std::abs(Sz_ref) / Sz_inc;
            sol.T = Sz_tra / Sz_inc;
            sol.A = 1.0 - sol.R - sol.T;
        }
        else {
            sol.R = 0; sol.T = 0; sol.A = 0;
        }

        // Calcolo angolo rifrazione reale (dalla direzione di fase)
        double theta_t = std::atan2(std::real(kp2.x), std::real(kp2.z));
        sol.theta_t_rad = theta_t ;

        return sol;
    }

  



    // Solve TE case using 2x2 linear system 
    Solution solve_TE() {
        // unknowns: Ar1y, Ap2y
        complexd Ap1y = Ap1.y;


        complexd K1 = m1.k; // Numero d'onda scalare mezzo 1
        complexd mu1 = m1.mu_abs;
        complexd mu2 = m2.mu_abs;

        double _kp1vz = kp1vz;
        double _kr1vz = kr1vz;

        complexd Ar1y = ((K1*_kp1vz*mu2)-(kp2.z*mu1))*Ap1y / ((kp2.z*mu1)-(K1*_kr1vz*mu2));

        complexd Ap2y = Ap1y + Ar1y;

        

        Solution sol;
        sol.Ar1 = Vector3D(complexd(0.0), Ar1y, complexd(0.0));
        sol.Ap2 = Vector3D(complexd(0.0), Ap2y, complexd(0.0));

        return calculate_power(sol);
    }





    Solution solve_TM() {
        // Mappatura variabili locali 
        complexd Ap1x = Ap1.x;
        complexd Ap1z = Ap1.z;

        complexd K1 = m1.k;      
        complexd mu1 = m1.mu_abs;
        complexd mu2 = m2.mu_abs;


        complexd _kp2x = kp2.x; // kp2x e kp2z come nelle formule
        complexd _kp2z = kp2.z;

        // Pre-calcolo quadrati per leggibilità e corrispondenza visiva alle formule
        complexd kp2x_sq = _kp2x * _kp2x;
        complexd kp2z_sq = _kp2z * _kp2z;
        complexd kr1vx_sq = kr1vx * kr1vx;
        complexd kr1vz_sq = kr1vz * kr1vz;


        	// --- Denominatore Ar1x ---
        complexd Denom = (-mu2 * _kp2z * K1 * kr1vx_sq - mu2 * _kp2z * K1 * kr1vz_sq + mu1  * kp2x_sq * kr1vz + kp2z_sq * mu1 * kr1vz);

            // --- Numeratore Ar1x ---
        complexd Num_Ar1x = (kp1vz * mu2 * _kp2z * K1 * kr1vz - mu1 * kp2x_sq * kr1vz - kp2z_sq * mu1 * kr1vz) * Ap1x - (kp1vx * K1 * mu2 * _kp2z * kr1vz) * Ap1z;

        	// --- Calcolo Ar1x ---
        complexd Ar1x = Num_Ar1x / Denom;

        complexd Ar1z = -Ar1x * (kr1vx / kr1vz);

        complexd Ap2x = Ap1x + Ar1x;

        complexd Ap2z = -Ap1x * (_kp2x / _kp2z) - Ar1x * (_kp2x / _kp2z);

        // Verifica denominatore
        if (abs(Denom) < 1e-15) {
            throw runtime_error("Dominatore nullo in solve_TM (configurazione singolare).");
        }


        // --- Costruzione Soluzione ---
        Solution sol;
        sol.Ar1 = Vector3D(Ar1x, complexd(0.0), Ar1z);
        sol.Ap2 = Vector3D(Ap2x, complexd(0.0), Ap2z);

        return calculate_power(sol);
    }
    
    
    
    
    
    
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    try {
        cout << "=== EM Interface Solver (Lossy Media Support) ===\n";
        cout << "Inserire frequenza f [Hz]: "; double f; cin >> f;
        double omega = 2.0 * PI * f;
   
        
        auto read_complex_rel = [&](const string& name)->complexd {
            double re, im;
            cout << "Inserire " << name << " parte reale: "; cin >> re;
            cout << "Inserire " << name << " parte immaginaria: "; cin >> im;
            return complexd(re, im);
            };


        cout << "\n--- Mezzo 1 ---\n";
        complexd eps1_r = read_complex_rel("eps_r1");
        complexd mu1_r = read_complex_rel("mu_r1");
        cout << "Inserire conduttivita sigma1 [S/m]: "; double sigma1; cin >> sigma1;

        cout << "\n--- Mezzo 2 ---\n";
        complexd eps2_r = read_complex_rel("eps_r2");
        complexd mu2_r = read_complex_rel("mu_r2");
        cout << "Inserire conduttivita sigma2 [S/m]: "; double sigma2; cin >> sigma2;

        cout << "\nAngolo di incidenza theta_i [radianti]: "; double theta_i_rad; cin >> theta_i_rad;
        cout << "Polarizzazione (TE = 1, TM = 0): "; int pol; cin >> pol;
        bool isTE = (pol == 1);

        Vector3D Ap1;
        if (isTE) {
            cout << "Inserire ampiezza Ap1_y (Re Im): "; double ar, ai; cin >> ar >> ai;
            Ap1 = Vector3D(0.0, complexd(ar, ai), 0.0);
        }
        else {
            cout << "Inserire ampiezza totale TM (Re Im): "; double ar, ai; cin >> ar >> ai;
            complexd A_total = complexd(ar, ai);
            double th_rad = theta_i_rad;
            // Proiezione su assi x,z per incidenza
            Ap1 = Vector3D(A_total * cos(th_rad), 0.0, A_total * -sin(th_rad));
        }

        EMInterfaceSolver solver(eps1_r, mu1_r, sigma1, eps2_r, mu2_r, sigma2, omega, theta_i_rad, Ap1, isTE);
        Solution sol = solver.solve();




        cout << "\n=== Risultati ===\n";
        cout << "Riflessa Ar1: " << cstr(sol.Ar1.x) << ", " << cstr(sol.Ar1.y) << ", " << cstr(sol.Ar1.z) << "\n";
        cout << "Trasmessa Ap2: " << cstr(sol.Ap2.x) << ", " << cstr(sol.Ap2.y) << ", " << cstr(sol.Ap2.z) << "\n";

        cout << fixed << setprecision(2);
        cout << "\nBilancio Energetico:\n";
        cout << "Riflettanza R: " << sol.R * 100.0 << " %\n";
        cout << "Trasmittanza T: " << sol.T * 100.0 << " %\n";
        cout << "Assorbanza A : " << sol.A * 100.0 << " %\n";
        cout << "Angolo rifrazione reale: " << sol.theta_t_rad << " rad\n";


    }
    catch (const exception& e) {
        cerr << "Errore critico: " << e.what() << endl;
    }
    return 0;
}