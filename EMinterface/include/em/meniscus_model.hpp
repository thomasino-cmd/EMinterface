#pragma once
#include <cmath>
#include <stdexcept>
#include "include/utils/constants.hpp"


namespace em {

    // Parametri (tutti in SI):
    // V   [m^3] volume soluzione (positivo)
    // Dm  [m]   profondità menisco (può essere <0 o >0)
    // Rin [m]   raggio interno cuvetta (positivo)
    struct MeniscusParams {
        double V{ 0.0 };
        double Dm{ 0.0 };
        double Rin{ 0.0 };
    };

    class ParaboloidMeniscus {
    public:
        explicit ParaboloidMeniscus(const MeniscusParams& p) : p_(p) {
            if (p_.V <= 0.0) throw std::runtime_error("Meniscus: V must be > 0");
            if (p_.Rin <= 0.0) throw std::runtime_error("Meniscus: Rin must be > 0");

            // a_pm = V/(pi Rin^2) - Dm/2
            apm_ = (p_.V / (utils::pi
                * p_.Rin * p_.Rin)) - (p_.Dm * 0.5);

            // coeff paraboloide: k = Dm / Rin^2  [1/m]
            k_ = p_.Dm / (p_.Rin * p_.Rin);
        }

        double Rin() const { return p_.Rin; }
        double Dm()  const { return p_.Dm; }
        double V()   const { return p_.V; }

        double apm() const { return apm_; }
        double k()   const { return k_; }

        // Superficie aria-salina: z1 = f2(x1,y1) = k*(x1^2 + y1^2) + a_pm
        double f2(double x1, double y1) const {
            return k_ * (x1 * x1 + y1 * y1) + apm_;
        }

        // Interfaccia piatta equivalente (stesso volume) quando Dm=0:
        // z_flat = V/(pi Rin^2) = a_pm + Dm/2
        double z_flat() const {
            return (p_.V / (utils::pi
                * p_.Rin * p_.Rin));
        }

        bool inside_cuvette(double x1, double y1, double eps = 0.0) const {
            return (x1 * x1 + y1 * y1) <= (p_.Rin + eps) * (p_.Rin + eps);
        }

    private:
        MeniscusParams p_;
        double apm_{ 0.0 };
        double k_{ 0.0 };
    };

} // namespace em
