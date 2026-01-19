	#pragma once
#include <vector>
#include <array>
#include "include/utils/constants.hpp"

#include <cmath>

#include <stdexcept>
#include "../geom/vec2.hpp"

namespace em {

// Mesh 2D su disco: nodi (x,y) + triangoli (indici)
struct Mesh2D {
    std::vector<geom::Vec2> nodes;
    std::vector<std::array<int,3>> tris;
};

// Mesh polare:
// - n_pr >=2 : numero punti lungo ciascun raggio (incluso centro e bordo)
// - n_rays >=3 : numero raggi/spicchi angolari
inline Mesh2D generate_disk_mesh(double Rf, int n_pr, int n_rays) {
    if (Rf <= 0.0) throw std::runtime_error("generate_disk_mesh: Rf must be > 0");
    if (n_pr < 2) throw std::runtime_error("generate_disk_mesh: n_pr must be >= 2");
    if (n_rays < 3) throw std::runtime_error("generate_disk_mesh: n_rays must be >= 3");

    Mesh2D m;

    // Nodo 0: centro
    m.nodes.push_back({0.0, 0.0});

    // Anelli 1..(n_pr-1)
    for(int ir=1; ir<=n_pr-1; ++ir) {
        double r = (double)ir / (double)(n_pr-1) * Rf;
        for(int k=0; k<n_rays; ++k) {
            double a = 2.0* utils::pi
                * (double)k/(double)n_rays;
            m.nodes.push_back({r*std::cos(a), r*std::sin(a)});
        }
    }

    auto idx = [&](int ir, int k)->int {
        // ir=0 => centro
        if(ir==0) return 0;
        // ir>=1: offset 1 + (ir-1)*n_rays + k
        return 1 + (ir-1)*n_rays + (k % n_rays);
    };

    // Triangoli tra centro e primo anello
    for(int k=0;k<n_rays;++k){
        m.tris.push_back({ idx(0,0), idx(1,k), idx(1,k+1) });
    }

    // Triangoli tra anelli consecutivi
    for(int ir=1; ir<=n_pr-2; ++ir){
        for(int k=0;k<n_rays;++k){
            int a = idx(ir, k);
            int b = idx(ir, k+1);
            int c = idx(ir+1, k);
            int d = idx(ir+1, k+1);

            // due triangoli per quad
            m.tris.push_back({a, c, d});
            m.tris.push_back({a, d, b});
        }
    }

    return m;
}

} // namespace em
