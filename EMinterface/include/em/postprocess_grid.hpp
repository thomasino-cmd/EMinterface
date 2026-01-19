#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <array>
#include "../geom/vec2.hpp"

namespace em {

    struct Grid2D {
        double x_min{ 0 }, x_max{ 0 };
        double y_min{ 0 }, y_max{ 0 };
        int Nx{ 0 }, Ny{ 0 };
        std::vector<double> val; // Nx*Ny

        Grid2D(double xmin, double xmax, double ymin, double ymax, int nx, int ny)
            : x_min(xmin), x_max(xmax), y_min(ymin), y_max(ymax), Nx(nx), Ny(ny) {

            if (Nx <= 0 || Ny <= 0) throw std::runtime_error("Grid2D: Nx,Ny must be > 0");
            if (!(x_max > x_min) || !(y_max > y_min)) throw std::runtime_error("Grid2D: invalid bounds");
            val.assign((size_t)Nx * (size_t)Ny, 0.0);
        }

        double dx() const { return (x_max - x_min) / (double)Nx; }
        double dy() const { return (y_max - y_min) / (double)Ny; }

        double& at(int ix, int iy) { return val[(size_t)iy * (size_t)Nx + (size_t)ix]; }
        double  at(int ix, int iy) const { return val[(size_t)iy * (size_t)Nx + (size_t)ix]; }

        static bool point_in_tri(const geom::Vec2& p,
            const geom::Vec2& a,
            const geom::Vec2& b,
            const geom::Vec2& c,
            double eps = 1e-14) {
            // barycentric sign method
            auto sgn = [&](const geom::Vec2& p1, const geom::Vec2& p2, const geom::Vec2& p3) {
                return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
                };
            double d1 = sgn(p, a, b);
            double d2 = sgn(p, b, c);
            double d3 = sgn(p, c, a);
            bool has_neg = (d1 < -eps) || (d2 < -eps) || (d3 < -eps);
            bool has_pos = (d1 > eps) || (d2 > eps) || (d3 > eps);
            return !(has_neg && has_pos);
        }

        void add_triangle_uniform(const std::array<geom::Vec2, 3>& t, double add_value) {
            // bounding box -> pixel range
            double xmin = std::min({ t[0].x,t[1].x,t[2].x });
            double xmax = std::max({ t[0].x,t[1].x,t[2].x });
            double ymin = std::min({ t[0].y,t[1].y,t[2].y });
            double ymax = std::max({ t[0].y,t[1].y,t[2].y });

            int ix0 = (int)std::floor((xmin - x_min) / dx()); ix0 = std::max(ix0, 0);
            int ix1 = (int)std::floor((xmax - x_min) / dx()); ix1 = std::min(ix1, Nx - 1);
            int iy0 = (int)std::floor((ymin - y_min) / dy()); iy0 = std::max(iy0, 0);
            int iy1 = (int)std::floor((ymax - y_min) / dy()); iy1 = std::min(iy1, Ny - 1);

            for (int iy = iy0; iy <= iy1; ++iy) {
                for (int ix = ix0; ix <= ix1; ++ix) {
                    double xc = x_min + (ix + 0.5) * dx();
                    double yc = y_min + (iy + 0.5) * dy();
                    if (point_in_tri({ xc,yc }, t[0], t[1], t[2])) {
                        at(ix, iy) += add_value;
                    }
                }
            }
        }

        void save_csv(const std::string& path) const {
            std::ofstream f(path);
            if (!f) throw std::runtime_error("Grid2D: cannot open output file");

            f << "x,y,value\n";
            for (int iy = 0; iy < Ny; ++iy) {
                for (int ix = 0; ix < Nx; ++ix) {
                    double xc = x_min + (ix + 0.5) * dx();
                    double yc = y_min + (iy + 0.5) * dy();
                    f << xc << "," << yc << "," << at(ix, iy) << "\n";
                }
            }
        }
    };

} // namespace em
