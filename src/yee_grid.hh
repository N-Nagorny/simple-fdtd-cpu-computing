#pragma once
#include <vector>
#include "rvlm/core/SolidArray3d.hh"

class YeeGrid {

    using Array = std::vector<float>;
    using SolidArray3d = rvlm::core::SolidArray3d<float>;

public:

    YeeGrid(int nx, int ny, int nz):
        x_Hx(nx),                   x_Ex(nx),
        y_Hx(ny),                   y_Ex(ny),
        z_Hx(nz),                   z_Ex(nz),
        x_Hy(nx),                   x_Ey(nx),
        y_Hy(ny),                   y_Ey(ny),
        z_Hy(nz),                   z_Ey(nz),
        x_Hz(nx),                   x_Ez(nx),
        y_Hz(ny),                   y_Ez(ny),
        z_Hz(nz),                   z_Ez(nz),
        mu_Hx(nx+1, ny+1, nz+1),          epsilon_Ex(nx+1, ny+1, nz+1),
        mu_Hy(nx+1, ny+1, nz+1),          epsilon_Ey(nx+1, ny+1, nz+1),
        mu_Hz(nx+1, ny+1, nz+1),          epsilon_Ez(nx+1, ny+1, nz+1),
        sigma_Hx(nx+1, ny+1, nz+1),       sigma_Ex(nx+1, ny+1, nz+1),
        sigma_Hy(nx+1, ny+1, nz+1),       sigma_Ey(nx+1, ny+1, nz+1),
        sigma_Hz(nx+1, ny+1, nz+1),       sigma_Ez(nx+1, ny+1, nz+1)
    {}

    // Magnetic field arays
    Array x_Hx;
    Array x_Hy;
    Array x_Hz;
    Array y_Hx;
    Array y_Hy;
    Array y_Hz;
    Array z_Hx;
    Array z_Hy;
    Array z_Hz;

    SolidArray3d mu_Hx;
    SolidArray3d mu_Hy;
    SolidArray3d mu_Hz;
    SolidArray3d sigma_Hx;
    SolidArray3d sigma_Hy;
    SolidArray3d sigma_Hz;

    // Electric field arays
    Array x_Ex;
    Array x_Ey;
    Array x_Ez;
    Array y_Ex;
    Array y_Ey;
    Array y_Ez;
    Array z_Ex;
    Array z_Ey;
    Array z_Ez;

    SolidArray3d epsilon_Ex;
    SolidArray3d epsilon_Ey;
    SolidArray3d epsilon_Ez;
    SolidArray3d sigma_Ex;
    SolidArray3d sigma_Ey;
    SolidArray3d sigma_Ez;
};
