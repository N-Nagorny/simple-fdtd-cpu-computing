#pragma once
#include <vector>
#include "rvlm/core/SolidArray3d.hh"

class YeeGrid {

    using Array = std::vector<float>;
    using SolidArray3d = rvlm::core::SolidArray3d<float>;

public:
    YeeGrid(int nx, int ny, int nz, float deltaT):
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
        sigma_Hz(nx+1, ny+1, nz+1),       sigma_Ez(nx+1, ny+1, nz+1),

        D_Hx(nx + 1, ny + 1, nz + 1),
        D_Hy(nx + 1, ny + 1, nz + 1),
        D_Hz(nx + 1, ny + 1, nz + 1),

        C_Ex(nx + 1, ny + 1, nz + 1),
        C_Ey(nx + 1, ny + 1, nz + 1),
        C_Ez(nx + 1, ny + 1, nz + 1),

        D_Ex(nx + 1, ny + 1, nz + 1),
        D_Ey(nx + 1, ny + 1, nz + 1),
        D_Ez(nx + 1, ny + 1, nz + 1),

        //delta_x_Hx(nx),
        delta_y_Hx(ny),
        delta_z_Hx(nz),
        delta_x_Hy(nx),
        //delta_y_Hy(ny),
        delta_z_Hy(nz),
        delta_x_Hz(nx),
        delta_y_Hz(ny),
        //delta_z_Hz(nz),

        //delta_x_Ex(nx),
        delta_y_Ex(ny),
        delta_z_Ex(nz),
        delta_x_Ey(nx),
        //delta_y_Ey(ny),
        delta_z_Ey(nz),
        delta_x_Ez(nx),
        delta_y_Ez(ny),
        //delta_z_Ez(nz),

        Hx(nx + 1, ny + 1, nz + 1),
        Hy(nx + 1, ny + 1, nz + 1),
        Hz(nx + 1, ny + 1, nz + 1),

        Ex(nx + 1, ny + 1, nz + 1),
        Ey(nx + 1, ny + 1, nz + 1),
        Ez(nx + 1, ny + 1, nz + 1)

    {
        mu_Hx.fill(0.0f);
        mu_Hy.fill(0.0f);
        mu_Hz.fill(0.0f);
        sigma_Hx.fill(0.0f);
        sigma_Hy.fill(0.0f);
        sigma_Hz.fill(0.0f);

        epsilon_Ex.fill(0.0f);
        epsilon_Ey.fill(0.0f);
        epsilon_Ez.fill(0.0f);
        sigma_Ex.fill(0.0f);
        sigma_Ey.fill(0.0f);
        sigma_Ez.fill(0.0f);

        D_Hx.fill(0.0f);
        D_Hy.fill(0.0f);
        D_Hz.fill(0.0f);

        C_Ex.fill(0.0f);
        C_Ey.fill(0.0f);
        C_Ez.fill(0.0f);

        D_Ex.fill(0.0f);
        D_Ey.fill(0.0f);
        D_Ez.fill(0.0f);

        Hx.fill(0.0f);
        Hy.fill(0.0f);
        Hz.fill(0.0f);

        Ex.fill(0.0f);
        Ey.fill(0.0f);
        Ez.fill(0.0f);

        delta_t = deltaT;
    }

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

    //Calculation vars
    SolidArray3d D_Hx;
    SolidArray3d D_Hy;
    SolidArray3d D_Hz;

    SolidArray3d C_Ex;
    SolidArray3d C_Ey;
    SolidArray3d C_Ez;

    SolidArray3d D_Ex;
    SolidArray3d D_Ey;
    SolidArray3d D_Ez;

    //Array delta_x_Hx;
    Array delta_y_Hx;
    Array delta_z_Hx;
    Array delta_x_Hy;
    //Array delta_y_Hy;
    Array delta_z_Hy;
    Array delta_x_Hz;
    Array delta_y_Hz;
    //Array delta_z_Hz;

    //Array delta_x_Ex;
    Array delta_y_Ex;
    Array delta_z_Ex;
    Array delta_x_Ey;
    //Array delta_y_Ey;
    Array delta_z_Ey;
    Array delta_x_Ez;
    Array delta_y_Ez;
    //Array delta_z_Ez;

    float delta_t;

    SolidArray3d Hx;
    SolidArray3d Hy;
    SolidArray3d Hz;

    SolidArray3d Ex;
    SolidArray3d Ey;
    SolidArray3d Ez;
};
