#pragma once
#include "rvlm/core/SolidArray3d.hh"

class YeeGrid {

    using Array = std::vector<float>;
    using SolidArray3d = rvlm::core::SolidArray3d<float>;

public:
    YeeGrid(int nx, int ny, int nz, float deltaT, float deltaX, float deltaY, float deltaZ):
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
        delta_x = deltaX;
        delta_y = deltaY;
        delta_z = deltaZ;
    }

    SolidArray3d mu_Hx;
    SolidArray3d mu_Hy;
    SolidArray3d mu_Hz;
    SolidArray3d sigma_Hx;
    SolidArray3d sigma_Hy;
    SolidArray3d sigma_Hz;

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

    float delta_x;
    float delta_y;
    float delta_z;
    float delta_t;

    SolidArray3d Hx;
    SolidArray3d Hy;
    SolidArray3d Hz;

    SolidArray3d Ex;
    SolidArray3d Ey;
    SolidArray3d Ez;
};
