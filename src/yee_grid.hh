#pragma once
#include <vector>
#include "rvlm/core/SolidArray3d.hh"

const float eps0 = 8.8541878E-12;
const float mu0  = 1.2566371E-6;

class YeeGrid {

    using Array = std::vector<float>;
    using SolidArray3d = rvlm::core::SolidArray3d<float>;

public:
    YeeGrid(int nx, int ny, int nz, float deltaT, float deltaX, float deltaY, float deltaZ):
        mu_Hx(nx, ny, nz),          epsilon_Ex(nx, ny, nz),
        mu_Hy(nx, ny, nz),          epsilon_Ey(nx, ny, nz),
        mu_Hz(nx, ny, nz),          epsilon_Ez(nx, ny, nz),
        sigma_Hx(nx, ny, nz),       sigma_Ex(nx, ny, nz),
        sigma_Hy(nx, ny, nz),       sigma_Ey(nx, ny, nz),
        sigma_Hz(nx, ny, nz),       sigma_Ez(nx, ny, nz),

        D_Hx(nx, ny, nz),
        D_Hy(nx, ny, nz),
        D_Hz(nx, ny, nz),

        C_Ex(nx, ny, nz),
        C_Ey(nx, ny, nz),
        C_Ez(nx, ny, nz),
        D_Ex(nx, ny, nz),
        D_Ey(nx, ny, nz),
        D_Ez(nx, ny, nz),

        Hx(nx, ny, nz),
        Hy(nx, ny, nz),
        Hz(nx, ny, nz),

        Ex(nx, ny, nz),
        Ey(nx, ny, nz),
        Ez(nx, ny, nz)

    {
        mu_Hx.fill(mu0);
        mu_Hy.fill(mu0);
        mu_Hz.fill(mu0);
        sigma_Hx.fill(0.0f);
        sigma_Hy.fill(0.0f);
        sigma_Hz.fill(0.0f);

        epsilon_Ex.fill(eps0);
        epsilon_Ey.fill(eps0);
        epsilon_Ez.fill(eps0);
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
