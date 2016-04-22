#ifndef YEE_GRID_H_
#define YEE_GRID_H_

#include <stddef.h>

struct yee_grid {
    size_t nx_Ex;
    size_t nx_Ey;
    size_t nx_Ez;
    size_t ny_Ex;
    size_t ny_Ey;
    size_t ny_Ez;
    size_t nz_Ex;
    size_t nz_Ey;
    size_t nz_Ez;
    size_t nx_Hx;
    size_t nx_Hy;
    size_t nx_Hz;
    size_t ny_Hx;
    size_t ny_Hy;
    size_t ny_Hz;
    size_t nz_Hx;
    size_t nz_Hy;
    size_t nz_Hz;

    // Electric field arays
    float *x_Ex;
    float *x_Ey;
    float *x_Ez;
    float *y_Ex;
    float *y_Ey;
    float *y_Ez;
    float *z_Ex;
    float *z_Ey;
    float *z_Ez;

    float *epsilon_Ex;
    float *epsilon_Ey;
    float *epsilon_Ez;
    float *sigmaE_Ex;
    float *sigmaE_Ey;
    float *sigmaE_Ez;

    // Magnetic field arays
    float *x_Hx;
    float *x_Hy;
    float *x_Hz;
    float *y_Hx;
    float *y_Hy;
    float *y_Hz;
    float *z_Hx;
    float *z_Hy;
    float *z_Hz;

    float *epsilon_Hx;
    float *epsilon_Hy;
    float *epsilon_Hz;
    float *sigmaH_Hx;
    float *sigmaH_Hy;
    float *sigmaH_Hz;
};

#endif
