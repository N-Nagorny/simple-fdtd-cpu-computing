#ifndef FDTD_FIELD_H_
#define FDTD_FIELD_H_

#include "rvlm/core/SolidArray3d.hh"

struct fdtd_field {

    float *Hx;
    float *Hy;
    float *Hz;

    float *Ex;
    float *Ey;
    float *Ez;

    float *D_Hx;
    float *D_Hy;
    float *D_Hz;

    float *C_Ex;
    float *C_Ey;
    float *C_Ez;

    float *D_Ex;
    float *D_Ey;
    float *D_Ez;

    float *delta_x_Ex;
    float *delta_y_Ex;
    float *delta_z_Ex;
    float *delta_x_Ey;
    float *delta_y_Ey;
    float *delta_z_Ey;
    float *delta_x_Ez;
    float *delta_y_Ez;
    float *delta_z_Ez;
    float *delta_x_Hx;
    float *delta_y_Hx;
    float *delta_z_Hx;
    float *delta_x_Hy;
    float *delta_y_Hy;
    float *delta_z_Hy;
    float *delta_x_Hz;
    float *delta_y_Hz;
    float *delta_z_Hz;

};

#endif
