#pragma once
#include "rvlm/core/Constants.hh"
#include "rvlm/core/SolidArray3d.hh"

template <typename valueType>
class YeeGrid {
public:

    using Const = rvlm::core::Constants<valueType>;

    typedef valueType ValueType;
    typedef rvlm::core::SolidArray3d<valueType> ArrayType;

public:
    YeeGrid(int nx, int ny, int nz,
            ValueType dt, ValueType dx, ValueType dy, ValueType dz)
        : delta_t(dt)
        , delta_x(dx)
        , delta_y(dy)
        , delta_z(dz)
        , Hx         (nx, ny, nz,   0)
        , Hy         (nx, ny, nz,   0)
        , Hz         (nx, ny, nz,   0)
        , Ex         (nx, ny, nz,   0)
        , Ey         (nx, ny, nz,   0)
        , Ez         (nx, ny, nz,   0)
        , mu_Hx      (nx, ny, nz,   Const::MU_0())
        , mu_Hy      (nx, ny, nz,   Const::MU_0())
        , mu_Hz      (nx, ny, nz,   Const::MU_0())
        , epsilon_Ex (nx, ny, nz,   Const::EPS_0())
        , epsilon_Ey (nx, ny, nz,   Const::EPS_0())
        , epsilon_Ez (nx, ny, nz,   Const::EPS_0())
        , sigma_Hx   (nx, ny, nz,   0)
        , sigma_Hy   (nx, ny, nz,   0)
        , sigma_Hz   (nx, ny, nz,   0)
        , sigma_Ex   (nx, ny, nz,   0)
        , sigma_Ey   (nx, ny, nz,   0)
        , sigma_Ez   (nx, ny, nz,   0)
        , D_Hx       (nx, ny, nz,   0)
        , D_Hy       (nx, ny, nz,   0)
        , D_Hz       (nx, ny, nz,   0)
        , C_Ex       (nx, ny, nz,   0)
        , C_Ey       (nx, ny, nz,   0)
        , C_Ez       (nx, ny, nz,   0)
        , D_Ex       (nx, ny, nz,   0)
        , D_Ey       (nx, ny, nz,   0)
        , D_Ez       (nx, ny, nz,   0)
    {}

public:
    ValueType delta_x;
    ValueType delta_y;
    ValueType delta_z;
    ValueType delta_t;

    ArrayType Hx;
    ArrayType Hy;
    ArrayType Hz;
    ArrayType Ex;
    ArrayType Ey;
    ArrayType Ez;

    ArrayType mu_Hx;
    ArrayType mu_Hy;
    ArrayType mu_Hz;
    ArrayType epsilon_Ex;
    ArrayType epsilon_Ey;
    ArrayType epsilon_Ez;
    ArrayType sigma_Hx;
    ArrayType sigma_Hy;
    ArrayType sigma_Hz;
    ArrayType sigma_Ex;
    ArrayType sigma_Ey;
    ArrayType sigma_Ez;

    ArrayType D_Hx;
    ArrayType D_Hy;
    ArrayType D_Hz;
    ArrayType C_Ex;
    ArrayType C_Ey;
    ArrayType C_Ez;
    ArrayType D_Ex;
    ArrayType D_Ey;
    ArrayType D_Ez;
};
