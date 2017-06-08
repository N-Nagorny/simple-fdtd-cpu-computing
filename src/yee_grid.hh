#pragma once
#include "rvlm/core/SolidArray3d.hh"
#include "Constants.hh"
#include "Dimensions.hh"

template <typename Y>
class YeeGrid {
public:

    using Const = Constants<Y>;
    using ValueType = Y;

    template <typename U>
    using ArrayType = rvlm::core::SolidArray3d<U>;

public:
    YeeGrid(int nx, int ny, int nz,
            Time<Y> dt, Length<Y> dx, Length<Y> dy, Length<Y> dz)
        : delta_t(dt)
        , delta_x(dx)
        , delta_y(dy)
        , delta_z(dz)
        , Hx         (nx, ny, nz,   MagneticIntensity<Y>::from_value(0))
        , Hy         (nx, ny, nz,   MagneticIntensity<Y>::from_value(0))
        , Hz         (nx, ny, nz,   MagneticIntensity<Y>::from_value(0))
        , Ex         (nx, ny, nz,   ElectricIntensity<Y>::from_value(0))
        , Ey         (nx, ny, nz,   ElectricIntensity<Y>::from_value(0))
        , Ez         (nx, ny, nz,   ElectricIntensity<Y>::from_value(0))
        , mu_Hx      (nx, ny, nz,   Const::MU_0())
        , mu_Hy      (nx, ny, nz,   Const::MU_0())
        , mu_Hz      (nx, ny, nz,   Const::MU_0())
        , epsilon_Ex (nx, ny, nz,   Const::EPS_0())
        , epsilon_Ey (nx, ny, nz,   Const::EPS_0())
        , epsilon_Ez (nx, ny, nz,   Const::EPS_0())
        , sigma_Hx   (nx, ny, nz,   MagneticLoss<Y>::from_value(0))
        , sigma_Hy   (nx, ny, nz,   MagneticLoss<Y>::from_value(0))
        , sigma_Hz   (nx, ny, nz,   MagneticLoss<Y>::from_value(0))
        , sigma_Ex   (nx, ny, nz,   ElectricConductivity<Y>::from_value(0))
        , sigma_Ey   (nx, ny, nz,   ElectricConductivity<Y>::from_value(0))
        , sigma_Ez   (nx, ny, nz,   ElectricConductivity<Y>::from_value(0))
        , D_Hx       (nx, ny, nz,   MagneticCurlCoefficient<Y>::from_value(0))
        , D_Hy       (nx, ny, nz,   MagneticCurlCoefficient<Y>::from_value(0))
        , D_Hz       (nx, ny, nz,   MagneticCurlCoefficient<Y>::from_value(0))
        , C_Ex       (nx, ny, nz,   Dimensionless<Y>(0))
        , C_Ey       (nx, ny, nz,   Dimensionless<Y>(0))
        , C_Ez       (nx, ny, nz,   Dimensionless<Y>(0))
        , D_Ex       (nx, ny, nz,   ElectricCurlCoefficient<Y>::from_value(0))
        , D_Ey       (nx, ny, nz,   ElectricCurlCoefficient<Y>::from_value(0))
        , D_Ez       (nx, ny, nz,   ElectricCurlCoefficient<Y>::from_value(0))
    {}

public:
    Length<ValueType>
        delta_x,
        delta_y,
        delta_z;

    Time<ValueType>
        delta_t;

    ArrayType<MagneticIntensity<Y>> Hx, Hy, Hz;
    ArrayType<ElectricIntensity<Y>> Ex, Ey, Ez;

    ArrayType<Permeability<Y>> mu_Hx, mu_Hy, mu_Hz;
    ArrayType<Permittivity<Y>> epsilon_Ex, epsilon_Ey, epsilon_Ez;
    ArrayType<ElectricConductivity<Y>> sigma_Ex, sigma_Ey, sigma_Ez;
    ArrayType<MagneticLoss<Y>> sigma_Hx, sigma_Hy, sigma_Hz;
    ArrayType<MagneticCurlCoefficient<Y>> D_Hx, D_Hy, D_Hz;
    ArrayType<ElectricCurlCoefficient<Y>> D_Ex, D_Ey, D_Ez;

    ArrayType<Dimensionless<Y>> C_Ex, C_Ey, C_Ez;
};
