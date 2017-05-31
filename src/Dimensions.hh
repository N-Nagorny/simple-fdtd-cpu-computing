#pragma once

#include <boost/units/quantity.hpp>
#include <boost/units/physical_dimensions.hpp>
#include <boost/units/systems/si.hpp>

template <typename T>
using Dimensionless =
    boost::units::quantity<boost::units::si::dimensionless, T>;

template <typename T>
using AngularVelocity =
    boost::units::quantity<boost::units::si::angular_velocity, T>;

template <typename T>
using Length = boost::units::quantity<boost::units::si::length, T>;

template <typename T>
using Time = boost::units::quantity<boost::units::si::time, T>;

template <typename T>
using Velocity = boost::units::quantity<boost::units::si::velocity, T>;

template <typename T>
using Permittivity = boost::units::quantity<boost::units::si::permittivity, T>;

template <typename T>
using Permeability = boost::units::quantity<boost::units::si::permeability, T>;

template <typename T>
using Resistance = boost::units::quantity<boost::units::si::resistance, T>;

using MagneticLossUnit = boost::units::divide_typeof_helper<
        boost::units::si::resistance,
        boost::units::si::length>::type;

template <typename T>
using MagneticLoss = boost::units::quantity<MagneticLossUnit, T>;

template <typename T>
using ElectricConductivity = boost::units::quantity<boost::units::si::conductivity, T>;

template <typename T>
using MagneticIntensity = boost::units::quantity<boost::units::si::magnetic_field_intensity, T>;

template <typename T>
using ElectricPotential = boost::units::quantity<boost::units::si::electric_potential, T>;

template <typename T>
using ElectricIntensity =
    boost::units::quantity<
        boost::units::unit<
            boost::units::derived_dimension<
                boost::units::length_base_dimension,   1,
                boost::units::mass_base_dimension,     1,
                boost::units::current_base_dimension, -1,
                boost::units::time_base_dimension,    -3>::type,
            boost::units::si::system>,
        T>;

using ElectricCurlCoefficientUnit =
    boost::units::divide_typeof_helper<
        boost::units::si::time,
        boost::units::si::permittivity>::type;

using MagneticCurlCoefficientUnit =
    boost::units::divide_typeof_helper<
        boost::units::si::time,
        boost::units::si::permeability>::type;

template <typename T>
using ElectricCurlCoefficient = boost::units::quantity<ElectricCurlCoefficientUnit, T>;

template <typename T>
using MagneticCurlCoefficient = boost::units::quantity<MagneticCurlCoefficientUnit, T>;
