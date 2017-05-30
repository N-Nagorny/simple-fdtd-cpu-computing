#pragma once
#include "Dimensions.hh"

namespace rvlm {
namespace core {
/**
 * This class acts as a namespace and contains various mathematical and
 * physycal universal constants. The physical one must be stated in the SI unit
 * system with dimention name specified in parentheses.
 * @see http://en.wikipedia.org/wiki/International_System_of_Units
 */
template <typename T>
class Constants {

public:
    /**
     * Euclid pi number (dimensionless quantity).
     * This constant is defines here because C++ standard library does not
     * guarantee to have @c M_PI defined.
     * @see http://stackoverflow.com/questions/1727881
     */
    static T PI() { return (T)3.141592653589793; }

    /**
     * Light's speed in vacuum (meter/second).
     * @see http://en.wikipedia.org/wiki/Speed_of_light
     */
    static Velocity<T> C() { return (T)299792501 * boost::units::si::meter_per_second; }

    /**
     * Electric constant or vacuum permittivity (farad/meter).
     * @see http://en.wikipedia.org/wiki/Vacuum_permittivity
     */
    static Permittivity<T> EPS_0() { return (T)8.8541878187E-12 * (boost::units::si::farad / boost::units::si::meter); }

    /**
     * Magnetic constant or vacuum permeability (Henry/meter).
     * @see http://en.wikipedia.org/wiki/Vacuum_permeability
     */
    static Permeability<T> MU_0() { return PI() * (T)4E-7 * (boost::units::si::henry / boost::units::si::meter); }

    /** @todo Find proper name and write documentation. */
    static Resistance<T> ETA_0() { return MU_0() * C(); }
};

} // namespace core
} // namespace rvlm
