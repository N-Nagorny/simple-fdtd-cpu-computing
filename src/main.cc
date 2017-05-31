#include <iostream>
#include <fstream>
#include <boost/format.hpp>
#include <boost/multiprecision/mpfr.hpp>

#include <boost/units/cmath.hpp>
#include <boost/units/get_dimension.hpp>
#include <boost/units/get_system.hpp>
#include "Dimensions.hh"

#include "Constants.hh"

#include "calc_field.hh"
#include "resistive_source.hh"

typedef boost::multiprecision::number<
    boost::multiprecision::mpfr_float_backend<250,
    boost::multiprecision::allocate_stack> > float100;

namespace rvlm {
namespace core {
    using namespace boost::multiprecision;

    template
    class Constants<float100>;
}
}

typedef rvlm::core::Constants<float> Const;

// power<N>()
// power<N, M>()
// peciprocal()
//

template <typename Unit, int N, int M=1>
using static_power_helper = boost::units::unit<
        typename boost::units::static_power<
            typename boost::units::get_dimension<Unit>::type,
            typename boost::units::static_rational<N, M>::type>::type,
        typename boost::units::get_system<Unit>::type>;

template <int N, typename Unit, typename Y>
auto mypow(boost::units::quantity<Unit, Y> const& x)
      -> boost::units::quantity<static_power_helper<Unit, N>, Y>
{
    using std::pow;
    return decltype(mypow<N>(x))::from_value(pow(x.value(), N));
}

template <int N, int M, typename Unit, typename Y>
auto mypow(boost::units::quantity<Unit, Y> const& x)
      -> boost::units::quantity<static_power_helper<Unit, N, M>, Y>
{
    using std::pow;
    return decltype(mypow<N, M>(x))::from_value(pow(x.value(), (Y)N/M));
}

    const int nx = 129;
    const int ny = 129;
    const int nz = 129;
    const int nl = 20;
    const Length<float> lambda = 0.05f * boost::units::si::meter;
    const boost::units::quantity<boost::units::si::angular_velocity, float> omega  = 2 * Const::PI() * boost::units::si::radian * Const::C() / lambda;
    const Length<float> dx = lambda / (float)nl;
    const Length<float> dy = dx;
    const Length<float> dz = dx;
    // TODO: FIx it.
    const Time<float> dt =
        mypow<-1>(Const::C() * mypow<1,2>(
            1.0f/(dx*dx) + 1.0f/(dy*dy) + 1.0f/(dz*dz))) *
        boost::units::quantity<boost::units::si::dimensionless, float>(0.5f);

template<typename valueType>
void dumpImage(std::string const& filename, YeeGrid<valueType> const& grid) {
    int nx = grid.Ez.getCountX();
    int ny = grid.Ez.getCountY();
    int z0 = grid.Ex.getCountZ() / 2;

    valueType maxEx = 0;
    for (int ix = 0; ix < nx; ++ix)
    for (int iy = 0; iy < ny; ++iy) {
        if ( ((nx/2 - 3) <= ix && ix <= (nx/2 + 3)) &&
             ((ny/2 - 3) <= iy && iy <= (ny/2 + 3)) )
            continue;


        valueType curAbs = abs(grid.Ez.at(ix, iy, z0));
        if (curAbs > maxEx)
            maxEx = curAbs;
    }

    std::cout << "MAX_EZ: " << maxEx << std::endl;

    std::ofstream output(filename, std::ios::binary);
    output << "P6\n" << nx << " " << ny << "\n" << 255 << "\n";

    valueType x = 0 / maxEx;
    if (boost::math::isnan(x)) {
        maxEx = 1;
        std::cout << "No NaN hack" << std::endl;
    }

    valueType two = 2;
    char color[3];
    for (int iy = 0; iy < ny; ++iy)
    for (int ix = 0; ix < nx; ++ix) {
        float val  = grid.Ez.at(ix, iy, z0);
        unsigned blue = 0;
        unsigned red  = 0;
        valueType level = std::abs(val) / maxEx * 255;

        //if (level > 0)
        //    std::cout << "level: " << level << std::endl;

        (val >= 0 ? blue : red) = (int)level;

//        std::cout << "red: " << red << " blue: " << blue << std::endl;
        color[0] = red;   // red
        color[1] = 0;     // green
        color[2] = blue;  // blue
        output.write(&color[0], 3);
    }
}

template<typename valueTypeLo, typename valueTypeHi>
void dumpDifference(
        std::string const& filename,
        YeeGrid<valueTypeLo> const& grid1,
        YeeGrid<valueTypeHi> const& grid2) {

    int nx = grid2.Ez.getCountX();
    int ny = grid2.Ez.getCountY();
    int z0 = grid2.Ex.getCountZ() / 2;

    valueTypeHi maxRelErr   = 0;
    valueTypeHi sumRelErr   = 0;
    valueTypeHi sumRelErrSq = 0;
    for (int ix = 0; ix < nx; ++ix)
    for (int iy = 0; iy < ny; ++iy) {
        valueTypeLo lo       = grid1.Ez.at(ix, iy, z0);
        valueTypeHi hi       = grid2.Ez.at(ix, iy, z0);
        valueTypeHi hiAbs    = (hi > 0 ? hi : -hi);
        valueTypeHi err      = hi - lo;
        valueTypeHi errAbs   = (err > 0 ? err : -err);
        valueTypeHi relErr   = errAbs / hiAbs;
        if (boost::math::isnan(relErr))
            relErr = 0;

        valueTypeHi relErrSq = relErr * relErr;

        sumRelErr   += relErr;
        sumRelErrSq += relErrSq;

        if (relErr > maxRelErr)
            maxRelErr = relErr;
    }

    valueTypeHi meanRelErr   = sumRelErr / (nx*ny);
    valueTypeHi meanRelErrSq = sumRelErrSq / (nx*ny);
    valueTypeHi dispersion   = meanRelErrSq - meanRelErr*meanRelErr;

    std::cout << "maxRelErr:    " << maxRelErr    << std::endl
              << "meanRelErr:   " << meanRelErr   << std::endl
              << "meanRelErrSq: " << meanRelErrSq << std::endl
              << "dispresion:   " << dispersion   << std::endl;

    std::ofstream output(filename, std::ios::binary);
    output << "P6\n" << nx << " " << ny << "\n" << 255 << "\n";

    valueTypeHi x = 0/maxRelErr;
    if (boost::math::isnan(x))
        maxRelErr = 1;

    char color[3];
    for (int iy = 0; iy < ny; ++iy)
    for (int ix = 0; ix < nx; ++ix) {
        valueTypeLo lo  = grid1.Ez.at(ix, iy, z0);
        valueTypeHi hi  = grid2.Ez.at(ix, iy, z0);

        unsigned red  = 0;
        valueTypeHi hiAbs    = (hi > 0 ? hi : -hi);
        valueTypeHi err      = hi - lo;
        valueTypeHi errAbs   = (err > 0 ? err : -err);
        valueTypeHi relErr   = errAbs / hiAbs;
        if (boost::math::isnan(relErr))
            relErr = 0;


        valueTypeHi level = relErr/maxRelErr * 1000;

        //if (level > 0)
        //    std::cout << "level: " << level << std::endl;

        red = (int)level;

        color[0] = red;   // red
        color[1] = 0;     // green
        color[2] = 0;  // blue
        output.write(&color[0], 3);
    }
}

template<typename valueType>
void setupGrid(YeeGrid<valueType>& grid) {
    int x0 = nx / 2;
    int y0 = ny / 2;
    int z0 = nz / 2;

    Permittivity<valueType> epsilonAntenna = rvlm::core::Constants<valueType>::EPS_0() * (valueType)10;
    ElectricConductivity<valueType> sigmaAntenna = ElectricConductivity<valueType>::from_value(10000);
    for (int i = 0; i < nl; ++i) {
        grid.epsilon_Ez.at(x0, y0, z0 + i) = epsilonAntenna;
        grid.epsilon_Ez.at(x0, y0, z0 - i) = epsilonAntenna;
        grid.sigma_Ez.at(x0, y0, z0 + i)   = sigmaAntenna;
        grid.sigma_Ez.at(x0, y0, z0 - i)   = sigmaAntenna;
    }
}

template <typename T>
ElectricPotential<T> voltage(Time<T> time) {
    return std::sin(omega*time) * 1000.0f * boost::units::si::volt;
}

int main(int argc, char *argv[]) {
    using YeeGridF = YeeGrid<float100>;
    using ResistiveSourceF = ResistiveSource<float100>;

    YeeGrid<float> grid1(nx, ny, nz, dt, dx, dy, dz);
    setupGrid(grid1);
    calcCoefs(grid1);

    //YeeGridF grid2(nx, ny, nz, dt, dx, dy, dz);
    //setupGrid(grid2);
    //calcCoefs(grid2);

    int x0 = nx / 2;
    int y0 = ny / 2;
    int z0 = nz / 2;
    //ResistiveSource<float> rsource1(x0, y0, z0, 10 * boost::units::si::ohm);
    //rsource1.calcCoefs(grid1);

    //ResistiveSourceF rsource2(x0, y0, z0, 10 * boost::units::si::ohm);
    //rsource2.calcCoefs(grid2);

    Time<float> time = 0 * boost::units::si::second;
    int iter = 0;
    while (iter < 100) {
        std::cout << "(1) Rocking iteration #" << iter << std::endl;
        calcH(grid1);

        //rsource1.resqueFields(grid1);
        calcE(grid1);
        //rsource1.updateFields(grid1, voltage<float>(time));

        //std::cout << "(2) Rocking iteration #" << iter << std::endl;
        //calcH(grid2);

        //rsource2.resqueFields(grid2);
        //calcE(grid2);
        //rsource2.updateFields(grid2, voltage<float100>(time));

        std::string fn = str(boost::format("fdtd_diff_%02d.ppm") % iter);
        //dumpDifference(fn, grid1, grid2);

        iter++;
        time = dt*(float)iter;
    }

    return 0;
}
