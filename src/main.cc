#include <iomanip>
#include <iostream>
#include <fstream>
#include <boost/format.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/io.hpp>
#include "Dimensions.hh"

#include "Constants.hh"
#include "BoostUnitsHelpers.hh"

#include "calc_field.hh"
#include "CurrentSource.hh"
#include "ConvolutionalPML.hh"

template <typename Y>
class Problem {
public:

    using Const = Constants<Y>;

    const int nx;
    const int ny;
    const int nz;
    const int nl;
    const Length<Y> lambda;
    const AngularVelocity<Y> omega;
    const Length<Y> dx;
    const Length<Y> dy;
    const Length<Y> dz;
    const Time<Y> dt;

    Problem()
        : nx(129)
        , ny(129)
        , nz(129)
        , nl(20)
        , lambda(Y(0.05) * boost::units::si::meter)
        , omega(Dimensionless<Y>(2)*Const::PI() * Const::C()/lambda * boost::units::si::radian)
        , dx(lambda / Y(nl))
        , dy(dx)
        , dz(dx)
        , dt(Dimensionless<Y>(0.5)/(Const::C() * usqrt(underlying_cast<Y>(Dimensionless<Y>(1)/(dx*dx) + Dimensionless<Y>(1)/(dy*dy) + Dimensionless<Y>(1)/(dz*dz))))) {}

    void setupGrid(YeeGrid<Y>& grid) {
        int x0 = nx / 2;
        int y0 = ny / 2;
        int z0 = nz / 2;

        Permittivity<Y> epsilonAntenna = Const::EPS_0() * Y(10);
        ElectricConductivity<Y> sigmaAntenna = ElectricConductivity<Y>::from_value(10000);
        for (int i = 0; i < nl; ++i) {
            grid.epsilon_Ez.at(x0, y0, z0 + i) = epsilonAntenna;
            grid.epsilon_Ez.at(x0, y0, z0 - i) = epsilonAntenna;
            grid.sigma_Ez.at(x0, y0, z0 + i)   = sigmaAntenna;
            grid.sigma_Ez.at(x0, y0, z0 - i)   = sigmaAntenna;
        }
    }

    Current<Y> signal(Time<Y> const& time) {
        return usin(underlying_cast<Y>(omega*time)) * Dimensionless<Y>(1000) * boost::units::si::ampere;
    }

    void main() {
        YeeGrid<Y> grid1(nx, ny, nz, dt, dx, dy, dz);
        setupGrid(grid1);
        calcCoefs(grid1);

//        ConvolutionalPML<Y> pml(&grid1,
//                make_triple(HalfOpenIndexRange(0, 8),
//                            HalfOpenIndexRange(0, ny),
//                            HalfOpenIndexRange(0, nz)),
//                                make_triple<Optional<AxialDirection>>(AxialDirection::positive, {}, {}),
//                make_triple<Y>  (0.001, 0, 0),
//                make_triple<Y>  (2.5,   0, 0));

        ConvolutionalPML<Y> pml2 {&grid1,
                { {0, 8},
                  {0, ny},
                  {0, nz}},
                {AxialDirection::negative, {}, {}},
                {Y(0.001), Y(0.001), Y(0.0)},
                {Y(2.5),   Y(2.5),   Y(0.0)}};

        pml2.setup();

        ConvolutionalPML<Y> pml3 {&grid1,
                { {0, nx},
                  {0, 8},
                  {0, nz}},
                {Optional<AxialDirection>(), AxialDirection::negative, {}},
                {Y(0.001), Y(0.001), Y(0.0)},
                {Y(2.5),   Y(2.5),   Y(0.0)}};

        pml3.setup();


        dumpImage("sigma_Ex.ppm", grid1.sigma_Ex);
        dumpImage("sigma_Ey.ppm", grid1.sigma_Ey);

        int x0 = nx / 2;
        int y0 = ny / 2;
        int z0 = nz / 2;
        for (Index iy = 0; iy < 30; ++iy) {
            std::cout << iy << '\t'
                      << grid1.sigma_Ey.at(0, iy, 0) << '\n';
        }
        CurrentSource<Y> source(x0, y0, z0);
        source.calcCoefs(grid1);

        Time<Y> time;
        time = Dimensionless<Y>(0) * boost::units::si::second;
        int iter = 0;
        while (iter < 200) {
            auto fmt = boost::format("field_%00d.ppm") % iter;
            dumpImage(fmt.str(), grid1.Ez);
            calcH(grid1);

            calcE(grid1);
            source.updateFields(grid1, signal(time));

            std::cout << iter << '\t'
                      << time << '\t'
                      << grid1.Ez.at(x0 + 30, y0, z0) << std::endl;

            iter++;
            time = Y(iter) * dt;
        }
    }

    template <typename Q>
    void dumpImage(std::string const& filename, Array<Q> const& field) {
        int nx = field.getCountX();
        int ny = field.getCountY();
        int z0 = field.getCountZ() / 2;

        using boost::units::abs;

        Q maxField = Q::from_value(Y(0));
        for (int ix = 0; ix < nx; ++ix)
        for (int iy = 0; iy < ny; ++iy) {
            if ( ((nx/2 - 3) <= ix && ix <= (nx/2 + 3)) &&
                 ((ny/2 - 3) <= iy && iy <= (ny/2 + 3)) )
                continue;

            Q curAbs = abs(field.at(ix, iy, z0));
            if (curAbs > maxField)
                maxField = curAbs;
        }

        std::ofstream output(filename, std::ios::binary);
        output << "P6\n" << nx << " " << ny << "\n" << 255 << "\n";

//        valueType x = 0 / maxEx;
//        if (boost::units::isnan(x)) {
//            maxEz = 1;
//            std::cout << "No NaN hack" << std::endl;
//        }

        Dimensionless<Y> two {2};
        Q zero = maxField - maxField;

        char color[3];
        for (int iy = 0; iy < ny; ++iy)
        for (int ix = 0; ix < nx; ++ix) {
            Q val  = field.at(ix, iy, z0);

            unsigned blue = 0;
            unsigned red  = 0;
            Dimensionless<Y> level = abs(val) / maxField;
            level *= Y(255);

            //if (level > 0)
            //    std::cout << "level: " << level << std::endl;

            int lvl = (int)level.value();
            (val >= zero ? blue : red) = lvl;

    //        std::cout << "red: " << red << " blue: " << blue << std::endl;
            color[0] = red;   // red
            color[1] = 0;     // green
            color[2] = blue;  // blue
            output.write(&color[0], 3);
        }
    }

};

int main(int argc, char *argv[]) {
    if (argc != 2)
        return 1;

    std::string mode = argv[1];
    if (mode == "float") {
        Problem<float> problem;
        problem.main();
        return 0;
    }

    if (mode == "float100") {
        using float100 = boost::multiprecision::number<
                boost::multiprecision::mpfr_float_backend<100,
                        boost::multiprecision::allocate_stack>>;

        Problem<float100> problem;
        problem.main();
        return 0;
    }

    if (mode == "float16") {
        using float16 = boost::multiprecision::number<
                boost::multiprecision::mpfr_float_backend<16,
                        boost::multiprecision::allocate_stack>>;

        Problem<float16> problem;
        problem.main();
        return 0;
    }

    return 1;
}

#if 0

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

#endif
