#include <iostream>
#include <fstream>
#include <boost/format.hpp>
#include <boost/multiprecision/mpfr.hpp>

#include "rvlm/core/Constants.hh"

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

    const int nx = 129;
    const int ny = 129;
    const int nz = 129;
    const int nl = 20;
    const float lambda = 0.05;
    const float omega  = 2 * Const::PI() * Const::C() / lambda;
    const float dx = lambda / nl;
    const float dy = dx;
    const float dz = dx;
    const float dt = 0.5/(Const::C() * std::sqrt(1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz)));

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

    valueType epsilonAntenna = 10;
    valueType sigmaAntenna   = 10000;
    for (int i = 0; i < nl; ++i) {
        grid.epsilon_Ez.at(x0, y0, z0 + i) = epsilonAntenna;
        grid.epsilon_Ez.at(x0, y0, z0 - i) = epsilonAntenna;
        grid.sigma_Ez.at(x0, y0, z0 + i)   = sigmaAntenna;
        grid.sigma_Ez.at(x0, y0, z0 - i)   = sigmaAntenna;
    }
}

float voltage(float time) {
    return std::sin(omega*time) * 1000;
}

int main(int argc, char *argv[]) {
    using YeeGridF = YeeGrid<float100>;
    using ResistiveSourceF = ResistiveSource<float100>;

    YeeGrid<float> grid1(nx, ny, nz, dt, dx, dy, dz);
    setupGrid(grid1);
    calcCoefs(grid1);

    YeeGridF grid2(nx, ny, nz, dt, dx, dy, dz);
    setupGrid(grid2);
    calcCoefs(grid2);

    int x0 = nx / 2;
    int y0 = ny / 2;
    int z0 = nz / 2;
    ResistiveSource<float> rsource1(x0, y0, z0, 10);
    rsource1.calcCoefs(grid1);

    ResistiveSourceF rsource2(x0, y0, z0, 10);
    rsource2.calcCoefs(grid2);

    float time = 0;
    int iter = 0;
    while (true) {
        std::cout << "(1) Rocking iteration #" << iter << std::endl;
        calcH(grid1);

        rsource1.resqueFields(grid1);
        calcE(grid1);
        rsource1.updateFields(grid1, voltage(time));

        std::cout << "(2) Rocking iteration #" << iter << std::endl;
        calcH(grid2);

        rsource2.resqueFields(grid2);
        calcE(grid2);
        rsource2.updateFields(grid2, voltage(time));

        std::string fn = str(boost::format("fdtd_diff_%02d.ppm") % iter);
        dumpDifference(fn, grid1, grid2);

        iter++;
        time = dt*iter;
    }

    return 0;
}
