#include "read_grid.hh"
#include "calc_field.hh"
#include "resistive_source.hh"
#include <iostream>
#include <fstream>
#include <boost/format.hpp>

const float c = 299792458.0f;
const float pi = 3.14159265358;

const int nx = 129;
const int ny = 129;
const int nz = 129;
const int nl = 20;
const float lambda = 0.05;
const float omega  = 2 * pi * c / lambda;
const float dx = lambda / nl;
const float dy = dx;
const float dz = dx;
const float dt = 0.5/(c * std::sqrt(1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz)));

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

        valueType curAbs = std::abs(grid.Ez.at(ix, iy, z0));
        if (curAbs > maxEx)
            maxEx = curAbs;
    }

    std::ofstream output(filename, std::ios::binary);
    output << "P6\n" << nx << " " << ny << "\n" << 255 << "\n";

    char color[3];
    for (int iy = 0; iy < ny; ++iy)
    for (int ix = 0; ix < nx; ++ix) {
        float val  = grid.Ez.at(ix, iy, z0);
        unsigned blue = 0;
        unsigned red  = 0;
        if (val >= 0)
            blue = std::log(val/maxEx + 1) / std::log(2) * 255;
        else
            red = std::log(-val/maxEx + 1) / std::log(2) * 255;

        color[0] = red;   // red
        color[1] = 0;     // green
        color[2] = blue;  // blue
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
    using YeeGridF = YeeGrid<float>;
    using ResistiveSourceF = ResistiveSource<float>;

    YeeGridF grid(nx, ny, nz, dt, dx, dy, dz);
    setupGrid(grid);
    calcCoefs(grid);

    dumpImage("test.ppm", grid);

    int x0 = nx / 2;
    int y0 = ny / 2;
    int z0 = nz / 2;
    ResistiveSourceF rsource(x0, y0, z0, 10);
    rsource.calcCoefs(grid);

    float time = 0;
    int iter = 0;
    while (true) {
        std::cout << "Rocking iteration #" << iter << std::endl;
        calcH(grid);

        rsource.resqueFields(grid);
        calcE(grid);
        rsource.updateFields(grid, voltage(time));

        std::string filename = str(boost::format("field_%04d.ppm") % iter);
        dumpImage(filename, grid);
        iter++;
        time = dt*iter;
    }

    return 0;
}
