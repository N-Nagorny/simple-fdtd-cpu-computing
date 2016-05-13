#include "read_grid.hh"
#include "calc_field.hh"
#include "resistive_source.hh"
#include <iostream>
#include <fstream>
#include <boost/format.hpp>

const float c = 299792458.0f;
const float pi = 3.14159265358;

const int nx = 100;
const int ny = 100;
const int nz = 100;
const int nl = 20;
const float lambda = 0.05;
const float omega  = 2 * pi * c / lambda;
const float dx = lambda / nl;
const float dy = dx;
const float dz = dx;
const float dt = 1/(c * std::sqrt(1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz)));

void dumpImage(std::string const& filename, YeeGrid& grid) {
    int nx = grid.Ez.getCountX();
    int ny = grid.Ez.getCountY();
    int z0 = grid.Ex.getCountZ() / 2;

    float maxEx = 0;
    for (int ix = 0; ix < nx; ++ix)
    for (int iy = 0; iy < ny; ++iy) {
        if ( ((nx/2 - 3) <= ix && ix <= (nx/2 + 3)) &&
             ((ny/2 - 3) <= iy && iy <= (nx/2 + 3)) )
            continue;

        float curAbs = std::abs(grid.Ez.at(ix, iy, z0));
        if (curAbs > maxEx)
            maxEx = curAbs;
    }

    std::cout << "MAX: " << maxEx << std::endl;

    std::ofstream output(filename, std::ios::binary);
    output << "P6\n" << nx << " " << ny << "\n" << 255 << "\n";

    unsigned char color[3];
    for (int iy = 0; iy < ny; ++iy)
    for (int ix = 0; ix < nx; ++ix) {
        float val  = grid.Ez.at(ix, iy, z0);
        unsigned char  blue = 0, red = 0;
        if (val >= 0)
            blue = std::log(val/maxEx + 1) / std::log(2) * 255;
        else
            red = std::log(-val/maxEx + 1) / std::log(2) * 255;


        color[0] = red;        /* red */
        color[1] = 0;        /* green */
        color[2] = blue;  /* blue */
        output.write(reinterpret_cast<char*>(&color[0]), 3);
    }
}

YeeGrid setupGrid() {

    YeeGrid grid(100, 100, 100, dt, dx, dy, dz);
    int x0 = nx / 2;
    int y0 = ny / 2;
    int z0 = nz / 2;

    float epsilonAntenna = 10;
    float sigmaAntenna = 10000000;
    for (int i = 0; i < nl; ++i) {
        grid.epsilon_Ez.at(x0, y0, z0 + i) = epsilonAntenna;
        grid.epsilon_Ez.at(x0, y0, z0 + i) = epsilonAntenna;
        grid.sigma_Ez.at(x0, y0, z0 - i) = sigmaAntenna;
        grid.sigma_Ez.at(x0, y0, z0 - i) = sigmaAntenna;
    }

    return grid;
}

float voltage(float time) {
    return std::sin(omega*time);
}

int main(int argc, char *argv[]) {
    YeeGrid grid = setupGrid();
    calcCoefs(grid);

    dumpImage("test.ppm", grid);

    int x0 = nx / 2;
    int y0 = ny / 2;
    int z0 = nz / 2;
    ResistiveSource rsource(x0, y0, z0, 10);

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
