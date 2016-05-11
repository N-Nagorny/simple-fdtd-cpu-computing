#include "read_grid.hh"
#include "calc_field.hh"
#include <iostream>

int main(int argc, char *argv[]) {
    YeeGrid grid = readGridData("test.nc", 60e-12);
    for (int ix = 0; ix < grid.mu_Hx.getCountX(); ix++)
    for (int iy = 0; iy < grid.mu_Hx.getCountY(); iy++)
    for (int iz = 0; iz < grid.mu_Hx.getCountZ(); iz++) {
        std::cout << grid.mu_Hx.at(ix, iy, iz) << std::endl;
    }
    calcCoefs(grid);
    for (int ix = 0; ix < grid.D_Hx.getCountX(); ix++)
    for (int iy = 0; iy < grid.D_Hx.getCountY(); iy++)
    for (int iz = 0; iz < grid.D_Hx.getCountZ(); iz++) {
        std::cout << grid.D_Hx.at(ix, iy, iz) << std::endl;
    }

    calcH(grid);
    calcE(grid);
    return 0;
}
