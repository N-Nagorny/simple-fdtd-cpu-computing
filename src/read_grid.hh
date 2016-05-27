#pragma once
#include <vector>
#include <string>
#include <netcdf>
#include "yee_grid.hh"

template <typename valueType>
static void readGridVar(rvlm::core::SolidArray3d<valueType>& dest,
                        netCDF::NcVar var,
                        int startx, int starty, int startz) {

    int nx = var.getDim(0).getSize();
    int ny = var.getDim(1).getSize();
    int nz = var.getDim(2).getSize();

    auto cursor = dest.getCursor(startx, starty, startz);

    std::vector<size_t> start(3), count(3);
    count[0] = 1;
    count[1] = 1;
    count[2] = nz;

    dest.fill(0.0f);

    std::vector<valueType> buf(nz);
    for (int ix = 0; ix < nx; ++ix, dest.cursorMoveToNextX(cursor)) {
        auto yCursor = cursor;
        for (int iy = 0; iy < ny; ++iy, dest.cursorMoveToNextY(yCursor)) {
            start[0] = ix;
            start[1] = iy;
            start[2] = 0;

            var.getVar(start, count, buf.data());

            auto zCursor = yCursor;
            for (int iz = 0; iz < nz; ++iz, dest.cursorMoveToNextZ(zCursor)) {
                dest.at(zCursor) = buf[iz];
            }
        }
    }
}

template <typename valueType>
YeeGrid<valueType> readGridData(std::string const& filename,
                                valueType deltaT,
                                valueType deltaX,
                                valueType deltaY,
                                valueType deltaZ) {

    netCDF::NcFile file(filename, netCDF::NcFile::read);
    int nx = file.getDim("grid_nx").getSize();
    int ny = file.getDim("grid_ny").getSize();
    int nz = file.getDim("grid_nz").getSize();

    YeeGrid<valueType> grid(nx, ny, nz, deltaT, deltaX, deltaY, deltaZ);

    readGridVar(grid.epsilon_Ex, file.getVar("world_epsilon_Ex"),  0, 0, 0);
    readGridVar(grid.epsilon_Ey, file.getVar("world_epsilon_Ey"),  0, 0, 0);
    readGridVar(grid.epsilon_Ez, file.getVar("world_epsilon_Ez"),  0, 0, 0);
    readGridVar(grid.sigma_Ex,   file.getVar("world_sigma_Ex"),    0, 0, 0);
    readGridVar(grid.sigma_Ey,   file.getVar("world_sigma_Ey"),    0, 0, 0);
    readGridVar(grid.sigma_Ez,   file.getVar("world_sigma_Ez"),    0, 0, 0);
    readGridVar(grid.mu_Hx,      file.getVar("world_mu_Hx"),       1, 1, 1);
    readGridVar(grid.mu_Hy,      file.getVar("world_mu_Hy"),       1, 1, 1);
    readGridVar(grid.mu_Hz,      file.getVar("world_mu_Hz"),       1, 1, 1);
    readGridVar(grid.sigma_Hx,   file.getVar("world_sigma_Hx"),    1, 1, 1);
    readGridVar(grid.sigma_Hy,   file.getVar("world_sigma_Hy"),    1, 1, 1);
    readGridVar(grid.sigma_Hz,   file.getVar("world_sigma_Hz"),    1, 1, 1);

    return grid;
}
