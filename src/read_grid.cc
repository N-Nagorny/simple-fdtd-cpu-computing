#include "read_grid.hh"
#include <vector>
#include <netcdf>

using namespace netCDF;

using Array = std::vector<float>;
using SolidArray3d = rvlm::core::SolidArray3d<float>;

static void readGridVar(SolidArray3d& destination, NcVar var,
                   int startx, int starty, int startz) {

    int nx = var.getDim(0).getSize();
    int ny = var.getDim(1).getSize();
    int nz = var.getDim(2).getSize();

    SolidArray3d::CursorType cursor;
    destination.cursorMoveTo(cursor, startx, starty, startz);

    std::vector<size_t> start(3), count(3);
    count[0] = 1;
    count[1] = 1;
    count[2] = nz;

    destination.fill(0.0f);

    std::vector<float> buf(nz);
    for (int ix = 0; ix < nx; ++ix, destination.cursorMoveToNextX(cursor)) {

	SolidArray3d::CursorType yCursor = cursor;
        for (int iy = 0; iy < ny; ++iy, destination.cursorMoveToNextY(yCursor)) {
            start[0] = ix;
            start[1] = iy;
            start[2] = 0;

            var.getVar(start, count, buf.data());

            SolidArray3d::CursorType zCursor = yCursor;
            for (int iz = 0; iz < nz; ++iz, destination.cursorMoveToNextZ(zCursor)) {
                destination.at(zCursor) = buf[iz];
            }
        }
    }
}

YeeGrid readGridData(std::string const& filename, float deltaT, float deltaX, float deltaY, float deltaZ) {
    NcFile file(filename, NcFile::read);
    int nx = file.getDim("grid_nx").getSize();
    int ny = file.getDim("grid_ny").getSize();
    int nz = file.getDim("grid_nz").getSize();

    YeeGrid grid(nx, ny, nz, deltaT, deltaX, deltaY, deltaZ);

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
