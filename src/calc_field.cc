#include "yee_grid.hh"

using Array = std::vector<float>;
using SolidArray3d = rvlm::core::SolidArray3d<float>;

void calcD(SolidArray3d& D, float deltaT,
           const SolidArray3d& perm, const SolidArray3d& sigma) {

    for (int ix = 0; ix < D.getCountX(); ix++)
    for (int iy = 0; iy < D.getCountY(); iy++)
    for (int iz = 0; iz < D.getCountZ(); iz++) {
        float& curD     = D.at(ix, iy, iz);
        float  curPerm  = perm.at(ix, iy, iz);
        float  curSigma = sigma.at(ix, iy, iz);
        float  subexpr  = deltaT / curPerm;
        curD = subexpr / (1 + curSigma * subexpr/2);
    }
}

void calcC(SolidArray3d& C, float deltaT,
           const SolidArray3d& perm, const SolidArray3d& sigma) {

    for (int ix = 0; ix < C.getCountX(); ix++)
    for (int iy = 0; iy < C.getCountY(); iy++)
    for (int iz = 0; iz < C.getCountZ(); iz++) {
        float& curC     = C.at(ix, iy, iz);
        float  curPerm  = perm.at(ix, iy, iz);
        float  curSigma = sigma.at(ix, iy, iz);
        float  subexpr  = curSigma * deltaT / (2 * curPerm);
        curC = (1 - subexpr) / (1 + subexpr);
    }
}

void calcDelta(Array& delta, Array const& coord) {
    size_t n = coord.size();

    delta.resize(n);
    for (int i = 0; i < n - 1; i++)
        delta[i] = coord[i+1] - coord[i];

    // TODO: Allow explicit passing last delta values.
    delta[n-1] = delta[n-2];
}

void calcCoefs(YeeGrid& grid) {

    calcD(grid.D_Hx, grid.delta_t, grid.mu_Hx, grid.sigma_Hx);
    calcD(grid.D_Hy, grid.delta_t, grid.mu_Hy, grid.sigma_Hy);
    calcD(grid.D_Hz, grid.delta_t, grid.mu_Hz, grid.sigma_Hz);

    calcC(grid.C_Ex, grid.delta_t, grid.epsilon_Ex, grid.sigma_Ex);
    calcC(grid.C_Ey, grid.delta_t, grid.epsilon_Ey, grid.sigma_Ey);
    calcC(grid.C_Ez, grid.delta_t, grid.epsilon_Ez, grid.sigma_Ez);

    calcD(grid.D_Ex, grid.delta_t, grid.epsilon_Ex, grid.sigma_Ex);
    calcD(grid.D_Ey, grid.delta_t, grid.epsilon_Ey, grid.sigma_Ey);
    calcD(grid.D_Ez, grid.delta_t, grid.epsilon_Ez, grid.sigma_Ez);

    calcDelta(grid.delta_y_Hx, grid.y_Hx);
    calcDelta(grid.delta_z_Hx, grid.z_Hx);
    calcDelta(grid.delta_x_Hy, grid.x_Hy);
    calcDelta(grid.delta_z_Hy, grid.z_Hy);
    calcDelta(grid.delta_x_Hz, grid.x_Hz);
    calcDelta(grid.delta_y_Hz, grid.y_Hz);

    calcDelta(grid.delta_y_Ex, grid.y_Ex);
    calcDelta(grid.delta_z_Ex, grid.z_Ex);
    calcDelta(grid.delta_x_Ey, grid.x_Ey);
    calcDelta(grid.delta_z_Ey, grid.z_Ey);
    calcDelta(grid.delta_x_Ez, grid.x_Ez);
    calcDelta(grid.delta_y_Ez, grid.y_Ez);
}

void calcH(YeeGrid& grid) {

    int nx = grid.Hx.getCountX();
    int ny = grid.Hx.getCountY();
    int nz = grid.Hx.getCountZ();

    for(int ix = 0; ix < nx-1; ix++)
    for(int iy = 0; iy < ny-1; iy++)
    for(int iz = 0; iz < nz-1; iz++) {
        float& curHx      = grid.Hx          .at(ix,   iy,   iz);
        float  curD_Hx    = grid.D_Hx        .at(ix,   iy,   iz);
        float  curEz1     = grid.Ez          .at(ix,   iy+1, iz);
        float  curEz0     = grid.Ez          .at(ix,   iy,   iz);
        float  curEy1     = grid.Ey          .at(ix,   iy,   iz+1);
        float  curEy0     = grid.Ey          .at(ix,   iy,   iz);
        float  curDeltaEz = grid.delta_y_Ez  .at(iy);
        float  curDeltaEy = grid.delta_z_Ey  .at(iz);

        curHx -= curD_Hx * ((curEz1 - curEz0) / curDeltaEz -
                            (curEy1 - curEy0) / curDeltaEy);
    }

    for(int ix = 0; ix < nx-1; ix++)
    for(int iy = 0; iy < ny-1; iy++)
    for(int iz = 0; iz < nz-1; iz++) {
        float& curHy      = grid.Hy          .at(ix,   iy,   iz);
        float  curD_Hy    = grid.D_Hy        .at(ix,   iy,   iz);
        float  curEx1     = grid.Ex          .at(ix,   iy,   iz+1);
        float  curEx0     = grid.Ex          .at(ix,   iy,   iz);
        float  curEz1     = grid.Ez          .at(ix+1, iy,   iz);
        float  curEz0     = grid.Ez          .at(ix,   iy,   iz);
        float  curDeltaEx = grid.delta_z_Ex  .at(iz);
        float  curDeltaEz = grid.delta_x_Ez  .at(ix);

        curHy -= curD_Hy * ((curEx1 - curEx0) / curDeltaEx -
                            (curEz1 - curEz0) / curDeltaEz);
    }

    for(int ix = 0; ix < nx-1; ix++)
    for(int iy = 0; iy < ny-1; iy++)
    for(int iz = 0; iz < nz-1; iz++) {
        float& curHz      = grid.Hz          .at(ix,   iy,   iz);
        float  curD_Hz    = grid.D_Hz        .at(ix,   iy,   iz);
        float  curEy1     = grid.Ey          .at(ix+1, iy,   iz);
        float  curEy0     = grid.Ey          .at(ix,   iy,   iz);
        float  curEx1     = grid.Ex          .at(ix,   iy+1, iz);
        float  curEx0     = grid.Ex          .at(ix,   iy,   iz);
        float  curDeltaEy = grid.delta_x_Ey  .at(ix);
        float  curDeltaEx = grid.delta_y_Ex  .at(iy);

        curHz -= curD_Hz * ((curEy1 - curEy0) / curDeltaEy -
                            (curEx1 - curEx0) / curDeltaEx);
    }
}

void calcE(YeeGrid& grid) {

    int nx = grid.Ex.getCountX();
    int ny = grid.Ex.getCountY();
    int nz = grid.Ex.getCountZ();

    for(int ix = 1; ix < nx; ix++)
    for(int iy = 1; iy < ny; iy++)
    for(int iz = 1; iz < nz; iz++) {
        float& curEx      = grid.Ex          .at(ix,   iy,   iz);
        float  curC_Ex    = grid.C_Ex        .at(ix,   iy,   iz);
        float  curD_Ex    = grid.D_Ex        .at(ix,   iy,   iz);
        float  curHz0     = grid.Hz          .at(ix,   iy,   iz);
        float  curHz1     = grid.Hz          .at(ix,   iy-1, iz);
        float  curHy0     = grid.Hy          .at(ix,   iy,   iz);
        float  curHy1     = grid.Hy          .at(ix,   iy,   iz-1);
        float  curDeltaHz = grid.delta_y_Hz  .at(iy-1);
        float  curDeltaHy = grid.delta_z_Hy  .at(iz-1);

        curEx = curC_Ex * curEx + curD_Ex * ((curHz0 - curHz1) / curDeltaHz -
                                             (curHy0 - curHy1) / curDeltaHy);
    }

    for(int ix = 1; ix < nx; ix++)
    for(int iy = 1; iy < ny; iy++)
    for(int iz = 1; iz < nz; iz++) {
        float& curEy      = grid.Ey          .at(ix,   iy,   iz);
        float  curC_Ey    = grid.C_Ey        .at(ix,   iy,   iz);
        float  curD_Ey    = grid.D_Ey        .at(ix,   iy,   iz);
        float  curHx0     = grid.Hx          .at(ix,   iy,   iz);
        float  curHx1     = grid.Hx          .at(ix,   iy,   iz-1);
        float  curHz0     = grid.Hz          .at(ix,   iy,   iz);
        float  curHz1     = grid.Hz          .at(ix-1, iy,   iz);
        float  curDeltaHx = grid.delta_z_Hx  .at(iz-1);
        float  curDeltaHz = grid.delta_x_Hz  .at(ix-1);

        curEy = curC_Ey * curEy + curD_Ey * ((curHx0 - curHx1) / curDeltaHx -
                                             (curHz0 - curHz1) / curDeltaHz);
    }

    for(int ix = 1; ix < nx; ix++)
    for(int iy = 1; iy < ny; iy++)
    for(int iz = 1; iz < nz; iz++) {
        float& curEz      = grid.Ez          .at(ix,   iy,   iz);
        float  curC_Ez    = grid.C_Ez        .at(ix,   iy,   iz);
        float  curD_Ez    = grid.D_Ez        .at(ix,   iy,   iz);
        float  curHy0     = grid.Hy          .at(ix,   iy,   iz);
        float  curHy1     = grid.Hy          .at(ix-1, iy,   iz);
        float  curHx0     = grid.Hx          .at(ix,   iy,   iz);
        float  curHx1     = grid.Hx          .at(ix,   iy-1, iz);
        float  curDeltaHy = grid.delta_x_Hy  .at(ix-1);
        float  curDeltaHx = grid.delta_y_Hx  .at(iy-1);

        curEz = curC_Ez * curEz + curD_Ez * ((curHy0 - curHy1) / curDeltaHy -
                                             (curHx0 - curHx1) / curDeltaHx);
    }
}
