#include "rvlm/core/SolidArray3d.hh"
#include "yee_grid.hh"

template <typename CurlCoefficientT, typename TimeT, typename PermittivityT, typename ConductivityT>
void calcD(rvlm::core::SolidArray3d<CurlCoefficientT>& D,
           TimeT deltaT,
           rvlm::core::SolidArray3d<PermittivityT> const& perm,
           rvlm::core::SolidArray3d<ConductivityT> const& sigma) {

    for (int ix = 0; ix < D.getCountX(); ix++)
    for (int iy = 0; iy < D.getCountY(); iy++)
    for (int iz = 0; iz < D.getCountZ(); iz++) {
        auto& curD     = D.at(ix, iy, iz);
        auto  curPerm  = perm.at(ix, iy, iz);
        auto  curSigma = sigma.at(ix, iy, iz);
        auto  subexpr  = deltaT / curPerm;
        curD = subexpr / (1 + curSigma * subexpr/2);
    }
}

template <typename CurlCoefficientT, typename TimeT, typename PermittivityT, typename ConductivityT>
void calcC(rvlm::core::SolidArray3d<CurlCoefficientT>& C,
           TimeT deltaT,
           rvlm::core::SolidArray3d<PermittivityT> const& perm,
           rvlm::core::SolidArray3d<ConductivityT> const& sigma) {

    for (int ix = 0; ix < C.getCountX(); ix++)
    for (int iy = 0; iy < C.getCountY(); iy++)
    for (int iz = 0; iz < C.getCountZ(); iz++) {
        auto& curC     = C.at(ix, iy, iz);
        auto  curPerm  = perm.at(ix, iy, iz);
        auto  curSigma = sigma.at(ix, iy, iz);
        auto  subexpr  = curSigma * deltaT / (2 * curPerm);
        curC = (1 - subexpr) / (1 + subexpr);
    }
}

template <typename valueType>
void calcCoefs(YeeGrid<valueType>& grid) {

    calcD(grid.D_Hx, grid.delta_t, grid.mu_Hx, grid.sigma_Hx);
    calcD(grid.D_Hy, grid.delta_t, grid.mu_Hy, grid.sigma_Hy);
    calcD(grid.D_Hz, grid.delta_t, grid.mu_Hz, grid.sigma_Hz);

    calcC(grid.C_Ex, grid.delta_t, grid.epsilon_Ex, grid.sigma_Ex);
    calcC(grid.C_Ey, grid.delta_t, grid.epsilon_Ey, grid.sigma_Ey);
    calcC(grid.C_Ez, grid.delta_t, grid.epsilon_Ez, grid.sigma_Ez);

    calcD(grid.D_Ex, grid.delta_t, grid.epsilon_Ex, grid.sigma_Ex);
    calcD(grid.D_Ey, grid.delta_t, grid.epsilon_Ey, grid.sigma_Ey);
    calcD(grid.D_Ez, grid.delta_t, grid.epsilon_Ez, grid.sigma_Ez);
}

template <typename T>
void calcH(YeeGrid<T>& grid) {

    int nx = grid.Hx.getCountX();
    int ny = grid.Hx.getCountY();
    int nz = grid.Hx.getCountZ();

    auto delta_x = grid.delta_x;
    auto delta_y = grid.delta_y;
    auto delta_z = grid.delta_z;

    for(int ix = 0; ix < nx-1; ix++)
    for(int iy = 0; iy < ny-1; iy++)
    for(int iz = 0; iz < nz-1; iz++) {
        auto& curHx      = grid.Hx          .at(ix,   iy,   iz);
        auto  curD_Hx    = grid.D_Hx        .at(ix,   iy,   iz);
        auto  curEz1     = grid.Ez          .at(ix,   iy+1, iz);
        auto  curEz0     = grid.Ez          .at(ix,   iy,   iz);
        auto  curEy1     = grid.Ey          .at(ix,   iy,   iz+1);
        auto  curEy0     = grid.Ey          .at(ix,   iy,   iz);

        curHx -= curD_Hx * ((curEz1 - curEz0) / delta_y -
                            (curEy1 - curEy0) / delta_z);
    }

    for(int ix = 0; ix < nx-1; ix++)
    for(int iy = 0; iy < ny-1; iy++)
    for(int iz = 0; iz < nz-1; iz++) {
        auto& curHy      = grid.Hy          .at(ix,   iy,   iz);
        auto  curD_Hy    = grid.D_Hy        .at(ix,   iy,   iz);
        auto  curEx1     = grid.Ex          .at(ix,   iy,   iz+1);
        auto  curEx0     = grid.Ex          .at(ix,   iy,   iz);
        auto  curEz1     = grid.Ez          .at(ix+1, iy,   iz);
        auto  curEz0     = grid.Ez          .at(ix,   iy,   iz);

        curHy -= curD_Hy * ((curEx1 - curEx0) / delta_z -
                            (curEz1 - curEz0) / delta_x);
    }

    for(int ix = 0; ix < nx-1; ix++)
    for(int iy = 0; iy < ny-1; iy++)
    for(int iz = 0; iz < nz-1; iz++) {
        auto& curHz      = grid.Hz          .at(ix,   iy,   iz);
        auto  curD_Hz    = grid.D_Hz        .at(ix,   iy,   iz);
        auto  curEy1     = grid.Ey          .at(ix+1, iy,   iz);
        auto  curEy0     = grid.Ey          .at(ix,   iy,   iz);
        auto  curEx1     = grid.Ex          .at(ix,   iy+1, iz);
        auto  curEx0     = grid.Ex          .at(ix,   iy,   iz);

        curHz -= curD_Hz * ((curEy1 - curEy0) / delta_x -
                            (curEx1 - curEx0) / delta_y);
    }
}

template <typename T>
void calcE(YeeGrid<T>& grid) {

    int nx = grid.Ex.getCountX();
    int ny = grid.Ex.getCountY();
    int nz = grid.Ex.getCountZ();

    auto delta_x = grid.delta_x;
    auto delta_y = grid.delta_y;
    auto delta_z = grid.delta_z;

    for(int ix = 1; ix < nx; ix++)
    for(int iy = 1; iy < ny; iy++)
    for(int iz = 1; iz < nz; iz++) {
        auto& curEx      = grid.Ex          .at(ix,   iy,   iz);
        auto  curC_Ex    = grid.C_Ex        .at(ix,   iy,   iz);
        auto  curD_Ex    = grid.D_Ex        .at(ix,   iy,   iz);
        auto  curHz0     = grid.Hz          .at(ix,   iy,   iz);
        auto  curHz1     = grid.Hz          .at(ix,   iy-1, iz);
        auto  curHy0     = grid.Hy          .at(ix,   iy,   iz);
        auto  curHy1     = grid.Hy          .at(ix,   iy,   iz-1);

        curEx = curC_Ex * curEx + curD_Ex * ((curHz0 - curHz1) / delta_y -
                                             (curHy0 - curHy1) / delta_z);
    }

    for(int ix = 1; ix < nx; ix++)
    for(int iy = 1; iy < ny; iy++)
    for(int iz = 1; iz < nz; iz++) {
        auto& curEy      = grid.Ey          .at(ix,   iy,   iz);
        auto  curC_Ey    = grid.C_Ey        .at(ix,   iy,   iz);
        auto  curD_Ey    = grid.D_Ey        .at(ix,   iy,   iz);
        auto  curHx0     = grid.Hx          .at(ix,   iy,   iz);
        auto  curHx1     = grid.Hx          .at(ix,   iy,   iz-1);
        auto  curHz0     = grid.Hz          .at(ix,   iy,   iz);
        auto  curHz1     = grid.Hz          .at(ix-1, iy,   iz);

        curEy = curC_Ey * curEy + curD_Ey * ((curHx0 - curHx1) / delta_z -
                                             (curHz0 - curHz1) / delta_x);
    }

    for(int ix = 1; ix < nx; ix++)
    for(int iy = 1; iy < ny; iy++)
    for(int iz = 1; iz < nz; iz++) {
        auto& curEz      = grid.Ez          .at(ix,   iy,   iz);
        auto  curC_Ez    = grid.C_Ez        .at(ix,   iy,   iz);
        auto  curD_Ez    = grid.D_Ez        .at(ix,   iy,   iz);
        auto  curHy0     = grid.Hy          .at(ix,   iy,   iz);
        auto  curHy1     = grid.Hy          .at(ix-1, iy,   iz);
        auto  curHx0     = grid.Hx          .at(ix,   iy,   iz);
        auto  curHx1     = grid.Hx          .at(ix,   iy-1, iz);

        curEz = curC_Ez * curEz + curD_Ez * ((curHy0 - curHy1) / delta_x -
                                             (curHx0 - curHx1) / delta_y);
    }
}
