#include "yee_grid.hh"

using Array = std::vector<float>;
using SolidArray3d = rvlm::core::SolidArray3d<float>;

void calcD(SolidArray3d& D,float deltaT,
           const SolidArray3d& perm, const SolidArray3d& sigma) {
    for(int ix = 0; ix < D.getCountX(); ix++)
        for(int iy = 0; iy < D.getCountY(); iy++)
            for (int iz = 0; iz < D.getCountZ(); iz++) {
                D.at(ix, iy, iz) = (deltaT / perm.at(ix, iy, iz)) /
                        (1 + (sigma.at(ix, iy, iz) * deltaT) /
                         (2 * perm.at(ix, iy, iz)));
            }

}

void calcC(SolidArray3d& E,float deltaT,
           const SolidArray3d& perm, const SolidArray3d& sigma) {
    for(int ix = 0; ix < E.getCountX(); ix++)
        for(int iy = 0; iy < E.getCountY(); iy++)
            for (int iz = 0; iz < E.getCountZ(); iz++) {
                E.at(ix, iy, iz) = (1 - (sigma.at(ix, iy, iz) * deltaT) /
                                    (2 * perm.at(ix, iy, iz))) /
                        (1 + (sigma.at(ix, iy, iz) * deltaT) /
                         (2 * perm.at(ix, iy, iz)));
            }

}

void calcDelta(Array& delta, Array& input) {
    for (int i = 0; i < delta.size() - 1; i++)
            delta[i] = input[i + 1] - input[i];
    delta[delta.size() - 1] = delta[delta.size() - 2];
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
    for(int i = 0; i < grid.Hx.getCountX(); i++)
    for(int j = 0; j < grid.Hx.getCountY(); j++)
    for(int k = 0; k < grid.Hx.getCountZ(); k++) {
        grid.Hx.at(i, j, k) -= grid.D_Hx.at(i, j, k) *
                ((grid.Ez.at(i, j + 1, k) - grid.Ez.at(i, j, k)) / grid.delta_y_Ez[j] -
                 (grid.Ey.at(i, j, k + 1) - grid.Ey.at(i, j, k)) / grid.delta_z_Ey[k]);
    }
    for(int i = 0; i < grid.Hy.getCountX(); i++)
    for(int j = 0; j < grid.Hy.getCountY(); j++)
    for(int k = 0; k < grid.Hy.getCountZ(); k++) {
        grid.Hy.at(i, j, k) -= grid.D_Hy.at(i, j, k) *
                ((grid.Ex.at(i, j, k + 1) - grid.Ex.at(i, j, k)) / grid.delta_z_Ex[i] -
                 (grid.Ez.at(i + 1, j, k) - grid.Ez.at(i, j, k)) / grid.delta_x_Ez[k]);
    }
    for(int i = 0; i < grid.Hz.getCountX(); i++)
    for(int j = 0; j < grid.Hz.getCountY(); j++)
    for(int k = 0; k < grid.Hz.getCountZ(); k++) {
        grid.Hz.at(i, j, k) -= grid.D_Hz.at(i, j, k) *
                ((grid.Ey.at(i + 1, j, k) - grid.Ey.at(i, j, k)) / grid.delta_x_Ey[j] -
                 (grid.Ex.at(i, j + 1, k) - grid.Ex.at(i, j, k)) / grid.delta_y_Ex[i]);
    }
}

void calcE(YeeGrid& grid) {
    for(int i = 0; i < grid.Ex.getCountX(); i++)
    for(int j = 0; j < grid.Ex.getCountY(); j++)
    for(int k = 0; k < grid.Ex.getCountZ(); k++) {
        grid.Ex.at(i, j, k) *= grid.C_Ex.at(i, j, k);
        grid.Ex.at(i, j, k) += grid.D_Ex.at(i, j, k) *
                ((grid.Hz.at(i, j, k) - grid.Hz.at(i, j - 1, k)) / grid.delta_y_Hz[j] -
                 (grid.Hy.at(i, j, k) - grid.Hy.at(i, j, k - 1)) / grid.delta_z_Hy[k]);
    }
    for(int i = 0; i < grid.Ey.getCountX(); i++)
    for(int j = 0; j < grid.Ey.getCountY(); j++)
    for(int k = 0; k < grid.Ey.getCountZ(); k++) {
        grid.Ey.at(i, j, k) *= grid.C_Ey.at(i, j, k);
        grid.Ey.at(i, j, k) += grid.D_Ey.at(i, j, k) *
                ((grid.Hx.at(i, j, k) - grid.Hx.at(i, j, k - 1)) / grid.delta_z_Hx[i] -
                 (grid.Hz.at(i, j, k) - grid.Hz.at(i - 1, j, k)) / grid.delta_x_Hz[j]);
    }
    for(int i = 0; i < grid.Ez.getCountX(); i++)
    for(int j = 0; j < grid.Ez.getCountY(); j++)
    for(int k = 0; k < grid.Ez.getCountZ(); k++) {
        grid.Ez.at(i, j, k) *= grid.C_Ez.at(i, j, k);
        grid.Ez.at(i, j, k) += grid.D_Ez.at(i, j, k) *
                ((grid.Hy.at(i, j, k) - grid.Hy.at(i - 1, j, k)) / grid.delta_x_Hy[j] -
                 (grid.Hx.at(i, j, k) - grid.Hx.at(i, j - 1, k)) / grid.delta_y_Hx[i]);
    }

}
