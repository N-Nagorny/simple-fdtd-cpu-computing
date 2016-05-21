#include "read_grid.hh"
#include "calc_field.hh"
#include "resistive_source.hh"
#include <iostream>
#include <fstream>
#include <boost/format.hpp>
#include "EasyCL.h"

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

void dumpImage(std::string const& filename, YeeGrid& grid) {
    int nx = grid.Ez.getCountX();
    int ny = grid.Ez.getCountY();
    int z0 = grid.Ex.getCountZ() / 2;

    float maxEx = 0;
    for (int ix = 0; ix < nx; ++ix)
    for (int iy = 0; iy < ny; ++iy) {
        if ( ((nx/2 - 3) <= ix && ix <= (nx/2 + 3)) &&
             ((ny/2 - 3) <= iy && iy <= (ny/2 + 3)) )
            continue;

        float curAbs = std::abs(grid.Ez.at(ix, iy, z0));
        if (curAbs > maxEx)
            maxEx = curAbs;
    }

    std::cout << "MAX: " << maxEx << std::endl;

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
        color[1] = 0;           // green
        color[2] = blue;  // blue
        output.write(&color[0], 3);
    }
}

void setupGrid(YeeGrid& grid) {
    int x0 = nx / 2;
    int y0 = ny / 2;
    int z0 = nz / 2;

    float epsilonAntenna = 10;
    float sigmaAntenna = 10000;
    for (int i = 0; i < nl; ++i) {
        grid.epsilon_Ez.at(x0, y0, z0 + i) = epsilonAntenna;
        grid.epsilon_Ez.at(x0, y0, z0 - i) = epsilonAntenna;
        grid.sigma_Ez.at(x0, y0, z0 + i) = sigmaAntenna;
        grid.sigma_Ez.at(x0, y0, z0 - i) = sigmaAntenna;
    }
}

float voltage(float time) {
    return std::sin(omega*time) * 1000;
}

int main_2(int argc, char *argv[]) {
    YeeGrid grid(nx, ny, nz, dt, dx, dy, dz);
    setupGrid(grid);
    calcCoefs(grid);

    dumpImage("test.ppm", grid);

    int x0 = nx / 2;
    int y0 = ny / 2;
    int z0 = nz / 2;
    ResistiveSource rsource(x0, y0, z0, 10);
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

int main(int argc, char *argv[]) {
    YeeGrid grid(nx, ny, nz, dt, dx, dy, dz);
    setupGrid(grid);
    calcCoefs(grid);
    ResistiveSource rsource(nx/2, ny/2, nz/2, 10);
    rsource.calcCoefs(grid);

    if( !EasyCL::isOpenCLAvailable() ) {
        std::cout << "opencl library not found" << std::endl;
        return -1;
    }

    EasyCL *cl = EasyCL::createForFirstGpu();

    CLWrapper *wEx = cl->wrap(grid.Ex.getTotalCount(), &grid.Ex.at(0,0,0));
    wEx->copyToDevice();
    CLWrapper *wEy = cl->wrap(grid.Ey.getTotalCount(), &grid.Ey.at(0,0,0));
    wEy->copyToDevice();
    CLWrapper *wEz = cl->wrap(grid.Ez.getTotalCount(), &grid.Ez.at(0,0,0));
    wEz->copyToDevice();

    CLWrapper *wHx = cl->wrap(grid.Hx.getTotalCount(), &grid.Hx.at(0,0,0));
    wHx->copyToDevice();
    CLWrapper *wHy = cl->wrap(grid.Hy.getTotalCount(), &grid.Hy.at(0,0,0));
    wHy->copyToDevice();
    CLWrapper *wHz = cl->wrap(grid.Hz.getTotalCount(), &grid.Hz.at(0,0,0));
    wHz->copyToDevice();

    CLWrapper *wC_Ex = cl->wrap(grid.C_Ex.getTotalCount(), &grid.C_Ex.at(0,0,0));
    wC_Ex->copyToDevice();
    CLWrapper *wC_Ey = cl->wrap(grid.C_Ey.getTotalCount(), &grid.C_Ey.at(0,0,0));
    wC_Ey->copyToDevice();
    CLWrapper *wC_Ez = cl->wrap(grid.C_Ez.getTotalCount(), &grid.C_Ez.at(0,0,0));
    wC_Ez->copyToDevice();

    CLWrapper *wD_Ex = cl->wrap(grid.D_Ex.getTotalCount(), &grid.D_Ex.at(0,0,0));
    wD_Ex->copyToDevice();
    CLWrapper *wD_Ey = cl->wrap(grid.D_Ey.getTotalCount(), &grid.D_Ey.at(0,0,0));
    wD_Ey->copyToDevice();
    CLWrapper *wD_Ez = cl->wrap(grid.D_Ez.getTotalCount(), &grid.D_Ez.at(0,0,0));
    wD_Ez->copyToDevice();


    CLWrapper *wD_Hx = cl->wrap(grid.D_Hx.getTotalCount(), &grid.D_Hx.at(0,0,0));
    wD_Hx->copyToDevice();
    CLWrapper *wD_Hy = cl->wrap(grid.D_Hy.getTotalCount(), &grid.D_Hy.at(0,0,0));
    wD_Hy->copyToDevice();
    CLWrapper *wD_Hz = cl->wrap(grid.D_Hz.getTotalCount(), &grid.D_Hz.at(0,0,0));
    wD_Hz->copyToDevice();


    std::cout << "fuck kernel" << std::endl;
    CLKernel *kernelE = cl->buildKernel("../src/kernels.cl", "calcE");

    kernelE->in((int)grid.Ex.getCountX());
    std::cout << "fuck kernel" << std::endl;
    kernelE->in((int)grid.Ex.getCountY());
    kernelE->in((int)grid.Ex.getCountZ());
    kernelE->in(grid.delta_x);
    kernelE->in(grid.delta_y);
    kernelE->in(grid.delta_z);
    std::cout << "fuck kernel" << std::endl;
    kernelE->out(wEx);
    kernelE->out(wEy);
    kernelE->out(wEz);
    kernelE->in(wHx);
    kernelE->in(wHy);
    kernelE->in(wHz);
    kernelE->in(wC_Ex);
    kernelE->in(wC_Ey);
    kernelE->in(wC_Ez);
    kernelE->in(wD_Ex);
    kernelE->in(wD_Ey);
    kernelE->in(wD_Ez);

    CLKernel *kernelH = cl->buildKernel("../src/kernels.cl", "calcH");

    kernelH->in((int)grid.Hx.getCountX());
    kernelH->in((int)grid.Hx.getCountY());
    kernelH->in((int)grid.Hx.getCountZ());
    kernelH->in(grid.delta_x);
    kernelH->in(grid.delta_y);
    kernelH->in(grid.delta_z);
    kernelH->in(wEx);
    kernelH->in(wEy);
    kernelH->in(wEz);
    kernelH->inout(wHx);
    kernelH->inout(wHy);
    kernelH->inout(wHz);
    kernelH->in(wD_Hx);
    kernelH->in(wD_Hy);
    kernelH->in(wD_Hz);

    float prevE = 0;
    CLWrapper *wprevE = cl->wrap(1, &prevE);
    wprevE->copyToDevice();

    CLKernel *rescueField = cl->buildKernel("../src/kernels.cl", "rescueField");
    rescueField->in(nx / 2)
               ->in(ny/2)
               ->in(nz/2)
               ->in(nx)
               ->in(ny)
               ->in(nz)
               ->in(wEz)
               ->out(wprevE);
    CLKernel *updateField = cl->buildKernel("../src/kernels.cl", "updateField");
    cl->storeKernel("updateField", updateField);



    std::vector<size_t> global_dims = { 128, 128, 128 },
                        local_dims  = { 8, 8, 8 };



    float time = 0;
    for (int iter = 0; iter < 300; ++iter) {
        time = iter * dt;

    CLKernel *upd = cl->getKernel("updateField");
        upd->in(nx / 2)
                   ->in(ny / 2)
                   ->in(nz / 2)
                   ->in(nx)
                   ->in(ny)
                   ->in(nz)
                   ->out(wEz)
                   ->in(wHx)
                   ->in(wHy)
                   ->in(wprevE)
                   ->in(rsource.getC())
                   ->in(rsource.getD())
                   ->in(10.f) // Ohm
                   ->in(voltage(time))
                   ->in(grid.delta_x)
                   ->in(grid.delta_y);


    std::cout << "step: " << iter << std::endl;
    std::cout << "  voltage: " << voltage(time) << std::endl;
    kernelH->run(3, global_dims.data(), local_dims.data());

    rescueField->run_1d(1,1);

    kernelE->run(3, global_dims.data(), local_dims.data());

    upd->run_1d(1, 1);
    cl->finish();

    wEz->copyToHost();
    std::string filename = str(boost::format("clfield_%04d.ppm") % iter);
    dumpImage(filename, grid);

    }

    delete kernelH;
    delete kernelE;

    return 0;
}
