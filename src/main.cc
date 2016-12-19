#include "read_grid.hh"
#include "calc_field.hh"
#include "resistive_source.hh"
#include <iostream>
#include <fstream>
#include <boost/format.hpp>
#include "EasyCL.h"
#include <chrono>
#include <cmath>

const float c = 299792458.0f;
const float pi = 3.14159265358;

int nx;
int ny;
int nz;
const int nl = 20;
const float lambda = 0.05;
const float omega  = 2 * pi * c / lambda;
const float dx = lambda / nl;
const float dy = dx;
const float dz = dx;
const float dt = 0.1/(c * std::sqrt(1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz)));

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

    //std::cout << "MAX: " << maxEx << std::endl;

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

void SetCube(YeeGrid& grid, float eps, float sigmaE, int Imin, int Imax, int Jmin, int Jmax, int Kmin, int Kmax)
{
    long int i, j, k;

    for(i=Imin;i<Imax;i++)
    {
        for(j=Jmin;j<Jmax;j++)
        {
            for(k=Kmin;k<Kmax;k++)
            {
                grid.epsilon_Ex.at(i, k, k) = eps;
                grid.epsilon_Ey.at(i, k, k) = eps;
                grid.epsilon_Ez.at(i, k, k) = eps;
                grid.sigma_Ex.at(i, k, k) = sigmaE;
                grid.sigma_Ey.at(i, k, k) = sigmaE;
                grid.sigma_Ey.at(i, k, k) = sigmaE;
            }
        }
    }
}

void setupGrid(YeeGrid& grid) {
    int x0 = nx / 2;
    int y0 = ny / 2;
    int z0 = nz / 2;
    const int w = 3;
    SetCube(grid, 10000, 10, x0-5, x0+5, y0-5, y0+5, z0-20, z0 + 20);
}

float voltage(float time) {
    return std::sin(omega*time) * 1000;
}

template<typename TimeT = std::chrono::milliseconds>
struct measure
{
    template<typename F, typename ...Args>
    static typename TimeT::rep execution(F func, Args&&... args)
    {
        auto start = std::chrono::system_clock::now();

        // Now call the function with all the parameters you need.
        func(std::forward<Args>(args)...);

        auto duration = std::chrono::duration_cast< TimeT>
                            (std::chrono::system_clock::now() - start);

        return duration.count();
    }
};


int cpu_main() {
    YeeGrid grid(nx, ny, nz, dt, dx, dy, dz);
    setupGrid(grid);
    calcCoefs(grid);

    //dumpImage("test.ppm", grid);

    int x0 = nx / 2;
    int y0 = ny / 2;
    int z0 = nz / 2;
    ResistiveSource rsource(x0, y0, z0, 10.0f);
    rsource.calcCoefs(grid);

    float time = 0;
    for (int iter = 0; iter < 300; ++iter) {

        std::cout << "Rocking iteration #" << iter << std::endl;
        std::cout << "calcH_time: "
                  << measure<std::chrono::milliseconds>::execution([&grid]() {calcH(grid);}) << std::endl;

        rsource.resqueFields(grid);
        std::cout << "calcE_time: "
                  << measure<std::chrono::milliseconds>::execution([&grid]() {calcE(grid);}) << std::endl;

        rsource.updateFields(grid, voltage(time));

        std::string filename = str(boost::format("field.%03d_%03d.ppm") % nx % iter);
        std::cout << "VAL: " << grid.Ez.at(x0 + 10, y0, z0) << std::endl;

        dumpImage(filename, grid);
        time = dt*iter;
    }

    return 0;
}

int gpu_main() {
    YeeGrid grid(nx, ny, nz, dt, dx, dy, dz);
    setupGrid(grid);
    calcCoefs(grid);

    int x0 = nx / 2;
    int y0 = ny / 2;
    int z0 = nz / 2;
    ResistiveSource rsource(x0, y0, z0, 10.0f);
    rsource.calcCoefs(grid);

    if( !EasyCL::isOpenCLAvailable() ) {
        std::cout << "opencl library not found" << std::endl;
        return -1;
    }

    EasyCL *cl = EasyCL::createForFirstGpu();
    cl->setProfiling(true);

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


    CLKernel *kernelE = cl->buildKernel("../src/kernels.cl", "calcE");
    CLKernel *kernelH = cl->buildKernel("../src/kernels.cl", "calcH");

    float prevE = 0;
    CLWrapper *wprevE = cl->wrap(1, &prevE);
    wprevE->copyToDevice();

    CLKernel *rescueField = cl->buildKernel("../src/kernels.cl", "rescueField");
    CLKernel *updateField = cl->buildKernel("../src/kernels.cl", "updateField");



    std::vector<size_t> global_dims = { nx - 1, ny - 1, nz - 1 },
                        local_dims  = { 8, 8, 8 };



    float time = 0;
    for (int iter = 0; iter < 300; ++iter) {
        time = iter * dt;


    std::cout << "step: " << iter << std::endl;
    std::cout << "  voltage: " << voltage(time) << std::endl;
    kernelH
            ->in((int)grid.Hx.getCountX())
            ->in((int)grid.Hx.getCountY())
            ->in((int)grid.Hx.getCountZ())
            ->in(grid.delta_x)
            ->in(grid.delta_y)
            ->in(grid.delta_z)
            ->in(wEx)
            ->in(wEy)
            ->in(wEz)
            ->inout(wHx)
            ->inout(wHy)
            ->inout(wHz)
            ->in(wD_Hx)
            ->in(wD_Hy)
            ->in(wD_Hz)
            ->run(3, global_dims.data(), local_dims.data());

    cl->finish();
    rescueField->in(nx / 2)
               ->in(ny/2)
               ->in(nz/2)
               ->in(nx)
               ->in(ny)
               ->in(nz)
               ->in(wEz)
               ->out(wprevE)
               ->run_1d(1, 1);
    //rsource.resqueFields(grid);

    kernelE->in((int)grid.Ex.getCountX())
            ->in((int)grid.Ex.getCountY())
            ->in((int)grid.Ex.getCountZ())
            ->in(grid.delta_x)
            ->in(grid.delta_y)
            ->in(grid.delta_z)
            ->out(wEx)
            ->out(wEy)
            ->out(wEz)
            ->in(wHx)
            ->in(wHy)
            ->in(wHz)
            ->in(wC_Ex)
            ->in(wC_Ey)
            ->in(wC_Ez)
            ->in(wD_Ex)
            ->in(wD_Ey)
            ->in(wD_Ez)
            ->run(3, global_dims.data(), local_dims.data());

    cl->finish();
    updateField->in(nx / 2)
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
                   ->in(grid.delta_y)
                   ->run_1d(1, 1);

    cl->finish();

    cl->dumpProfiling();
    wEz->copyToHost();

    //rsource.updateFields(grid, voltage(time));
    //wEz->copyToDevice();

    std::string filename = str(boost::format("CLfield.%03d_%03d.ppm") % nx % iter);
    std::cout << time << '\t' << grid.Ez.at(x0 + 10, y0, z0) << std::endl;
    dumpImage(filename, grid);

    }

    delete kernelH;
    delete kernelE;

    return 0;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cout << "Invalid command line parameters number.\n";
        return -1;
    }
    std::string work_mode = argv[1];
    nx = atoi(argv[2]);
    ny = atoi(argv[2]);
    nz = atoi(argv[2]);


    if (work_mode == "cpu")
        return cpu_main();
    if (work_mode == "gpu")
        return gpu_main();
}
