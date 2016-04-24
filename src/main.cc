#include "read_grid.hh"

int main(int argc, char *argv[]) {
    YeeGrid grid = readGridData("test.nc", 60e-12);

    return 0;
}
