#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <errno.h>
#include "yee_grid.h"
#include "read_grid.h"

int main(int argc, char *argv[]) {
    struct yee_grid *grid = NULL;
    int i;
    int retcode = 1;

    grid = read_fdtd_field("input/filename_%s.txt");
    if (!grid)
        goto e_out;

    for (i = 0; i < grid->nx_Ex; i++) {
        printf("%f\n", grid->x_Ex[i]);
    }

    printf("\n%zu\n", grid->nx_Ex);

    retcode = 0;

e_out:
    free_fdtd_struct(grid);
    free(grid);
    return retcode;
}
