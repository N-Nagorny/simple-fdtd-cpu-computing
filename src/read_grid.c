#include <stdlib.h>
#include <stdio.h>
#include "yee_grid.h"
#include "read_grid.h"

#define FDTD_FILENAME_BUF_SIZE 1024

static void init_fdtd_struct(struct yee_grid *grid) {
    grid->x_Ex = NULL;
    grid->x_Ey = NULL;
    grid->x_Ez = NULL;
    grid->y_Ex = NULL;
    grid->y_Ey = NULL;
    grid->y_Ez = NULL;
    grid->z_Ex = NULL;
    grid->z_Ey = NULL;
    grid->z_Ez = NULL;
    grid->epsilon_Ex = NULL;
    grid->epsilon_Ey = NULL;
    grid->epsilon_Ez = NULL;
    grid->sigmaE_Ex  = NULL;
    grid->sigmaE_Ey  = NULL;
    grid->sigmaE_Ez  = NULL;

    grid->x_Hx = NULL;
    grid->x_Hy = NULL;
    grid->x_Hz = NULL;
    grid->y_Hx = NULL;
    grid->y_Hy = NULL;
    grid->y_Hz = NULL;
    grid->z_Hx = NULL;
    grid->z_Hy = NULL;
    grid->z_Hz = NULL;
    grid->epsilon_Hx = NULL;
    grid->epsilon_Hy = NULL;
    grid->epsilon_Hz = NULL;
    grid->sigmaH_Hx  = NULL;
    grid->sigmaH_Hy  = NULL;
    grid->sigmaH_Hz  = NULL;
}

void free_fdtd_struct(struct yee_grid *grid) {

    if (grid == NULL)
        return;

    free(grid->x_Ex);
    free(grid->x_Ey);
    free(grid->x_Ez);
    free(grid->y_Ex);
    free(grid->y_Ey);
    free(grid->y_Ez);
    free(grid->z_Ex);
    free(grid->z_Ey);
    free(grid->z_Ez);
    free(grid->epsilon_Ex);
    free(grid->epsilon_Ey);
    free(grid->epsilon_Ez);
    free(grid->sigmaE_Ex);
    free(grid->sigmaE_Ey);
    free(grid->sigmaE_Ez);

    free (grid->x_Hx);
    free (grid->x_Hy);
    free (grid->x_Hz);
    free (grid->y_Hx);
    free (grid->y_Hy);
    free (grid->y_Hz);
    free (grid->z_Hx);
    free (grid->z_Hy);
    free (grid->z_Hz);
    free(grid->epsilon_Hx);
    free(grid->epsilon_Hy);
    free(grid->epsilon_Hz);
    free(grid->sigmaH_Hx);
    free(grid->sigmaH_Hy);
    free(grid->sigmaH_Hz);
}

static float *read_floats_alloc(
        char *fnbuf, size_t fnbuf_size,
        char *fnformat, char *identifier, size_t dim) {

    float *result = NULL;
    int wcount;
    size_t actually_read;

    wcount = snprintf(fnbuf, fnbuf_size, fnformat, identifier);
    if (wcount < 0 || wcount > fnbuf_size)
        goto e_out;

    result = malloc(sizeof(float) * dim);
    if (result == NULL)
        goto e_out;

    if (!read_values_float(fnbuf, result, dim, &actually_read))
        goto e_out;

    if (actually_read != dim)
        goto e_out;

    return result;

e_out:
    free(result);
    return NULL;
}

struct yee_grid *read_fdtd_field(char *fnformat) {
    size_t actually_read;
    size_t dims[18];
    char fnbuf[FDTD_FILENAME_BUF_SIZE];
    int const fnbuf_size = FDTD_FILENAME_BUF_SIZE;

    struct yee_grid *result = malloc(sizeof(struct yee_grid));
    init_fdtd_struct(result);

    int wcount = snprintf(&fnbuf[0], fnbuf_size, fnformat, "dims");
    if (wcount < 0 || wcount > fnbuf_size)
        goto e_out;

    if (!read_values_sizet(fnbuf, &dims[0], 18, &actually_read))
        goto e_out;

    if (actually_read != 18)
        goto e_out;

    result->nx_Ex = dims[0];
    result->nx_Ey = dims[1];
    result->nx_Ez = dims[2];
    result->ny_Ex = dims[3];
    result->ny_Ey = dims[4];
    result->ny_Ez = dims[5];
    result->nz_Ex = dims[6];
    result->nz_Ey = dims[7];
    result->nz_Ez = dims[8];
    result->nx_Hx = dims[9];
    result->nx_Hy = dims[10];
    result->nx_Hz = dims[11];
    result->ny_Hx = dims[12];
    result->ny_Hy = dims[13];
    result->ny_Hz = dims[14];
    result->nz_Hx = dims[15];
    result->nz_Hy = dims[16];
    result->nz_Hz = dims[17];

    result->x_Ex = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "x_Ex", result->nx_Ex);
    if (!result->x_Ex)
        goto e_out;

    result->x_Ey = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "x_Ey", result->nx_Ey);
    if (!result->x_Ey)
        goto e_out;

    result->x_Ez = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "x_Ez", result->nx_Ez);
    if (!result->x_Ez)
        goto e_out;

    result->y_Ex = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "y_Ex", result->ny_Ex);
    if (!result->y_Ex)
        goto e_out;

    result->y_Ey = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "y_Ey", result->ny_Ey);
    if (!result->y_Ey)
        goto e_out;

    result->y_Ez = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "y_Ez", result->ny_Ez);
    if (!result->y_Ez)
        goto e_out;

    result->z_Ex = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "z_Ex", result->nz_Ex);
    if (!result->z_Ex)
        goto e_out;

    result->z_Ey = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "z_Ey", result->nz_Ey);
    if (!result->z_Ey)
        goto e_out;

    result->z_Ez = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "z_Ez", result->nz_Ez);
    if (!result->z_Ez)
        goto e_out;

    result->epsilon_Ex = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "epsilon_Ex",
                result->nx_Ex * result->ny_Ex * result->nz_Ex);
    if (!result->epsilon_Ex)
        goto e_out;

    result->epsilon_Ey = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "epsilon_Ey",
                result->nx_Ey * result->ny_Ey * result->nz_Ey);
    if (!result->epsilon_Ey)
        goto e_out;

    result->epsilon_Ez = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "epsilon_Ez",
                result->nx_Ez * result->ny_Ez * result->nz_Ez);
    if (!result->epsilon_Ez)
        goto e_out;

    result->sigmaE_Ex = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "sigmaE_Ex",
                result->nx_Ex * result->ny_Ex * result->nz_Ex);
    if (!result->sigmaE_Ex)
        goto e_out;

    result->sigmaE_Ey = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "sigmaE_Ey",
                result->nx_Ey * result->ny_Ey * result->nz_Ey);
    if (!result->sigmaE_Ey)
        goto e_out;

    result->sigmaE_Ez = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "sigmaE_Ez",
                result->nx_Ez * result->ny_Ez * result->nz_Ez);
    if (!result->sigmaE_Ez)
        goto e_out;

    //
    result->x_Hx = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "x_Hx", result->nx_Hx);
    if (!result->x_Hx)
        goto e_out;

    result->x_Hy = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "x_Hy", result->nx_Hy);
    if (!result->x_Hy)
        goto e_out;

    result->x_Hz = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "x_Hz", result->nx_Hz);
    if (!result->x_Hz)
        goto e_out;

    result->y_Hx = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "y_Hx", result->ny_Hx);
    if (!result->y_Hx)
        goto e_out;

    result->y_Hy = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "y_Hy", result->ny_Hy);
    if (!result->y_Hy)
        goto e_out;

    result->y_Hz = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "y_Hz", result->ny_Hz);
    if (!result->y_Hz)
        goto e_out;

    result->z_Hx = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "z_Hx", result->nz_Hx);
    if (!result->z_Hx)
        goto e_out;

    result->z_Hy = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "z_Hy", result->nz_Hy);
    if (!result->z_Hy)
        goto e_out;

    result->z_Hz = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "z_Hz", result->nz_Hz);
    if (!result->z_Hz)
        goto e_out;

    result->epsilon_Hx = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "epsilon_Hx",
                result->nx_Hx * result->ny_Hx * result->nz_Hx);
    if (!result->epsilon_Hx)
        goto e_out;

    result->epsilon_Hy = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "epsilon_Hy",
                result->nx_Hy * result->ny_Hy * result->nz_Hy);
    if (!result->epsilon_Hy)
        goto e_out;

    result->epsilon_Hz = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "epsilon_Hz",
                result->nx_Hz * result->ny_Hz * result->nz_Hz);
    if (!result->epsilon_Hz)
        goto e_out;

    result->sigmaH_Hx = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "sigmaH_Hx",
                result->nx_Hx * result->ny_Hx * result->nz_Hx);
    if (!result->sigmaH_Hx)
        goto e_out;

    result->sigmaH_Hy = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "sigmaH_Hy",
                result->nx_Hy * result->ny_Hy * result->nz_Hy);
    if (!result->sigmaH_Hy)
        goto e_out;

    result->sigmaH_Hz = read_floats_alloc(
                fnbuf, fnbuf_size, fnformat, "sigmaH_Hz",
                result->nx_Hz * result->ny_Hz * result->nz_Hz);
    if (!result->sigmaH_Hz)
        goto e_out;

    return result;

e_out:

    free_fdtd_struct(result);
    free(result);
    return NULL;
}
