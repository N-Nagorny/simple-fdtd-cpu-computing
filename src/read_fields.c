#include <stdlib.h>
#include <stdio.h>
#include "fdtd_fields.h"

#define FDTD_FILENAME_BUF_SIZE 1024

static void init_fdtd_struct(struct fdtd_fields *fields) {
    fields->x_Ex = NULL;
    fields->x_Ey = NULL;
    fields->x_Ez = NULL;
    fields->y_Ex = NULL;
    fields->y_Ey = NULL;
    fields->y_Ez = NULL;
    fields->z_Ex = NULL;
    fields->z_Ey = NULL;
    fields->z_Ez = NULL;
    fields->epsilon_Ex = NULL;
    fields->epsilon_Ey = NULL;
    fields->epsilon_Ez = NULL;
    fields->sigmaE_Ex  = NULL;
    fields->sigmaE_Ey  = NULL;
    fields->sigmaE_Ez  = NULL;

    fields->x_Hx = NULL;
    fields->x_Hy = NULL;
    fields->x_Hz = NULL;
    fields->y_Hx = NULL;
    fields->y_Hy = NULL;
    fields->y_Hz = NULL;
    fields->z_Hx = NULL;
    fields->z_Hy = NULL;
    fields->z_Hz = NULL;
    fields->epsilon_Hx = NULL;
    fields->epsilon_Hy = NULL;
    fields->epsilon_Hz = NULL;
    fields->sigmaH_Hx  = NULL;
    fields->sigmaH_Hy  = NULL;
    fields->sigmaH_Hz  = NULL;
}

void free_fdtd_struct(struct fdtd_fields *fields) {

    if (fields == NULL)
        return;

    free(fields->x_Ex);
    free(fields->x_Ey);
    free(fields->x_Ez);
    free(fields->y_Ex);
    free(fields->y_Ey);
    free(fields->y_Ez);
    free(fields->z_Ex);
    free(fields->z_Ey);
    free(fields->z_Ez);
    free(fields->epsilon_Ex);
    free(fields->epsilon_Ey);
    free(fields->epsilon_Ez);
    free(fields->sigmaE_Ex);
    free(fields->sigmaE_Ey);
    free(fields->sigmaE_Ez);

    free (fields->x_Hx);
    free (fields->x_Hy);
    free (fields->x_Hz);
    free (fields->y_Hx);
    free (fields->y_Hy);
    free (fields->y_Hz);
    free (fields->z_Hx);
    free (fields->z_Hy);
    free (fields->z_Hz);
    free(fields->epsilon_Hx);
    free(fields->epsilon_Hy);
    free(fields->epsilon_Hz);
    free(fields->sigmaH_Hx);
    free(fields->sigmaH_Hy);
    free(fields->sigmaH_Hz);
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

struct fdtd_fields *read_fdtd_field(char *fnformat) {
    size_t actually_read;
    size_t dims[18];
    char fnbuf[FDTD_FILENAME_BUF_SIZE];
    int const fnbuf_size = FDTD_FILENAME_BUF_SIZE;

    struct fdtd_fields *result = malloc(sizeof(struct fdtd_fields));
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
