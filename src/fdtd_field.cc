#include "fdtd_field.h"
#include "yee_grid.h"
#include <stdlib.h>

struct fdtd_field *init_grid_alloc (struct yee_grid *grid) {

    struct fdtd_field *field = NULL;
    field = malloc(sizeof(struct fdtd_field));
    if (!field)
        goto e_out;

    field->Hx = NULL;
    field->Hy = NULL;
    field->Hz = NULL;

    field->Ex = NULL;
    field->Ey = NULL;
    field->Ez = NULL;

    field->D_Hx = NULL;
    field->D_Hy = NULL;
    field->D_Hz = NULL;

    field->C_Ex = NULL;
    field->C_Ey = NULL;
    field->C_Ez = NULL;

    field->D_Ex = NULL;
    field->D_Ey = NULL;
    field->D_Ez = NULL;

    field->Hx = malloc((grid->nx_Hx + 1) * (grid->ny_Hx + 1) * (grid->nz_Hx + 1);
    if (!field->Hx)
        goto e_out;
    field->Hy = malloc((grid->nx_Hy + 1) * (grid->ny_Hy + 1) * (grid->nz_Hy + 1);
    if (!field->Hy)
        goto e_out;
    field->Hz = malloc((grid->nx_Hz + 1) * (grid->ny_Hz + 1) * (grid->nz_Hz + 1);
    if (!field->Hz)
        goto e_out;

    field->Ex = malloc((grid->nx_Ex + 1) * (grid->ny_Ex + 1) * (grid->nz_Ex + 1);
    if (!field->Ex)
        goto e_out;
    field->Ey = malloc((grid->nx_Ey + 1) * (grid->ny_Ey + 1) * (grid->nz_Ey + 1);
    if (!field->Ey)
        goto e_out;
    field->Ez = malloc((grid->nx_Ez + 1) * (grid->ny_Ez + 1) * (grid->nz_Ez + 1);
    if (!field->Ez)
        goto e_out;

    field->D_Hx = malloc((grid->nx_Hx + 1) * (grid->ny_Hx + 1) * (grid->nz_Hx + 1);
    if (!field->D_Hx)
        goto e_out;
    field->D_Hy = malloc((grid->nx_Hy + 1) * (grid->ny_Hy + 1) * (grid->nz_Hy + 1);
    if (!field->D_Hy)
        goto e_out;
    field->D_Hz = malloc((grid->nx_Hz + 1) * (grid->ny_Hz + 1) * (grid->nz_Hz + 1);
    if (!field->D_Hz)
        goto e_out;

    field->C_Ex = malloc((grid->nx_Hx + 1) * (grid->ny_Hx + 1) * (grid->nz_Hx + 1);
    if (!field->C_Ex)
        goto e_out;
    field->C_Ey = malloc((grid->nx_Hy + 1) * (grid->ny_Hy + 1) * (grid->nz_Hy + 1);
    if (!field->C_Ey)
        goto e_out;
    field->C_Ez = malloc((grid->nx_Hz + 1) * (grid->ny_Hz + 1) * (grid->nz_Hz + 1);
    if (!field->C_Ez)
        goto e_out;

    field->D_Ex = malloc((grid->nx_Hx + 1) * (grid->ny_Hx + 1) * (grid->nz_Hx + 1);
    if (!field->D_Ex)
        goto e_out;
    field->D_Ey = malloc((grid->nx_Hy + 1) * (grid->ny_Hy + 1) * (grid->nz_Hy + 1);
    if (!field->D_Ey)
        goto e_out;
    field->D_Ez = malloc((grid->nx_Hz + 1) * (grid->ny_Hz + 1) * (grid->nz_Hz + 1);
    if (!field->D_Ez)
        goto e_out;


    field->delta_x_Ex = malloc(grid->nx_Ex + 1);
    if (!field->delta_x_Ex)
        goto e_out;
    field->delta_y_Ex = malloc(grid->ny_Ex + 1);
    if (!field->delta_y_Ex)
        goto e_out;
    field->delta_z_Ex = malloc(grid->nz_Ex + 1);
    if (!field->delta_z_Ex)
        goto e_out;
    field->delta_x_Ey = malloc(grid->nx_Ey + 1);
    if (!field->delta_x_Ey)
        goto e_out;
    field->delta_y_Ey = malloc(grid->ny_Ey + 1);
    if (!field->delta_y_Ey)
        goto e_out;
    field->delta_z_Ey = malloc(grid->nz_Ey + 1);
    if (!field->delta_z_Ey)
        goto e_out;
    field->delta_x_Ez = malloc(grid->nx_Ez + 1);
    if (!field->delta_x_Ez)
        goto e_out;
    field->delta_y_Ez = malloc(grid->ny_Ez + 1);
    if (!field->delta_y_Ez)
        goto e_out;
    field->delta_z_Ez = malloc(grid->nz_Ez + 1);
    if (!field->delta_z_Ez)
        goto e_out;

    field->delta_x_Hx = malloc(grid->nx_Hx + 1);
    if (!field->delta_x_Hx)
        goto e_out;
    field->delta_y_Hx = malloc(grid->ny_Hx + 1);
    if (!field->delta_y_Hx)
        goto e_out;
    field->delta_z_Hx = malloc(grid->nz_Hx + 1);
    if (!field->delta_z_Hx)
        goto e_out;
    field->delta_x_Hy = malloc(grid->nx_Hy + 1);
    if (!field->delta_x_Hy)
        goto e_out;
    field->delta_y_Hy = malloc(grid->ny_Hy + 1);
    if (!field->delta_y_Hy)
        goto e_out;
    field->delta_z_Hy = malloc(grid->nz_Hy + 1);
    if (!field->delta_z_Hy)
        goto e_out;
    field->delta_x_Hz = malloc(grid->nx_Hz + 1);
    if (!field->delta_x_Hz)
        goto e_out;
    field->delta_y_Hz = malloc(grid->ny_Hz + 1);
    if (!field->delta_y_Hz)
        goto e_out;
    field->delta_z_Hz = malloc(grid->nz_Hz + 1);
    if (!field->delta_z_Hz)
        goto e_out;

    return field;

e_out:
    return NULL;

}

void free_grid (struct fdtd_field *field) {

    if (field == NULL)
        return;

    free(field->Hx);
    free(field->Hy);
    free(field->Hz);

    free(field->Ex);
    free(field->Ey);
    free(field->Ez);

    free(field->D_Hx);
    free(field->D_Hy);
    free(field->D_Hz);

    free(field->C_Ex);
    free(field->C_Ey);
    free(field->C_Ez);

    free(field->D_Ex);
    free(field->D_Ey);
    free(field->D_Ez);


    free(field->delta_x_Ex);
    free(field->delta_y_Ex);
    free(field->delta_z_Ex);
    free(field->delta_x_Ey);
    free(field->delta_y_Ey);
    free(field->delta_z_Ey);
    free(field->delta_x_Ez);
    free(field->delta_y_Ez);
    free(field->delta_z_Ez);
    free(field->delta_x_Hx);
    free(field->delta_y_Hx);
    free(field->delta_z_Hx);
    free(field->delta_x_Hy);
    free(field->delta_y_Hy);
    free(field->delta_z_Hy);
    free(field->delta_x_Hz);
    free(field->delta_y_Hz);
    free(field->delta_z_Hz);

}

void calc_coefs(struct yee_grid *grid, struct fdtd_field *field, float delta_t) {
    int ix = 0;
    int iy = 0;
    int iz = 0;


    for (ix = 0; ix < grid->nx_Ex; ix++) {
        field->delta_x_Ex[ix] = grid->x_Ex[ix + 1] - grid->x_Ex[ix];
    }

    for (iy = 0; iy < grid->ny_Ex; iy++) {
        field->delta_y_Ex[iy] = grid->y_Ex[iy + 1] - grid->y_Ex[iy];
    }

    for (iz = 0; iz < grid->nz_Ex; iz++) {
        field->delta_z_Ex[iz] = grid->z_Ex[iz + 1] - grid->z_Ex[iz];
    }

    for (ix = 0; ix < grid->nx_Ey; ix++) {
        field->delta_x_Ey[ix] = grid->x_Ey[ix + 1] - grid->x_Ey[ix];
    }

    for (iy = 0; iy < grid->ny_Ey; iy++) {
        field->delta_y_Ey[iy] = grid->y_Ey[iy + 1] - grid->y_Ey[iy];
    }

    for (iz = 0; iz < grid->nz_Ey; iz++) {
        field->delta_z_Ey[iz] = grid->z_Ey[iz + 1] - grid->z_Ey[iz];
    }

    for (ix = 0; ix < grid->nx_Ez; ix++) {
        field->delta_x_Ez[ix] = grid->x_Ez[ix + 1] - grid->x_Ez[ix];
    }

    for (iy = 0; iy < grid->ny_Ez; iy++) {
        field->delta_y_Ez[iy] = grid->y_Ez[iy + 1] - grid->y_Ez[iy];
    }

    for (iz = 0; iz < grid->nz_Ez; iz++) {
        field->delta_z_Ez[iz] = grid->z_Ez[iz + 1] - grid->z_Ez[iz];
    }

    for (ix = 0; ix < grid->nx_Hx; ix++) {
        field->delta_x_Hx[ix] = grid->x_Hx[ix + 1] - grid->x_Hx[ix];
    }

    for (iy = 0; iy < grid->ny_Hx; iy++) {
        field->delta_y_Hx[iy] = grid->y_Hx[iy + 1] - grid->y_Hx[iy];
    }

    for (iz = 0; iz < grid->nz_Hx; iz++) {
        field->delta_z_Hx[iz] = grid->z_Hx[iz + 1] - grid->z_Hx[iz];
    }

    for (ix = 0; ix < grid->nx_Hy; ix++) {
        field->delta_x_Hy[ix] = grid->x_Hy[ix + 1] - grid->x_Hy[ix];
    }

    for (iy = 0; iy < grid->ny_Hy; iy++) {
        field->delta_y_Hy[iy] = grid->y_Hy[iy + 1] - grid->y_Hy[iy];
    }

    for (iz = 0; iz < grid->nz_Hy; iz++) {
        field->delta_z_Hy[iz] = grid->z_Hy[iz + 1] - grid->z_Hy[iz];
    }

    for (ix = 0; ix < grid->nx_Hz; ix++) {
        field->delta_x_Hz[ix] = grid->x_Hz[ix + 1] - grid->x_Hz[ix];
    }

    for (iy = 0; iy < grid->ny_Hz; iy++) {
        field->delta_y_Hz[iy] = grid->y_Hz[iy + 1] - grid->y_Hz[iy];
    }

    for (iz = 0; iz < grid->nz_Hz; iz++) {
        field->delta_z_Hz[iz] = grid->z_Hz[iz + 1] - grid->z_Hz[iz];
    }

    for (ix = 0; ix < grid->nx_Hx; ix++)
    for (iy = 0; iy < grid->ny_Hx; iy++)
    for (iz = 0; iz < grid->nz_Hx; iz++) {
        size_t idx    = ix * grid->nx_Hx * grid->ny_Hx + iy * grid->nz_Hx + iz;
        float epsilon = grid->epsilon_Hx[idx];
        float sigma   = grid->sigmaH_Hx[idx];
        field->D_Hx[idx] = - (delta_t / epsilon) /
                (1 + sigma * delta_t / (2 * epsilon));

        field->Hx[idx] = 0;
    }

    for (ix = 0; ix < grid->nx_Hy; ix++)
    for (iy = 0; iy < grid->ny_Hy; iy++)
    for (iz = 0; iz < grid->nz_Hy; iz++) {
        size_t idx    = ix * grid->nx_Hx * grid->ny_Hx + iy * grid->nz_Hx + iz;
        float epsilon = grid->epsilon_Hy[idx];
        float sigma   = grid->sigmaH_Hy[idx];
        field->D_Hy[idx] = - (delta_t / epsilon) /
                (1 + sigma * delta_t / (2 * epsilon));
        field->Hy[idx] = 0;
    }

    for (ix = 0; ix < grid->nx_Hz; ix++)
    for (iy = 0; iy < grid->ny_Hz; iy++)
    for (iz = 0; iz < grid->nz_Hz; iz++) {
        size_t idx    = ix * grid->nx_Hx * grid->ny_Hx + iy * grid->nz_Hx + iz;
        float epsilon = grid->epsilon_Hz[idx];
        float sigma   = grid->sigmaH_Hz[idx];
        field->D_Hz[idx] = - (delta_t / epsilon) /
                (1 + sigma * delta_t / (2 * epsilon));
        field->Hz[idx] = 0;
    }

    for (ix = 0; ix < grid->nx_Ex; ix++)
    for (iy = 0; iy < grid->ny_Ex; iy++)
    for (iz = 0; iz < grid->nz_Ex; iz++) {
        size_t idx    = ix * grid->nx_Hx * grid->ny_Hx + iy * grid->nz_Hx + iz;
        float epsilon = grid->epsilon_Ex[idx];
        float sigma   = grid->sigmaE_Ex[idx];
        field->C_Ex[idx] = (1 - sigma * delta_t / (2 * epsilon)) /
                (1 + sigma * delta_t / (2 * epsilon));
        field->D_Ex[idx] = - (delta_t / epsilon) /
                (1 + sigma * delta_t / (2 * epsilon));
        field->Ex[idx] = 0;
    }

    for (ix = 0; ix < grid->nx_Ey; ix++)
    for (iy = 0; iy < grid->ny_Ey; iy++)
    for (iz = 0; iz < grid->nz_Ey; iz++) {
        size_t idx    = ix * grid->nx_Hx * grid->ny_Hx + iy * grid->nz_Hx + iz;
        float epsilon = grid->epsilon_Ey[idx];
        float sigma   = grid->sigmaE_Ey[idx];
        field->C_Ey[idx] = (1 - sigma * delta_t / (2 * epsilon)) /
                (1 + sigma * delta_t / (2 * epsilon));
        field->D_Ey[idx] = - (delta_t / epsilon) /
                (1 + sigma * delta_t / (2 * epsilon));
        field->Ey[idx] = 0;
    }

    for (ix = 0; ix < grid->nx_Ez; ix++)
    for (iy = 0; iy < grid->ny_Ez; iy++)
    for (iz = 0; iz < grid->nz_Ez; iz++) {
        size_t idx    = ix * grid->nx_Hx * grid->ny_Hx + iy * grid->nz_Hx + iz;
        float epsilon = grid->epsilon_Ez[idx];
        float sigma   = grid->sigmaE_Ez[idx];
        field->C_Ez[idx] = (1 - sigma * delta_t / (2 * epsilon)) /
                (1 + sigma * delta_t / (2 * epsilon));
        field->D_Ez[idx] = - (delta_t / epsilon) /
                (1 + sigma * delta_t / (2 * epsilon));
        field->Ez[idx] = 0;
    }

}

void calc_fieldH(struct yee_grid *grid, struct fdtd_field *field) {

    int ix = 0;
    int iy = 0;
    int iz = 0;

    for (ix = 0; ix < grid->nx_Hx; ix++)
    for (iy = 0; iy < grid->ny_Hx; iy++)
    for (iz = 0; iz < grid->nz_Hx; iz++) {

        int idx_Hx = ix * grid->nx_Hx * grid->ny_Hx + iy * grid->nz_Hx + iz;
        int idx_Ey = ix * grid->nx_Ey * grid->ny_Ey + iy * grid->nz_Ey + iz;
        int idx_Ez = ix * grid->nx_Ez * grid->ny_Ez + iy * grid->nz_Ez + iz;

        int idx_Ez1 = ix * grid->nx_Ez * grid->ny_Ez + (iy + 1) * grid->nz_Ez + iz;
        int idx_Ey1 = ix * grid->nx_Ez * grid->ny_Ez + iy * grid->nz_Ez + (iz + 1);

        field->Hx[idx_Hx] += field->D_Hx[idx_Hx] *
                (((field->Ez[idx_Ez1] - field->Ez[idx_Ez]) / field->delta_y_Ez[iy]) -
                 ((field->Ey[idx_Ey1] - field->Ey[idx_Ey]) / field->delta_z_Ey[iz]));
    }

    for (ix = 0; ix < grid->nx_Hy; ix++)
    for (iy = 0; iy < grid->ny_Hy; iy++)
    for (iz = 0; iz < grid->nz_Hy; iz++) {

        int idx_Hy = ix * grid->nx_Hy * grid->ny_Hy + iy * grid->nz_Hy + iz;
        int idx_Ex = ix * grid->nx_Ex * grid->ny_Ex + iy * grid->nz_Ex + iz;
        int idx_Ez = ix * grid->nx_Ez * grid->ny_Ez + iy * grid->nz_Ez + iz;

        int idx_Ez1 = (ix + 1) * grid->nx_Ez * grid->ny_Ez + iy * grid->nz_Ez + iz;
        int idx_Ex1 = ix * grid->nx_Ez * grid->ny_Ez + iy * grid->nz_Ez + (iz + 1);

        field->Hy[idx_Hy] += field->D_Hy[idx_Hy] *
                (((field->Ex[idx_Ex1] - field->Ex[idx_Ex]) / field->delta_z_Ex[iz]) -
                 ((field->Ez[idx_Ez1] - field->Ez[idx_Ez]) / field->delta_x_Ez[ix]));
    }

    for (ix = 0; ix < grid->nx_Hz; ix++)
    for (iy = 0; iy < grid->ny_Hz; iy++)
    for (iz = 0; iz < grid->nz_Hz; iz++) {

        int idx_Hz = ix * grid->nx_Hz * grid->ny_Hz + iy * grid->nz_Hz + iz;
        int idx_Ex = ix * grid->nx_Ex * grid->ny_Ex + iy * grid->nz_Ex + iz;
        int idx_Ey = ix * grid->nx_Ey * grid->ny_Ey + iy * grid->nz_Ey + iz;

        int idx_Ey1 = (ix + 1) * grid->nx_Ey * grid->ny_Ey + iy * grid->nz_Ey + iz;
        int idx_Ex1 = ix * grid->nx_Ex * grid->ny_Ex + (iy + 1) * grid->nz_Ex + iz;

        field->Hz[idx_Hz] += field->D_Hz[idx_Hz] *
                (((field->Ey[idx_Ey1] - field->Ey[idx_Ey]) / field->delta_x_Ey[ix]) -
                 ((field->Ex[idx_Ex1] - field->Ex[idx_Ex]) / field->delta_y_Ex[iy]));
    }
}

void calc_fieldE(struct yee_grid *grid, struct fdtd_field *field) {

    int ix = 0;
    int iy = 0;
    int iz = 0;

    for (ix = 0; ix < grid->nx_Ex; ix++)
    for (iy = 0; iy < grid->ny_Ex; iy++)
    for (iz = 0; iz < grid->nz_Ex; iz++) {

        int idx_Ex = ix * grid->nx_Ex * grid->ny_Ex + iy * grid->nz_Ex + iz;
        int idx_Hy = ix * grid->nx_Hy * grid->ny_Hy + iy * grid->nz_Hy + iz;
        int idx_Hz = ix * grid->nx_Hz * grid->ny_Hz + iy * grid->nz_Hz + iz;

        int idx_Hz1 = ix * grid->nx_Hz * grid->ny_Hz + (iy - 1) * grid->nz_Hz + iz;
        int idx_Hy1 = ix * grid->nx_Hz * grid->ny_Hz + iy * grid->nz_Hz + (iz - 1);

        field->Ex[idx_Ex] += field->C_Ex[idx_Ex] *
                (((field->Hz[idx_Hz] - field->Hz[idx_Hz1]) / field->delta_y_Hz[iy]) -
                 ((field->Hy[idx_Hy] - field->Hy[idx_Hy1]) / field->delta_z_Hy[iz]));
    }

    for (ix = 0; ix < grid->nx_Ey; ix++)
    for (iy = 0; iy < grid->ny_Ey; iy++)
    for (iz = 0; iz < grid->nz_Ey; iz++) {

        int idx_Ey = ix * grid->nx_Ey * grid->ny_Ey + iy * grid->nz_Ey + iz;
        int idx_Hx = ix * grid->nx_Hx * grid->ny_Hx + iy * grid->nz_Hx + iz;
        int idx_Hz = ix * grid->nx_Hz * grid->ny_Hz + iy * grid->nz_Hz + iz;

        int idx_Hz1 = (ix - 1) * grid->nx_Hz * grid->ny_Hz + iy * grid->nz_Hz + iz;
        int idx_Hx1 = ix * grid->nx_Hx * grid->ny_Hx + iy * grid->nz_Hx + (iz - 1);

        field->Ey[idx_Ey] += field->C_Ey[idx_Ey] *
                (((field->Hx[idx_Hx] - field->Hx[idx_Hx1]) / field->delta_z_Hx[iz]) -
                 ((field->Hz[idx_Hz] - field->Hz[idx_Hz1]) / field->delta_x_Hz[ix]));
    }

    for (ix = 0; ix < grid->nx_Ez; ix++)
    for (iy = 0; iy < grid->ny_Ez; iy++)
    for (iz = 0; iz < grid->nz_Ez; iz++) {

        int idx_Ez = ix * grid->nx_Ez * grid->ny_Ez + iy * grid->nz_Ez + iz;
        int idx_Hx = ix * grid->nx_Hx * grid->ny_Hx + iy * grid->nz_Hx + iz;
        int idx_Hy = ix * grid->nx_Hy * grid->ny_Hy + iy * grid->nz_Hy + iz;

        int idx_Hy1 = (ix - 1) * grid->nx_Hy * grid->ny_Hy + iy * grid->nz_Hy + iz;
        int idx_Hx1 = ix * grid->nx_Hx * grid->ny_Hx + (iy - 1) * grid->nz_Hx + iz;

        field->Ez[idx_Ez] += field->C_Ez[idx_Ez] *
                (((field->Hy[idx_Hy] - field->Hy[idx_Hy1]) / field->delta_x_Hy[ix]) -
                 ((field->Hx[idx_Hx] - field->Hx[idx_Hx1]) / field->delta_y_Hx[iy]));
    }
}
