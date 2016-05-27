
__kernel void calcE(int nx, int ny, int nz, float delta_x, float delta_y, float delta_z,
    __global float *Ex,
    __global float *Ey,
    __global float *Ez,
    __global float *Hx,
    __global float *Hy,
    __global float *Hz,
    __global float *C_Ex,
    __global float *C_Ey,
    __global float *C_Ez,
    __global float *D_Ex,
    __global float *D_Ey,
    __global float *D_Ez) {

    int ix = get_global_id(0) + 1;
    int iy = get_global_id(1) + 1;
    int iz = get_global_id(2) + 1;

    int idx = ix * ny * nz + iy * nz + iz;

    int idx010 = ix * ny * nz + (iy - 1) * nz + iz;
    int idx001 = ix * ny * nz + iy * nz + (iz - 1);
    int idx100 = (ix - 1) * ny * nz + iy * nz + iz;

        Ex[idx] = C_Ex[idx] * Ex[idx] + D_Ex[idx] * ((Hz[idx] - Hz[idx010]) / delta_y -
                                             (Hy[idx] - Hy[idx001]) / delta_z);

        Ey[idx] = C_Ey[idx] * Ey[idx] + D_Ey[idx] * ((Hx[idx] - Hx[idx001]) / delta_z -
                                                 (Hz[idx] - Hz[idx100]) / delta_x);

        Ez[idx] = C_Ez[idx] * Ez[idx] + D_Ez[idx] * ((Hy[idx] - Hy[idx100]) / delta_x -
                                                (Hx[idx] - Hx[idx010]) / delta_y);

}

__kernel void calcH(int nx, int ny, int nz, float delta_x, float delta_y, float delta_z,
    __global float *Ex,
    __global float *Ey,
    __global float *Ez,
    __global float *Hx,
    __global float *Hy,
    __global float *Hz,
    __global float *D_Hx,
    __global float *D_Hy,
    __global float *D_Hz) {

    int ix = get_global_id(0);
    int iy = get_global_id(1);
    int iz = get_global_id(2);

    int idx = ix * ny * nz + iy * nz + iz;

    int idx010 = ix * ny * nz + (iy + 1) * nz + iz;
    int idx001 = ix * ny * nz + iy * nz + (iz + 1);
    int idx100 = (ix + 1) * ny * nz + iy * nz + iz;

        Hx[idx] -= D_Hx[idx] * ((Ez[idx010] - Ez[idx]) / delta_y -
                            (Ey[idx001] - Ey[idx]) / delta_z);

        Hy[idx] -= D_Hy[idx] * ((Ex[idx001] - Ex[idx]) / delta_z -
                            (Ez[idx100] - Ez[idx]) / delta_x);


        Hz[idx] -= D_Hz[idx] * ((Ey[idx100] - Ey[idx]) / delta_x -
                            (Ex[idx010] - Ex[idx]) / delta_y);



}

__kernel void rescueField(int ix, int iy, int iz,
                          int nx, int ny, int nz,
                          __global float *Ez,
                          __global float *prevE) {
    int idx = ix * ny * nz + iy * nz + iz;
    *prevE = Ez[idx];
}

__kernel void updateField(int ix, int iy, int iz,
                          int nx, int ny, int nz,
                          __global float *Ez,
                          __global float *Hx,
                          __global float *Hy,
                          __global float *prevE,
                          float C, float D, float R,
                          float voltage,
                          float delta_x,
                          float delta_y) {
    int idx = ix * ny * nz + iy * nz + iz;

    int idx010 = ix * ny * nz + (iy - 1) * nz + iz;
    int idx100 = (ix - 1) * ny * nz + iy * nz + iz;

    Ez[idx] = C * (*prevE) + D * ((Hy[idx] - Hy[idx100]) / delta_x -
                                     (Hx[idx] - Hx[idx010]) / delta_y +
                                     voltage / (R * delta_x * delta_y));
}
