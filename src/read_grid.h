#ifndef READ_GRID_H_
#define READ_GRID_H_

void free_fdtd_struct(struct yee_grid *grid);

struct yee_grid *read_fdtd_field(char *filename_format);

#endif
