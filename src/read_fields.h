#ifndef READ_FIELDS_H_
#define READ_FIELDS_H_

static void free_fdtd_struct(struct fdtd_fields *fields);

struct fdtd_fields *read_fdtd_field(char *filename_format);

#endif
