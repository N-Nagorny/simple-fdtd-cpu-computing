#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <errno.h>
#include "fdtd_fields.h"
#include "read_fields.h"

int main(int argc, char *argv[]) {
    struct fdtd_fields *fields = NULL;
    int i;
    int retcode = 1;

    fields = read_fdtd_field("input/filename_%s.txt");
    if (!fields)
        goto e_out;

    for (i = 0; i < fields->nx_Ex; i++) {
        printf("%f\n", fields->x_Ex[i]);
    }

    printf("\n%zu\n", fields->nx_Ex);

    retcode = 0;

e_out:
    free_fdtd_struct(fields);
    free(fields);
    return retcode;
}
