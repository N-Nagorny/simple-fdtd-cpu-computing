#include <errno.h>
#include <stdio.h>
#include "read_values.h"

bool read_values_sizet(
		char const *filename,
		size_t     *out_values,
		size_t      values_count,
		size_t     *out_really_read) {

	return read_values(filename, "%zu",
			out_values, sizeof(*out_values),
			values_count, out_really_read);
}

bool read_values_float(
		char const *filename,
		float      *out_values,
		size_t      values_count,
		size_t     *out_really_read) {

	return read_values(filename, "%f",
			out_values, sizeof(*out_values),
			values_count, out_really_read);
}

bool read_values(
		char const *filename,
		char const *format,
		void       *out_values,
                size_t      values_item_size,
		size_t      values_count,
		size_t     *out_really_read) {

	bool success = true;
	FILE *fp = NULL;

	fp = fopen(filename, "r");
	if (fp == NULL)
		return false;

	errno = 0;

	char *out_values_bytes = out_values;

	int i = 0;
	int res;
	do {
		if (i >= values_count) {
			success = true;
			goto e_out;
		}

		res = fscanf(fp, format, &out_values_bytes[i*values_item_size]);
		if (res == EOF) {
			success = (errno == 0);
			goto e_out;
		}

		if (res == 0) {
			success = false;
			goto e_out;
		}

		++i;
	} while (1);

e_out:
	if (fp != NULL) {
		int res = fclose(fp);
		if (res != 0)
			success = false;
	}

	*out_really_read = i;
	return success;
}

