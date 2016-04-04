#ifndef READ_VALUES_H_
#define READ_VALUES_H_

#include <stddef.h>
#include <stdbool.h>

bool read_values(
		char const *filename,
		char const *format,
		void       *out_values,
                size_t      values_item_size,
		size_t      values_count,
		size_t     *out_really_read);

bool read_values_sizet(
		char const *filename,
		size_t     *out_values,
		size_t      values_count,
		size_t     *out_really_read);

bool read_values_float(
		char const *filename,
		float      *out_values,
		size_t      values_count,
		size_t     *out_really_read);

#endif
