#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "read_values.h"

static void test_read_ints();
static void test_read_floats();

int main() {
	test_read_ints();
	test_read_floats();
}

static void test_read_ints() {
	size_t buf[10];
	size_t n = 10;
	size_t actual_n = 0;

        bool okay = read_values_sizet(
		"test_read_values_ints.txt", &buf[0], n, &actual_n);

	assert(okay);
	assert(actual_n == 4);
	assert(buf[0] == 100500);
	assert(buf[1] == 100501);
	assert(buf[2] == 100502);
	assert(buf[3] == 100542);
}

static void test_read_floats() {
	const float eps = 10e-5;
	float buf[10];
	size_t n = 10;
	size_t actual_n = 0;

        bool okay = read_values_float(
		"test_read_values_floats.txt", &buf[0], n, &actual_n);

	assert(okay);
	assert(actual_n == 4);
	assert(fabs(buf[0] - 100500) < eps);
	assert(fabs(buf[1] - 100501) < eps);
	assert(fabs(buf[2] - 100502) < eps);
	assert(fabs(buf[3] - 100542) < eps);
}

