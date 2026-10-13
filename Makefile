.DEFAULT_GOAL := all
CC = gcc
FLINT_PREFIX ?= /usr/local
CPPFLAGS += -include precision.h -I$(FLINT_PREFIX)/include/flint -I$(FLINT_PREFIX)/include -Isrc
CFLAGS = -O2 -g -std=c11 -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -MMD -MP -MF $@.d
LDLIBS = -L$(FLINT_PREFIX)/lib -Wl,-rpath,$(FLINT_PREFIX)/lib -lflint -lmpfr -lgmp -lm
EC = $(addprefix build/,$(addsuffix .eo,basic_functions complex_functions constant_functions methods clausen_cache))
DC = $(addprefix build/,$(addsuffix .co,basic_functions complex_functions constant_functions methods clausen_cache))
build/%.eo: existence/%.c
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -Iexistence -c $< -o $@
build/%.co: convexity/%.c
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -Iconvexity -c $< -o $@
build/verify build/norm build/scalars build/test_evaluators: build/%: existence/%.c $(EC)
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -Iexistence $(filter %.eo %.co,$^) $< $(LDLIBS) -o $@
build/derivative_cusp build/derivative_bounds: build/%: convexity/%.c $(DC)
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -Iconvexity $(filter %.eo %.co,$^) $< $(LDLIBS) -o $@
build/local_t1: convexity/root_local_t1.c convexity/derivative_cusp.c $(DC)
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -DDERIVATIVE_CUSP_LIBRARY -Iconvexity convexity/root_local_t1.c convexity/derivative_cusp.c $(DC) $(LDLIBS) -o $@
build/local_t2: convexity/local_t2.c convexity/derivative_cusp.c $(DC)
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -DDERIVATIVE_CUSP_LIBRARY -Iconvexity convexity/local_t2.c convexity/derivative_cusp.c $(DC) $(LDLIBS) -o $@
build/operators: convexity/derivative_operators.c convexity/derivative_profile.c $(DC)
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -Iconvexity convexity/derivative_operators.c convexity/derivative_profile.c $(DC) $(LDLIBS) -o $@
build/matrix: src/matrix_driver.c build/basic_functions.eo build/clausen_cache.eo
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -Iexistence src/matrix_driver.c build/basic_functions.eo build/clausen_cache.eo $(LDLIBS) -o $@
build/sign: src/sign.c src/kernel.c src/ball.c
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) src/sign.c src/kernel.c src/ball.c $(LDLIBS) -o $@
build/test_refine: src/test_refine.c src/refine.c
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LDLIBS) -o $@
.PHONY: all verify check clean
all: build/verify build/norm build/scalars build/test_evaluators build/derivative_cusp build/derivative_bounds build/local_t1 build/local_t2 build/operators build/matrix build/sign build/test_refine build/uniqueness build/closure build/test_kernel
build/uniqueness: src/uniqueness.c src/kernel.c src/ball.c src/refine.c
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) src/uniqueness.c src/kernel.c src/ball.c src/refine.c $(LDLIBS) -o $@
build/closure: src/closure.c src/ball.c
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) src/closure.c src/ball.c $(LDLIBS) -o $@
build/test_kernel: src/test_kernel.c src/kernel.c src/ball.c
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) src/test_kernel.c src/kernel.c src/ball.c $(LDLIBS) -o $@

# Precision/header changes invalidate all proof executables and objects.
$(EC) $(DC): precision.h
build/verify build/norm build/scalars build/test_evaluators build/derivative_cusp build/derivative_bounds build/local_t1 build/local_t2 build/operators build/matrix build/sign build/test_refine build/uniqueness build/closure build/test_kernel: precision.h
-include $(wildcard build/*.d)
verify:
	python3 run.py
check:
	python3 run.py --check-existing
clean:
	rm -rf build

$(EC) $(DC): $(wildcard existence/*.h convexity/*.h)
build/sign build/test_kernel build/uniqueness build/closure: $(wildcard src/*.h)

all: build/supplemental_scalars
build/supplemental_scalars: src/supplemental_scalars.c $(EC) precision.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -Iexistence src/supplemental_scalars.c $(EC) $(LDLIBS) -o $@
