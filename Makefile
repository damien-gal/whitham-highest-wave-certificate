.DEFAULT_GOAL := all
CC = gcc
FLINT_PREFIX ?= /usr/local
CPPFLAGS += -include precision.h -I$(FLINT_PREFIX)/include/flint -I$(FLINT_PREFIX)/include -Icore -Isrc
CFLAGS = -O2 -g -std=c11 -Wall -Wextra -Wno-unused-parameter -MMD -MP -MF $@.d
LDLIBS = -L$(FLINT_PREFIX)/lib -Wl,-rpath,$(FLINT_PREFIX)/lib -lflint -lmpfr -lgmp -lm

# Shared evaluators (core/) are compiled once and linked into every program
# that needs them.
CORE = $(addprefix build/core_,$(addsuffix .o,basic_functions complex_functions constant_functions clausen_cache))
build/core_%.o: core/%.c
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# Existence programs (approximate solution, defect, T0 norm, scalar checks).
build/verify build/norm build/scalars build/test_evaluators: build/%: existence/%.c $(CORE)
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(CORE) $(LDLIBS) -o $@
build/supplemental_scalars: src/supplemental_scalars.c $(CORE)
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(CORE) $(LDLIBS) -o $@

# Convexity programs (derivative and curvature bounds, operator norms).
build/derivative_cusp build/derivative_bounds: build/%: convexity/%.c $(CORE)
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -Iconvexity $< $(CORE) $(LDLIBS) -o $@
build/local_t1: convexity/root_local_t1.c convexity/derivative_cusp.c $(CORE)
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -DDERIVATIVE_CUSP_LIBRARY -Iconvexity convexity/root_local_t1.c convexity/derivative_cusp.c $(CORE) $(LDLIBS) -o $@
build/local_t2: convexity/local_t2.c convexity/derivative_cusp.c $(CORE)
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -DDERIVATIVE_CUSP_LIBRARY -Iconvexity convexity/local_t2.c convexity/derivative_cusp.c $(CORE) $(LDLIBS) -o $@
build/operators: convexity/derivative_operators.c convexity/derivative_profile.c $(CORE)
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -Iconvexity convexity/derivative_operators.c convexity/derivative_profile.c $(CORE) $(LDLIBS) -o $@

# Uniqueness programs (kernel sign, moments, exact refinement, comparison).
build/matrix: src/matrix_driver.c build/core_basic_functions.o build/core_clausen_cache.o
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) src/matrix_driver.c build/core_basic_functions.o build/core_clausen_cache.o $(LDLIBS) -o $@
build/sign: src/sign.c src/kernel.c src/ball.c
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) src/sign.c src/kernel.c src/ball.c $(LDLIBS) -o $@
build/uniqueness: src/uniqueness.c src/kernel.c src/ball.c src/refine.c
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) src/uniqueness.c src/kernel.c src/ball.c src/refine.c $(LDLIBS) -o $@
build/closure: src/closure.c src/ball.c
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) src/closure.c src/ball.c $(LDLIBS) -o $@
build/test_refine: src/test_refine.c src/refine.c
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LDLIBS) -o $@
build/test_kernel: src/test_kernel.c src/kernel.c src/ball.c
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) src/test_kernel.c src/kernel.c src/ball.c $(LDLIBS) -o $@

PROGRAMS = build/verify build/norm build/scalars build/test_evaluators build/supplemental_scalars \
	build/derivative_cusp build/derivative_bounds build/local_t1 build/local_t2 build/operators \
	build/matrix build/sign build/uniqueness build/closure build/test_refine build/test_kernel

.PHONY: all verify check clean
all: $(PROGRAMS)

# Precision and header changes invalidate every object and program.
$(CORE) $(PROGRAMS): precision.h $(wildcard core/*.h) $(wildcard src/*.h) $(wildcard convexity/*.h)
-include $(wildcard build/*.d)

verify:
	python3 run.py
check:
	python3 run.py --check-existing
clean:
	rm -rf build
