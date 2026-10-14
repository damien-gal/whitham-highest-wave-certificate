.DEFAULT_GOAL := all
CC = gcc
FLINT_PREFIX ?= /usr/local
CPPFLAGS += -include precision.h -I$(FLINT_PREFIX)/include/flint -I$(FLINT_PREFIX)/include -Icore -Iconvexity -Iuniqueness
CFLAGS = -O2 -g -std=c11 -Wall -Wextra -Wno-unused-parameter -MMD -MP -MF $@.d
LDLIBS = -L$(FLINT_PREFIX)/lib -Wl,-rpath,$(FLINT_PREFIX)/lib -lflint -lmpfr -lgmp -lm

# Libraries: the sources of core/, and the non-program sources of convexity/
# and uniqueness/. Each is compiled once.
CORE = $(addprefix build/core_,$(addsuffix .o,basic_functions complex_functions constant_functions clausen_cache ball))
CUSP = build/convexity_cusp.o
UNIQUENESS = build/uniqueness_kernel.o build/uniqueness_refine.o
build/core_%.o: core/%.c
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@
build/convexity_%.o: convexity/%.c
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@
build/uniqueness_%.o: uniqueness/%.c
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# Programs: one per source file, named after it.
EXISTENCE = defect norm_t0 scalars supplemental_scalars test_evaluators
CONVEXITY = derivative_cusp derivative_compact local_t1 local_t2 operators closure
UNIQUENESS_PROGRAMS = sign matrix comparison test_kernel test_refine
PROGRAMS = $(addprefix build/,$(EXISTENCE) $(CONVEXITY) $(UNIQUENESS_PROGRAMS))
$(addprefix build/,$(EXISTENCE)): build/%: existence/%.c $(CORE)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(CORE) $(LDLIBS) -o $@
$(addprefix build/,$(CONVEXITY)): build/%: convexity/%.c $(CUSP) $(CORE)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(CUSP) $(CORE) $(LDLIBS) -o $@
$(addprefix build/,$(UNIQUENESS_PROGRAMS)): build/%: uniqueness/%.c $(UNIQUENESS) $(CORE)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(UNIQUENESS) $(CORE) $(LDLIBS) -o $@

.PHONY: all verify check clean
all: $(PROGRAMS)

# A precision or header change rebuilds everything.
$(CORE) $(CUSP) $(UNIQUENESS) $(PROGRAMS): precision.h $(wildcard core/*.h convexity/*.h uniqueness/*.h)
-include $(wildcard build/*.d)

verify:
	python3 run.py
check:
	python3 run.py --check-existing
clean:
	rm -rf build
