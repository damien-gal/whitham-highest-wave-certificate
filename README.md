# Computer-assisted proof code: Uniqueness, convexity and sharp asymptotics for Whitham's highest wave

This directory contains the code for the computer-assisted parts of the paper
*Uniqueness, convexity and sharp asymptotics for Whitham's highest wave*. It is
distributed with the paper as arXiv ancillary files; in the authors' working
repository it is the directory `certificate/`, beside the paper source.

## Authors

- Alberto Enciso (ICMAT, Madrid), <aenciso@icmat.es>
- Damien Galant (Brown University), <damien_galant@brown.edu>
- Javier Gómez-Serrano (Brown University), <javier_gomez_serrano@brown.edu>
- Bruno Vergara, <bvergara@ub.edu>

## Dependencies

All proof arithmetic is written in C and uses the ball arithmetic of
[FLINT](https://flintlib.org) version 3 (which contains Arb), at a common
working precision of 128 bits fixed in `precision.h`. The reference
environment is:

- Debian 12 (x86_64), GCC 12.2 (C11 with signed `__int128`), GNU Make;
- FLINT 3.2.1 built from source with its bundled dependencies GMP 6 and MPFR 4,
  installed under `/usr/local`;
- Python 3.11 (standard library only) to orchestrate builds, runs and hashes.

To install FLINT 3.2.1 from source on Debian or Ubuntu:

```bash
sudo apt-get install build-essential libgmp-dev libmpfr-dev
curl -LO https://github.com/flintlib/flint/releases/download/v3.2.1/flint-3.2.1.tar.gz
tar xzf flint-3.2.1.tar.gz && cd flint-3.2.1
./configure --prefix=/usr/local && make -j && sudo make install && sudo ldconfig
```

Another installation prefix is passed as `FLINT_PREFIX=<prefix> make verify`.

## Organization

- `run.py` is the orchestrator: it compiles every program, runs the stages below
  with at most two processes, checks their printed contracts and coverage
  records, binds the exact kernel-sign inputs to Appendix C of the paper, and
  writes `artifacts/certificate.json` together with one log per stage.
- `core/` contains the libraries shared by all programs: the special
  functions, the approximate solution and its derivatives, the cusp
  expansions with their remainders and the scalar constants, and a scratch
  arena for ball formulas.
- `existence/` contains the programs of the existence part: the defect bound
  (`defect.c`), the operator norm of $T_0$ (`norm_t0.c`), the scalar conditions
  (`scalars.c`, `supplemental_scalars.c`) and the evaluator checks.
- `convexity/` contains the programs of the convexity part: the derivative and
  curvature bounds near the cusp (`derivative_cusp.c`) and on the compact
  interval (`derivative_compact.c`), the local $T_1$ and $T_2$ estimates, the
  operator norms of $T_1$ and $T_2$ (`operators.c`) and the scalar assembly of
  the correction bounds (`closure.c`); `cusp.c` is their shared library.
- `uniqueness/` contains the programs of the uniqueness part: the kernel-sign
  certificate (`sign.c`), the kernel moment matrices (`matrix.c`) and the
  a priori enclosures with the comparison (`comparison.c`), with the kernel
  arithmetic (`kernel.c`) and the exact integer refinement (`refine.c`) as
  libraries, plus their regression tests.
- `inputs/` holds the exact inputs: the approximation coefficients
  (`approximation.json`), the kernel-sign polynomials and targets
  (`kernel-sign.json`), the closure constants (`closure.json`) and the 720-node
  refinement grid (`grid720.bin`, binary64: a geometric progression from
  $10^{-3}$ to $\pi$ in double precision). Each file states where it comes from;
  `import-manifest.json` records the provenance of the grid.
- `PROOF.md` maps every stage to the statements of the paper and describes
  the constructions; `import-manifest.json` records the provenance of every
  source file.

The stages and the bounds they assert are those of the summary table of
Appendix B of the paper (its last column names the stage):

| Stage | Asserts |
|---|---|
| `existence` | defect $\|F_0\|_\infty \le 10^{-7}$ on $[0.05,\pi]$, with the cusp bounds of Lemma 5.6 |
| `scalars`, `supplemental_scalars` | $C_B < 0.997362$, the speed and $\nu$ enclosures, the local $T_0$ criterion |
| `T0` | $\|T_0\| \le 0.998$ on the compact cells |
| `derivative_cusp`, `derivative_compact` | $\|F_1\|_\infty$, $\|F_2\|_\infty$, the coefficient bounds $\alpha_1,\alpha_2,\alpha_f$ and the curvature bound |
| `local1`, `local2`, `T1`, `T2` | $\|T_1\|,\|T_2\| \le 0.998$ near the cusp and on the compact cells |
| `closure` | the fixed-point radii and the curvature margin |
| `sign` | the kernel-zero brackets of Lemma C.1 |
| `matrix`, `uniqueness` | the kernel moments, the a priori enclosures of Section 8 and the comparison $Q < 0.982$ |
| `evaluators`, `kernel_tests`, `refinement_tests`, `operator_regression` | regression checks of the evaluators and of the exact refinement |
| `reject_*` | deliberate rejection tests: each must return UNRESOLVED |

Every decimal printed in the paper is asserted by one of these stages.

## Reproducing the proof

From this directory:

```bash
make verify
make check
```

`verify` compiles everything, runs all stages and writes `artifacts/certificate.json`
with `"status": "PASS"` only when every stage, coverage record and rejection test
succeeds; any failure yields `UNRESOLVED`. `check` re-validates an existing
certificate against the present sources, inputs, binaries and logs without
re-running the arithmetic. Both commands read the kernel-sign block of Appendix C
from `../paper/appendices/c-uniqueness-kernel-estimates.tex` when that file is
present (it is, in the repository; in the ancillary package the binding is
reported as absent and the exact targets are taken from `inputs/kernel-sign.json`).

## Runtimes

Wall-clock times of the stages in a complete `verify` run on a single machine:
Intel Core i7-11800H (8 cores, 16 threads), 15 GB RAM, Debian 12, with the
reference environment above. Two stages run concurrently.

| Stage | Wall-clock |
|---|---|
| `existence` | 159 s |
| `derivative_compact` | 137 s |
| `T0` | 21 s |
| `T1` | 30 s |
| `T2` | 317 s |
| `matrix` | 289 s |
| `uniqueness` | 2 s |
| `sign` | 0.5 s |
| all other stages | under 1 s each |
| **Total (two workers)** | **486 s** |

## Provenance

The evaluators of `core/` descend from the authors' original C/Arb
implementation of the existence part (2023-2024): the working precision was
unified, routines that no program calls were removed and the working notes
were turned into documentation. The programs of `existence/`, `convexity/`
and `uniqueness/` were written in 2026 for this package; the kernel-sign
certificate, the kernel arithmetic and the uniqueness comparison port the
authors' earlier Python implementation to C. The code written in 2026 was produced with the
assistance of large language models, as stated in the paper, and was checked
by the authors. `import-manifest.json` lists, for every file, its origin and
the edits made since it entered the package.

## License

This code is released under the MIT License; see [`LICENSE`](LICENSE).

## Technical notes

Every FLINT ball call uses the 128-bit working precision of `precision.h`. The
exact lattice refinement keeps signed 128-bit integers with coefficient scale
`2^-42` and profile scale `2^-32`; the stored binary64 grid and the outward matrix
endpoints are exact dyadic inputs, not floating-point acceptance tests. The
input files require little-endian IEEE binary64, which the uniqueness verifier
checks. A successful run records the compiler, the hashes of the sources,
binaries, inputs, logs and linked libraries, and the exact Arb balls of every
printed quantity (`ARB_DUMP` lines in the logs).

The precision was chosen as follows. The earlier implementations used 160, 192
and 256 bits for headroom; at 96 bits one three-image remainder coefficient
obtained by subtraction could not be proved positive, and at 128 bits every
full-domain cover and inequality passes with the manuscript's targets. This is a
demonstrated sufficient precision, not a minimal one. There is no automatic
precision escalation.

The certificate is a reproducible record of one run. It is not a formal proof
object and does not replace the mathematical review of the paper; its scope is
stated in the `scope`, `human_review` and `remaining` fields of
`artifacts/certificate.json`.
