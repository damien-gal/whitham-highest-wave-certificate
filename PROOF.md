# C-port proof contract and review map

AI-written implementation/notes, 10 October 2026, authorized by Damien's request
to unify the CAP in C and use a common lower precision. No production wording
is authorized by this code request. The pending evaluator-text proposal must
be revised after the numerical wave, not applied with its old mixed precisions.

## Retained construction and convexity

The files in `existence/` and `convexity/` are identified in the import manifest.
Their mathematical algorithms and exact approximation are retained. The only
source edits on import replace hardcoded working precision by CAP_PREC. The
unused old drivers are omitted, not newly validated. The analytic contracts are
`archive/existence/PROOF.md` and the derivative/compact/local operator and
assembly proofs in `review/convexity-wave-2026-10-10/` (paths as of the
2026-10-12 restructuring; older records cite `trusted/` and `audit/`).

At 128 bits a passed interval comparison remains a rigorous comparison; an
insufficient enclosure triggers subdivision or rejection. Keeping the previous
truncations is safe because their tails are enclosed, not neglected. Every
component is rerun; old PASS values are not promoted to new precision.

The exact targets are unchanged: A=2.7, delta0=1e-7, q0=q1=q2=.998, r=6e-5;
delta1=5e-7, delta2=1e-5, alpha1=.32, alpha2=1.383, alphaf=.414, curvature .089.
`src/closure.c` ports the derivative-radius assembly from Python Fraction to
Arb enclosures at the same working precision. Strict positive margins and the
printed outward radius bounds are checked in C. It is conditional on the other
component passes; the runner supplies their conjunction.

The supplemental scalar driver also checks the retained .1 hat bounds and
the printed .142 local T0 criterion at the same 128-bit precision.

## Ported kernel bounds

`src/kernel.c` ports the finite Clausen/Fourier sums and analytic tails in
`archive/cap/src/whitham_cap/validated_kernel.py`. It retains direct interval formulas;
polynomial evaluation uses a centered Taylor shift, with interval coefficients.
The shift need not be exact for validity: every coefficient is enclosed.
`src/sign.c` verifies the 26 exact printed root-polynomial coefficients and all
seven targets on 64 cusp cells and a recursively subdivided 384-cell compact
cover ending at 22/7 > pi. Both children must succeed after every split.
There is no sampled root interpolation. The normalized cusp estimates and
mean-value compact estimates are the ones in
`cap/obligations/section-4-kernel-sign-proof.md`.

`src/ball.c` is a single-thread scratch arena for Arb formula evaluation. Each
operation returns a fresh initialized ball. Pointers remain stable; a reset
releases only later temporary values. Persistent values are copied to owned
Arb objects before reset. No arithmetic is performed in Python, no ball is
silently replaced by its midpoint for proof comparisons, and unsupported or
nonfinite operations fail subsequent required strict tests. The array limits
fail explicitly. Lattice conversion takes a directed Arb endpoint, scales by
an exact power of two, and uses directed integer rounding with a fit check.

The three-image second-derivative remainder retains the old subtraction formula
for its 48 finite coefficients and the same geometric tail. The 96-bit failed
positivity test is retained as a precision diagnostic. Rewriting that expression
could potentially allow fewer bits but is unnecessary to obtain a uniform,
successful 128-bit proof path.

## Ported uniqueness

`src/matrix_driver.c` is the maintained matrix driver with only working precision
changed. `src/refine.c` and its exhaustive small-domain test are unchanged.
`src/uniqueness.c` ports both Python validators. It reads fresh outward moments,
checks their sizes/order/conservation, replaces first-cell moments with weighted
bounds, verifies universal upper/lower seeds (the universal W coefficient, the
constant A of Section 8, is asserted below 0.855 as stated there and in Appendix C,
besides below 1; added 2026-10-12), runs 240 coupled cusp/nodal updates,
and validates the cusp comparison plus every one of the 719 compact intervals.
The same grid-index selection is used in matrix generation and consumption.
That selection chooses a valid pivot node; no binary-double comparison decides
a mathematical inequality. The final grid point is exact pi in both stages.

The mathematical proof is already integrated in §§8–9 and Appendix C; the frozen
side proof is `archive/uniqueness/side-proof.tex`. In particular this port retains:

- the monotone lower envelope at the cusp, without presuming u/sqrt(x) monotonic;
- the continuous-versus-lattice minimization loss;
- the strict whole-interval exclusion test, not unsupported root selection;
- checked signed-128-bit intermediates;
- the regular remainders M<.083 and M3<.067;
- the constructed-wave error 6e-5*x in comparison denominators;
- the derivative cover between grid nodes and the entire punctured cusp cell.

The imported integer proof implementation is not made more rigorous by a C
port; its preexisting analytic justification and independent human review remain
necessary. Regression checks are supplementary, not premises of the proof.

## Validation and provenance

Before a release PASS, run all uniform components and deliberate impossible
bounds, the 300 exhaustive discrete minimization tests, and the C helper/kernel
checks. Differential diagnostics compare 20 C kernel values/derivatives with
the old Python implementation, and compare every integer nodal enclosure when
the C port is supplied with the historical moment file. These old-data checks
are explicitly diagnostic; the actual certificate uses freshly generated
128-bit matrices and a fresh 128-bit kernel-sign verification.

The runner accepts a conjunction of exact stage contracts and binds local
sources, inputs, executable binaries, logs, arithmetic headers/libraries and
regenerated moments. It checks current manuscript root coefficients and targets
separately from recording the whole-file hash for historical traceability.
Failure before/during compilation or numerical work invalidates the current
attempt; there is no fixed old existence-certificate dependency. This repairs
the previous package's two identified replay blockers for the new entrypoint.
The historical runners/artifacts remain frozen, with their original limitations.
