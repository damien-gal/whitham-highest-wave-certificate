# Proof contract and review map

This file states what each program of the package certifies, under which
analytic premises, and how the runner binds the pieces together. The
mathematical arguments are those of the paper *Uniqueness, convexity and sharp
asymptotics for Whitham's highest wave*; nothing here replaces them, and the
certificate does not replace the mathematical review of the paper.

## Scope and arithmetic

All proof arithmetic is C code calling the ball arithmetic of FLINT (Arb) at
the single working precision `CAP_PREC = 128` bits of `precision.h`, together
with the exact signed 128-bit integer refinement of `uniqueness/refine.c`. The
orchestrator `run.py` uses only the Python standard library: it compiles the
programs, runs the stages, checks their printed contracts and coverage records,
and hashes sources, inputs, binaries, logs and the linked libraries. It never
performs proof arithmetic and never accepts a numerical sample as evidence.

A passed ball comparison at 128 bits is a rigorous comparison. Where a cover
is built, an enclosure that is too wide triggers a subdivision; everywhere
else it yields the outcome UNRESOLVED. Every truncated series is accompanied
by an enclosed remainder; no tail is neglected. The precision was chosen as
follows: at 96 bits one remainder coefficient of the three-image second
derivative bound, obtained by subtraction, could not be proved positive, while
at 128 bits every cover and inequality passes with the targets of the paper.
This is a demonstrated sufficient precision, not a minimal one; there is no
automatic precision escalation.

## Stages

| Stage | Program | Certifies | Paper |
|---|---|---|---|
| `existence` | `existence/defect.c` | the defect bound $\|F_0\|_\infty\le10^{-7}$: cusp expansions on $(0,0.05]$, recursive mean-value cover on $[0.05,\pi]$; the printed constants $\lambda_->0.604$ and the cusp defect bound $1.55\cdot10^{-8}$ | lemma on $F_0$, Appendix B |
| `scalars` | `existence/scalars.c` | $C_B<0.997362$ by a bisection enclosure of $t^*$; the hat constants at $\varepsilon=0.05$; the local $T_0$ criterion with the computed constants; the fixed-point constants $A,\delta,q,r$; $0.9455<\nu<0.9471$ and $0.766<\mu<0.770$ | lemma on $u_0$ near the cusp, existence proposition, wave-speed proposition |
| `supplemental_scalars` | `existence/supplemental_scalars.c` | the printed hat bounds $0.217$, $0.24$, $0.142$ for $M=5,7$ at $\varepsilon=0.1$, including $0.2395495$ and $0.1418003$; the printed local $T_0$ residual $<-0.034$ | scalar estimates of Appendix B |
| `T0` | `existence/norm_t0.c` | $\|T_0\|\le0.998$ on the compact cells, through the primitive of the row norm | proposition on $T_0$ |
| `derivative_cusp` | `convexity/derivative_cusp.c` | $F_1$, $F_2$, $\alpha_1$, $\alpha_2$, $\alpha_f$ and $-xu_0''$ on $(0,0.05]$; the printed $C<0.035061$, $0.114$ and $-0.70$ | lemmas on $F_1,F_2$, on the coefficient bounds and on the curvature |
| `derivative_compact` | `convexity/derivative_compact.c` | the same quantities on $[0.05,\pi]$ by a recursive cover | same |
| `local1`, `local2` | `convexity/local_t1.c`, `convexity/local_t2.c` | the inequalities (ineqT1), (ineqT2) at $\varepsilon=0.002$, residuals $<-0.196$ and $<-0.188$ | propositions on $T_1$, $T_2$ |
| `T1`, `T2` | `convexity/operators.c` | $\|T_1\|,\|T_2\|\le0.998$ on $[0.002,\pi]$ by recursive covers, one LEAF line per cell | same |
| `closure` | `convexity/closure.c` | the radii of the fixed-point argument assembled from the exact targets | propositions on the corrections |
| `sign` | `uniqueness/sign.c` | the kernel-zero brackets: the 26 printed coefficients, the seven targets, 64 cusp cells and a recursively subdivided compact cover ending at $22/7>\pi$ | lemma on the kernel zero, Appendix C |
| `matrix` | `uniqueness/matrix.c` | outward binary64 enclosures of the kernel moments and of the seeds on the 720-node grid | Section 8 |
| `uniqueness` | `uniqueness/comparison.c` | the universal seeds ($A<0.855$), 240 coupled refinement steps, the cusp comparison and the 719 compact cells, $Q<0.982$ | Sections 8-9, Appendix C |
| `evaluators`, `kernel_tests`, `refinement_tests`, `operator_regression` | tests | regression checks of the evaluators, the kernel code and the exact refinement; not premises of the proof | |
| `reject_*` | the programs above | deliberate rejection tests with impossible targets; each must end UNRESOLVED | |

The last column of the summary table of Appendix B names, for every printed
bound, the stage that asserts it.

## Exact inputs

The approximation $u_0$ is the exact rational function of equation (u0): its
exponents and coefficients are the terminating decimals of
`inputs/approximation.json`, set verbatim in `core/basic_functions.c`
(`run.py` checks that the two agree). The targets of the fixed-point argument
are the exact values of `inputs/closure.json`: $A=2.7$, $\delta_0=10^{-7}$,
$q_0=q_1=q_2=0.998$, $r=6\cdot10^{-5}$, $\delta_1=5\cdot10^{-7}$,
$\delta_2=10^{-5}$, $\alpha_1=0.32$, $\alpha_2=1.383$, $\alpha_f=0.414$ and the
curvature budget $0.089$; `run.py` checks that `convexity/closure.c` declares the
same values. The kernel-sign polynomials and the seven targets of the lemma on
the kernel zero are the exact decimals of `inputs/kernel-sign.json`; `run.py`
checks that `uniqueness/kernel.c` declares the same coefficients and, when the paper
source is present next to the package, that the lemma in Appendix C states
the same polynomials and targets. The refinement grid is the exact binary64
file `inputs/grid720.bin`, whose last node is the double nearest to $\pi$.

## Constructions

**Evaluators (`core/`).** The Clausen functions are evaluated through the
series of Appendix B with the stated tail bound; every Fourier sum with the
weights $\sqrt{\tanh n}$ is truncated after $n=32$ and its tail is bounded by
the Fourier-tail lemma. Remainders of the cusp expansions are enclosed by the
explicit bounds of Appendix B, with the derivatives of remainders bounded from
the differentiated convergent series, never by differentiating an inequality.

**Kernel (`uniqueness/kernel.c`).** The regular part of $K$ is a series of Clausen
coefficients with an analytic tail; polynomial values on a cell use a
centered Taylor shift with interval coefficients (the shift need not be exact
for validity, because every coefficient is enclosed). The kernel-sign
certificate covers closed cells recursively, and both children of a split
must pass; there is no sampled root interpolation. The three-image second
derivative bound keeps 48 finite coefficients and the geometric tail stated in
Appendix C.

**Scratch arena (`core/ball.c`).** A single-thread arena for ball formulas:
every operation returns a fresh initialized ball, pointers remain stable, and
a reset releases only later temporaries. No ball is replaced by its midpoint
in a proof comparison; unsupported or nonfinite operations make the following
strict test fail; the array limits fail explicitly. Lattice conversion takes a
directed endpoint, scales by an exact power of two and rounds in the stated
direction with a fit check.

**Comparison (`uniqueness/comparison.c`, `uniqueness/refine.c`).** The port keeps the
monotone lower envelope at the cusp (without assuming $u/\sqrt x$ monotone),
the loss between continuous and lattice minimization, the strict
whole-interval exclusion test, checked signed 128-bit intermediates, the
regular remainders $M<0.083$ and $M_3<0.067$, the constructed-wave error
$6\cdot10^{-5}x$ in the comparison denominators, and the derivative cover
between grid nodes and over the entire punctured cusp cell. The same grid-index
selection is used when the moments are generated and when they are consumed;
it chooses a valid pivot node and no binary double comparison decides a
mathematical inequality.

## What the certificate does not do

The analytic arguments of the paper (the fixed-point scheme, the monotone
envelope and comparison arguments, the derivation of every inequality the
stages check) are premises, not outputs, of the computation. The regression
stages compare values with independent evaluations and are supplementary.
An aggregate PASS records one reproducible run; the fields `scope`,
`human_review` and `remaining` of `artifacts/certificate.json` state what is
left to human review.

## Validation and provenance

A release run executes every stage, the deliberate rejection tests and the 300
exhaustive discrete-minimization tests. The runner accepts only the conjunction
of the exact stage contracts and binds the local sources, inputs, executable
binaries, logs, arithmetic headers and libraries and the regenerated moment
file. A failure before or during compilation or numerical work invalidates the
attempt. The provenance of the sources is recorded in `import-manifest.json`.
