# Sand v0.4.0: literature review, mathematical audit and research assessment

Prepared 27 October 2026, following work on 26–27 October (Asia/Shanghai).
This is a post-release assessment of the published theorem. It is neither a
new release nor a claim of human peer review or proof-assistant verification.
The main supporting files are the [literature ledger](LITERATURE.md),
[historical errata](ERRATA.md), and [execution evidence](evidence/).

**The main finding is that the registered v0.4.0 theorem survives this review.**
I found no substantive defect in its analytic assembly or exact certification
for the frozen canonical XO model. Both the research snapshot and the public
verification package passed their six checks, and a new independent algebraic
check supports the most delicate tail identities. Two genuine mathematical
corrections are needed in earlier M0 exposition, together with one readily
repairable inequality in the v0.3 robustness proof. None is a dependency of
the v0.4.0 main theorem. Correctness evidence is substantially stronger than
finite simulation; originality and significance remain separate questions.

The most defensible mathematical contribution is a uniform comparison of all
integer times for a growing family of discrete systems, with a fixed parameter
interval. It is a narrow theorem about a prescribed trajectory. The rotation
theory behind its periodic dynamics is largely established, and the certificate
tools are established. The next step should be focused external mathematical
review of a compact presentation, with further work selected for structural
insight rather than another version or a larger table of certified cases.

**The baseline was recovered from files and remote records.** The public
[v0.4.0 release](https://github.com/XuZhilin2007/sand-redistribution-study/releases/tag/v0.4.0)
was published on 26 October 2026 at 13:11:28 China time. The public annotated
tag object is `14144ecaee268068bf14c740f4173f1fcfa56048`, pointing to commit
`2d8820db3fd65733c9137c955fb1acc0b5466cc6`. These were checked through GitHub
read APIs. The private scientific freeze is
`9e2f0fc46787e4206f5729a9976d81c5ff235dc7`; the current research checkout was
`2c337bb`. The private repository's own v0.4.0 tag resolves to `e05f4c1`, a
later publication-record commit, so those tag names must not be used
interchangeably across the two repositories.

The public verification ZIP already present in the workspace has SHA-256
`20328551989ac1fe7ff7bd5d5b6c28a9eac52bf5d7b8dcbd11b3ac92866481f9`.
That hash matches the release asset digest freshly read from GitHub. Its
extracted tree was reverified here. This review did not claim a new anonymous
download; the earlier download provenance remains the earlier release record.
Public formatting adaptations, optional Git provenance and preserved original
sources are documented in
[V040_PUBLICATION](https://github.com/XuZhilin2007/sand-redistribution-study/blob/v0.4.0/docs/V040_PUBLICATION.md).
The local frozen README still describes a local freeze; that is a navigation
issue resolved by the dated publication record, not evidence that v0.4.0 was
unpublished. The frozen README was preserved.

The relevant history is: the September S3 boundary reduction; v0.2.0's
extension to 0<α<1 and its finite-window distinction; v0.3.0's α=1/4
all-K minimum, canonical word tables, robustness results and explicit
literature dictionary; then v0.4.0's real parameter interval. The v0.4 research
log also records repaired earlier sign tests, omitted times 0 and 1, and an
invalid continuity argument at ties. Those superseded PASS records are not
proofs of the final claim. The final independent checker reconstructs the
mathematics from the witness words and definitions; it does not simply accept
their status fields.

**The exact theorem has the following meaning.** For integer K≥6, use

\[
 p_{0,i}=(2K-2i+1)/K^2,\quad F_j=\sum_{i\le j}p_i,\quad D_j=F_j-j/K,
 \quad s=\min\arg\max_{1\le j\le K}D_j.
\]

The respray weights are \(g_i=G(i/K)-G((i-1)/K)\), where G is linear with
successive slopes 1/2, 3/2, 2, 0 and breaks 1/3, 2/3, 5/6. On an active
step, \(p'_i=(1-\alpha\mathbf1_{i\le s})p_i+\alpha F_sg_i\). The activity
test is \(\max D_j>10^{-12}\); the canonical invariant proves it at every
time. There is thus a definite infinite trajectory with this tie rule.

For every real \(\alpha\in J=[1999/8000,2001/8000]\),

\[
 \forall K\in\mathbb Z_{\ge6}\ \forall t\in\mathbb Z_{\ge0},\qquad
 U(p_t)\ge U(p_2),\qquad U(p)=K\sum_i(p_i-1/K)^2.
\]

This compares stopping times on an already fixed feedback orbit. It does not
optimize the selected prefix, the respray rule or the removal fraction. It
does not assert monotone loss after time 2, uniformization of the state, or
an optimum under other metrics. Exact interval arithmetic covers irrational
as well as rational α; it does not require representing every real α by a
finite rational trajectory. No physical calibration follows from this model.

**The proof has a genuine route through its infinite quantifiers.** The table
states which mathematical arguments, rather than which numerical samples,
cover each region.

| Obligation | Evidence reviewed | Assessment |
|---|---|---|
| No stopping; frozen suffix; fixed mass \(F_m=c\) | Theorem A §2: \(G_j=1\) for j≥m, strictly dominated canonical suffix; \(D_m=b(1-b)\ge10/121\) | Direct invariant induction; no infinite simulation assumption |
| Permanent boundary residence | Theorem A §§3–4 and S3: internal-coordinate barriers, endpoint comparison, separate K=6 argument | The smallest-index convention excludes an interior tie at first boundary selection; the exception is handled |
| True first two selectors for large K | Prefix lemma for K≥33; global comparison of both continuous peaks and grid rounding loss | Global dominance is proved, not inferred from a local quadratic maximum |
| Entry by state time 2 | D₂ short proof: all internal coordinates below 3/25, \(D_{2,\ell}\ge247/1980\); small-K witnesses | Every maximizer is a boundary index, with a strict internal margin |
| Times 0 and 1 | W0: density comparison, fixed-cutoff L1 contraction, cell variance, Bernstein bounds; K≥145 plus finite complement | The adaptive argmax map is not assumed contractive; cutoff differences are treated separately |
| All t≥3 for K≥204 | Exact two-branch residual identity, quadratic completion, Gram sign, 48 rational boxes | Boxes cover six residues, all h=1/K≤1/204, all true rounding errors and all α∈J |
| All times for 6≤K≤203 | Independent Bernstein checker reconstructs prefix, branch permissions, finite losses and an infinite tail | 198 needed K, plus five redundant large anchors; 203 total target K, 207 pieces and 80,246 sign obligations |

The two K ranges meet at 203/204. W0's analytic cutoff at 145 is compatible
with that split; it does not leave K=145..203 or K=204 untreated. The small-K
certificates cover all real α in the closed common interval, not just its
endpoints and midpoint. At K=16,32,96,160, the finite checker uses a closed
left piece and an open-left right piece with each branch's own U₂. K=32's
three-way first-step tie is included. For larger K, the prefix rounding
parameter covers its true ties, including K≡32 mod 64 at α=1/4. The proof
does not assume continuity of the state across a selector change.

After residence, with \(x=p_m\), \(\rho=1-\alpha\) and \(z=\rho^{t-2}\),

\[
 p_{t,i}=\frac{c-x_t}{1-q}g_i+z e_i\quad(i<m),\qquad \sum_{i<m}e_i=0.
\]

Both boundary branches multiply the same residual e by ρ. Expanding the
loss gives

\[
 U_t=Ax^2+Bx+C+2z(dx+f)+Ez^2,
 \qquad U_t-U_2=A\left(x+\frac{B+2dz}{2A}\right)^2+H(z).
\]

The factor 2 in the mixed term is essential and is present. With the scaled
moments \(C_g=K\sum g_i^2\), \(Z=K\sum g_ie_i\), \(E=K\sum e_i^2\),

\[
 A=K+C_g/(1-q)^2>0,\quad d=-Z/(1-q),\quad
 H_2=E-d^2/A\ge0
\]

by the Gram inequality. The exact boxes establish
\(H'(\rho)<-3/25\) and \(H(\rho)>1/100\). Hence H is decreasing on
[0,ρ], and every integer t≥3 has \(U_t-U_2>1/100\) for K≥204.
Future periodic words are unnecessary for this argument. For small K the
finite checker licenses an initial word and then proves a square-completion
bound for all sufficiently small z; this is the infinite tail missing from
a mere finite trajectory search.

The nonphysical points introduced by treating h and the rounding error as
independent box variables are harmless over-approximations of the formulas.
Gram positivity is needed on actual states, where it is proved from vectors;
it is not being asserted for arbitrary fictitious residuals in a box.

**The verification has several distinct levels of independence.** Rerunning
the release tools establishes reproducibility, not a second derivation. The
finite auditor is stronger than a replay of a generator: it rebuilds global
argmax inequalities, both formal vector updates, loss coefficients and tail
conditions, using a different Bernstein sign method from the generation
kernel. Its use of common definitions remains a legitimate common premise.
Theorem A's residence argument and the meaning of the polynomial/interval
algorithms remain analytic and software trust obligations. No Python PASS
turns them into proof-assistant theorems or human validation.

The new [audit_v040_review.py](../../tools/audit_v040_review.py) addresses a
particular remaining weakness: the release checks its tail moment identities
against finitely many full vectors. This new check independently sums bin
averages in seven arithmetic-progression runs, including the two mixed
kernel cells, and compares the result coefficient by coefficient with direct
affine-density integration and cell covariance corrections. All four symbols
h, α, η and γ remain formal. There are 18 zero polynomial differences across
six residues. Controls detect a reversed covariance sign and a missing mixed
term. This supports the universal algebraic identities, rather than merely
checking more points. It still relies on the documented grid alignment and
does not replace the positivity enclosures.

| Actual run | Result | What it establishes |
|---|---|---|
| Research snapshot `verify_v040.py` | PASS, 6/6 checks, 37 frozen files unchanged, 49.675 s | Reproduction of the frozen scientific proof checks |
| Public ZIP tree `verify_v040.py` | PASS, 6/6 checks, 78 frozen files unchanged, 58.114 s | Public package reproduces the same theorem checks |
| New formal moment comparison | PASS, 18 identities in \(\mathbb Q[h,\alpha,\eta,\gamma]\) | Independent bin-sum derivation agrees identically with the tail formulas |
| New integer-mass replay | PASS, 93 cases through t=60 | A different arithmetic representation checks actual updates, selectors and moments, including ±10⁻¹² near six ties; finite evidence only |
| New conjugacy checks | PASS, 21 parameter cases, including exact x=h | Checks the algebraic dictionary and branch convention; general translation also derived above |
| Scope and historical counterexamples | All reproduced exactly | Restricts possible generalizations and establishes the errata |

The evidence JSON contains full commands, actual exit codes, source hashes,
parameter lists and limitations. The final new auditor takes about one second
on this environment. Its finite tests are not presented as a proof for all K
or all t. The historical 610-assertion v0.3 suite and all its large certificate
generation jobs were not rerun; this review does not claim a complete new
audit of every supplementary v0.3 theorem. The B/B-II and C proof scopes were
inspected, and T7/G were checked at their reduction and literature interfaces.
The α=1/4 minimum itself is included in the fully rechecked v0.4 interval.

**Theorem-to-literature comparison gives a restrained originality judgment.**
The source identifiers below resolve to the titles, publication records,
exact premises and access levels in [LITERATURE.md](LITERATURE.md).

| Sand result or method | Established theory / exact limitation | Current classification |
|---|---|---|
| U as nonuniformity; fixed-selector L1 contraction | U is Pearson chi-square divergence; stochastic-operator contraction/data processing [L11] | Known metric and elementary specialization |
| Support invariant and finite boundary entry | Kernel geometry and elementary barriers; interval rotation theorems start after the needed reduction | Specific derivation for this model; priority not established |
| CYCLE rotation rationality for algebraic α | Exact conjugacy to the unequal-slope family [L1, Theorem 2; L2, Corollary 5.3] | Direct consequence of existing theory after the checked dictionary |
| Actual attracting periodic cycle | Rational rotation needs the appropriate endpoint case [L1, Theorem 3]; generic results [L5–L8] do not prove every canonical parameter | Inherited qualitative theory with certificate obligations remaining |
| Canonical explicit words, K thresholds and endpoint exclusions | Rotation/symbolic mechanisms are established [L1–L4,L9]; those theorems do not list Sand's integer lattice | Technically distinct explicit classification; broader significance uncertain |
| Residual decomposition and quadratic tail certificate | Linear invariant-subspace algebra, Gram inequality and completing a square; safety and quadratic-certificate precedents [L12–L14] | Standard tools assembled into a useful proof for this model; the method alone is not novel |
| Two-update minimum, uniformly in K and α∈J | No implication established from the reviewed asymptotic, averaging, sandpile or ergodic-average results [L1–L17] | Potential independent contribution; novelty remains unresolved |
| Fixed-zone M0 internal minimum | Quadratic minimization along a geometric sequence, with the discrete correction in E1 | Elementary structural result with a corrected quantifier |

In particular, each fixed-selector matrix is stochastic, but it generally
does not preserve the uniform target. Data processing transforms both
probability laws; it cannot be used to make U monotone while holding the
target fixed. Also the feedback map itself is discontinuous at ties. At an
adjacent global tie \(p_{s+1}=h\),

\[
 T_{s+1}(p)-T_s(p)=\alpha h(g-e_{s+1}),
\]

where \(e_{s+1}\) is a coordinate unit vector. Fixed-selector contraction
therefore does not prove contraction or continuity of adaptive feedback.
This matters for perturbations of the canonical initial state, which already
has an initial tie for every odd K. It also explains why a naive continuity
argument cannot replace the existing tie enclosures.

The full system is piecewise linear in p, or piecewise affine in simplex
coordinates, with conserved suffix directions. It is not automatically a
strict contraction on the full simplex. On the scalar CYCLE interval the two
slopes are strictly contracting for each fixed K, but \(1-\alpha q\to1\) as
K grows. Generic contraction or stability theorems do not supply uniform
constants for that degenerating family. These observations identify why the
all-K comparison requires additional work, without asserting that such work
has no antecedent elsewhere.

**The consequential corrections concern historical scope, not a failure of
the current main proof.** [ERRATA.md](ERRATA.md) gives derivations and exact
counterexamples. In fixed-zone M0, \(c_n>c_f\) guarantees a continuous path
vertex, whereas strict improvement at a positive integer time is equivalent
to \(z_*<(1+\rho)/2\). The legal K=6 distribution stated there has
\(c_n>c_f\), yet α=1/4 gives time 0 as the unique global integer optimum.
The same distribution improves at time 1 when α=1/100. Thus the historical
α-independent existence claim cannot be retained for integer operations.

The exact iid sampling formula is
\(\mathbb E U_{\rm M0}(\widehat p)=U_{\rm M0}(p)+(1-\sum p_i^2)/N\),
not the constant uniform noise term added to U. The old numerical noise
estimates may remain close, but their displayed exact equality is false.
Finally, B's proof displays a slope bound in the wrong direction; its needed
negative sign follows immediately from w<C<1. This is a repairable line in
a supporting historical proof, with a different severity from the M0 claim.
The addendum preserves all original files and distinguishes these cases.

The review also exactly reconfirmed the published limitations. Outside J,
K=9, α=1/5 gives \(U_3-U_2=-83303117/5978711250\); K=6, α=7/24 gives
\(U_1-U_2=-161/11943936\). On the canonical K=9, α=1/4 orbit,
\(TV_3-TV_2=-337/78732\), with \(TV=\frac12\sum|p_i-1/K|\).
For the existing P3 perturbation that transfers 1/10 from bin 5 to bin 7 at
K=9, preserving F_m and the suffix,
\(U_3-U_2=-210896134981/9917964518400\). These finite exact inequalities
disprove the corresponding unrestricted generalizations. They do not disprove
smaller, carefully quantified neighborhoods or suitable classes of metrics.

**The significance is currently that of a rigorous case study.** Uniformity
in a growing dimension and in a real interval makes the theorem mathematically
more informative than a single-trajectory observation. Its assumptions remain
highly selected: one initial profile, one support cutoff and piecewise density,
one tie convention and one quadratic loss. The interval is certified, not
sharp, and the proof constants are sufficient conditions. There is no evidence
here of a general principle about real sand, a new rotation theory, or a broad
redistribution universality class. A specialist may judge the current result
suitable for a short mathematical note or an illustrative computer-assisted
case study; the evidence does not justify predicting publication priority or
venue acceptance.

**Three research directions deserve priority, ordered by the knowledge they
could add.** They are proposed questions, not newly initiated theorem projects.

1. **Identify a structural class with a uniform early minimum.** Ask for
   checkable conditions on the initial density and supported respray kernel
   that ensure entry into a low-dimensional invariant family and attainment
   of the minimum within a fixed number of steps independent of K. Connect
   the residual representation to switched positive operators and safety
   certificates [L11–L14]. A meaningful advance would include a nontrivial
   family beyond XO, a verifiable finite-time criterion and examples where
   that criterion's hypotheses fail. Merely restating square completion as a
   general theorem would add little. The difficulty is proving residence and
   uniform prefix/loss margins despite discontinuous argmax selection and a
   scalar slope tending to one. Pursue it if existing inequalities can be
   expressed in a few stable geometric/moment conditions and verified for a
   second substantive family. Reassess it if those conditions only encode
   the original XO proof or small perturbations require an ever-growing
   prefix length. This has the greatest potential significance and the
   greatest risk of over-generalization.

2. **Characterize the component of valid removal fractions containing 1/4.**
   Define \(\mathcal W=\{\alpha\in(0,1):\forall K\ge6,\forall t\ge0,
   U(p_t)\ge U(p_2)\}\), with the frozen model. Determine the component
   containing J, its endpoint inclusion, and its binding K/time comparisons.
   The existing values near 0.2461548 (K=10) and 0.2916454 (K=6) are finite
   search candidates, not proved all-K endpoints. A meaningful advance would
   prove which finite constraints bind, rule out missed time minima and
   reentrant parameter components, and certify the remaining infinite K
   uniformly. The major difficulty is selector bifurcation: at α=16/63,
   two first-step peaks compete and the documented U₂ continuum limit has
   different subsequences. A single continuous branch cannot cross that
   point. Rotation plateau boundaries [L1–L2] concern asymptotic dynamics;
   they need not be loss-minimum boundaries. Pursue this after independently
   isolating the candidate loss roots and checking that no additional branch
   binds; abandon a simple-interval conjecture if an exact reentry or new
   binding comparison appears. This is concrete and feasible, with moderate
   value if it explains the boundary mechanisms rather than just enlarging J.

3. **Test survival under perturbations that break the exact support invariant.**
   A useful question is whether a quantified early or near-optimal stopping
   result persists for \(g^\varepsilon=(1-\varepsilon)g^{XO}+\varepsilon u\),
   or for a specified neighborhood of initial densities, with the metric and
   observation scale declared. This connects operator perturbation and
   contraction stability [L7–L8,L11], but late periodic-orbit stability is not
   early-minimum stability. Full-support leakage destroys the frozen suffix;
   the ε→0 and infinite-time limits may not commute. The existing P3 and TV
   counterexamples already exclude broad claims. A meaningful advance is an
   explicit regime with a uniform bound on stopping regret or a precise
   discontinuity/fragility result. Specify the norm and its K scaling before
   testing an initial-state radius, and handle every accessible tie branch.
   Pursue the direction if direct perturbation bounds control both the prefix
   and the infinite tail; switch to a fragility result if exact arbitrarily
   small perturbations destroy the property. Another fixed finite scan alone
   would not settle it.

**The recommended immediate decision is external review, supported by the
existing compact note and this addendum.** Send the mathematical note only
after its author checks that the new historical corrections are acknowledged
where those historical claims are reused. A reviewer in interval dynamics
should be asked to assess the reduction, endpoint conventions and originality;
a reviewer familiar with rigorous numerics should assess the analytic
dependencies and exact certification. Their central questions should be
whether the residence induction and moment assembly are sound, whether the
uniform transient comparison is already implied by a known theorem, and
whether the narrow result has sufficient interest as a note. This review
does not contact anyone or make claims about their availability.

Preparing a concise paper draft is reasonable as a means of making those
questions reviewable; it should lead with the uniform transient theorem,
credit the existing rotation theory, separate the certificate appendix, and
describe the result as a model theorem. A broad v0.5 phase or further theorem
counting should wait for that assessment. If a specialist identifies a direct
general antecedent, reposition Sand as an explicit application with useful
certificates. If the theorem is valid but too tailored, direction 1 provides
a more substantive reason to continue than another parameter scan.

**Changes made by this review are additive and reproducible.** Added this
assessment, the primary-source literature ledger, the errata addendum, the
independent review auditor, and new evidence records. No released scientific
file, historical report, numerical witness, model rule, Git tag or Git history
was edited. No commit, publication or external message was made. The
pre-existing untracked ZCODE documents and scratch material were retained.

To reproduce the additional audit from the repository root:

```powershell
.venv/Scripts/python.exe -B tools/audit_v040_review.py --report scratch/new-v040-review.json
```

Choose a new report name; the tool refuses to overwrite an existing evidence
record. The original `python -B tools/verify_v040.py` command still verifies
the frozen baseline independently of these new files. Remaining limitations
are substantive: selective literature coverage, no external human review of
this audit, no proof-assistant checking, no exhaustive audit of every older
supplement, no sharp parameter window, and no evidence establishing broad
novelty or physical applicability. The current evidence supports the main
theorem within its stated model and supports seeking a specialist's judgment
about its mathematical value.
