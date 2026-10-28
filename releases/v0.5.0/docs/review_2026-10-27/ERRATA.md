# Post-release mathematical errata identified in the v0.4.0 review

Review conducted 26–27 October 2026 (Asia/Shanghai). These corrections apply to
historical M0 material carried in the public v0.4.0 tree. They do **not** change
the canonical adaptive XO model or invalidate its registered v0.4.0 theorem.
The released files and their history have been left intact. Read this addendum
alongside the historical documents when reusing their claims.

**E1. The fixed-zone M0 criterion describes a continuous path vertex; it is
insufficient for a positive integer stopping time.** The affected passages are
[M0_GENERAL_Q, introduction and §A.5](../M0_GENERAL_Q.md), and repetitions in
the September checkpoints and research log. The public tagged version was also
checked directly at
[M0_GENERAL_Q §A.5](https://github.com/XuZhilin2007/sand-redistribution-study/blob/v0.4.0/docs/M0_GENERAL_Q.md).
The statements connecting the criterion to integer-time improvement require
the correction below. The formulas for the path and its vertex remain valid.

Use the historical M0 normalization

\[
 U_{\rm M0}(p)=\sum_i(p_i-1/K)^2,
\]

which differs by a positive factor K from the later XO loss. Let Q be the
initial/respray probability in the fixed selected zone, R=1−Q, and
\(S_n=\sum_{i\in n}q_i^2\), \(S_f=\sum_{i\in f}q_i^2\). Under the stated
assumptions \(0<Q<1\), \(0<\alpha\le1\), set

\[
 \rho=1-\alpha R\in(0,1),\qquad
 A=S_n+\frac{Q^2S_f}{R^2}>0,\qquad
 z_* = \frac{QS_f}{R^2S_n+Q^2S_f}>0.
\]

The exact path is \(z_t=\rho^t\), and its loss satisfies

\[
 U(z)=S_nz^2+S_f\left(\frac{1-Qz}{R}\right)^2-\frac1K
     =A(z-z_*)^2+C.
\]

Consequently \(S_n/Q>S_f/R\) is equivalent to \(z_*<1\). This is the
correct continuous-path condition. For the actual integer-time process the
precise condition for **some strictly better positive time** is instead

\[
 \exists t\in\mathbb Z_{\ge1}:U(\rho^t)<U(1)
 \quad\Longleftrightarrow\quad z_*<\frac{1+\rho}{2}.
\]

Indeed,

\[
 U(\rho)-U(1)=A(\rho-1)(\rho+1-2z_*).
\]

If this is negative, t=1 already improves on t=0. If it is nonnegative, then
\(\rho\le2z_*-1<z_*\). Every later \(\rho^t\) lies still farther to the left
of the vertex and has strictly larger loss than U(ρ). Equality gives a tie
between times 0 and 1, with all later times worse. Thus the equivalence covers
the whole infinite trajectory, rather than just a finite search.

When \(0<z_*<1\), the minimizing integer is obtained by comparing the two
times bracketing \(\log(z_*)/\log(\rho)\); comparing the losses is necessary.
Rounding that real-valued time to the nearest integer is not generally the
same operation. For \(z_*\ge1\), time 0 is optimal. The degenerate α=0 case
is a constant trajectory and is separate.

A legal, near-biased counterexample to the original integer-time reading is

\[
 K=6,\qquad n=\{1,2,3\},\qquad
 q=(17/100,17/100,17/100,49/300,49/300,49/300),\qquad \alpha=1/4.
\]

This can also be realized by a piecewise constant density with values 51/50
on [0,1/2] and 49/50 on (1/2,1]. It obeys every relevant fixed-zone M0
assumption. Exact arithmetic gives

\[
 Q=51/100,\quad c_n=17/100>49/300=c_f,\quad z_*=50/51,\quad\rho=351/400,
\]

\[
 U_0=1/15000,\qquad U_1=4405801/2400000000,\qquad
 U_1-U_0=1415267/800000000>0.
\]

Since \(\rho<z_*\), all subsequent losses exceed U₁. Time 0 is the unique
integer-time optimum despite the continuous interior vertex. With the **same
q** and α=1/100,

\[
 U_1-U_0=-14578333/500000000000<0.
\]

The continuous vertex is independent of α; existence of improvement at an
integer time can depend on α. The six previously tested distribution cases
do not establish the general integer-time claim. This correction does not
dispute their reported numerical outcomes.

The same density gives the same counterexample on every even K, including
M0's baseline resolution K=100: its loss is
\(U_K(z)=(51z/50-1)^2/K\). Changing K only multiplies all losses by a
positive constant, so this issue does not depend on using a small grid.

**E2. The exact iid histogram loss expectation has a distribution-dependent
variance term.** The affected equality is in
[M0_MECHANISM §2.7](../M0_MECHANISM.md), also verified in the
[public tagged file](https://github.com/XuZhilin2007/sand-redistribution-study/blob/v0.4.0/docs/M0_MECHANISM.md).
For N independent particles with common bin law p, write their empirical
histogram as \(\widehat p\). Unbiasedness and the multinomial variances give

\[
 \mathbb E\,U_{\rm M0}(\widehat p)
 = U_{\rm M0}(p)+\sum_i\operatorname{Var}(\widehat p_i)
 = U_{\rm M0}(p)+\frac{1-\sum_i p_i^2}{N}
 =\left(1-\frac1N\right)U_{\rm M0}(p)+\frac{1-1/K}{N}.
\]

The historical expression omits \(-U_{\rm M0}(p)/N\). Its noise term is
exact at uniformity and an upper bound elsewhere; it can be used as a nearby
noise-floor approximation when that approximation is declared.

An exact two-particle enumeration with K=2 and p=(3/4,1/4) gives

\[
 U_{\rm M0}(p)=1/8,\qquad
 \mathbb E\,U_{\rm M0}(\widehat p)=5/16,
\]

whereas the historical formula gives 3/8. The three outcomes with counts
0, 1, 2 have probabilities 1/16, 6/16, 9/16 and losses 1/2, 0, 1/2.
This is an exhaustive exact calculation, not Monte Carlo disagreement.
For the K=100, N=100000 canonical near-biased initial law, the omitted term
is about \(3.333\times10^{-8}\); the historical experimental noise estimates
can therefore remain numerically reasonable while the claimed equality is
false. Fixed-zone M0 evolves particles independently, so the iid calculation
applies at each fixed time under iid initialization. It should not be carried
over to an interacting adaptive particle system without checking that premise.

**E3. A displayed slope bound in the earlier robustness proof has the wrong
inequality direction, with an immediate repair.** In
[V03_B_ROBUSTNESS_PROOF §2(3)](../V03_B_ROBUSTNESS_PROOF.md), the argument uses
\(w\ge29/42\) and displays \(w/2-1<29/84-1<0\). The first inequality does
not follow from that lower bound. For example w=9/10 violates it. The proof
also has \(w<C<1\), which instead gives \(w/2-1<-1/2<0\), exactly the sign
needed. The other two slopes are positive by the displayed lower bound on w.
The theorem's conclusion does not need the erroneous numerical upper bound;
the monotonicity step is repaired by this one-line argument. This is a proof
presentation defect, distinct from E1's false integer-time implication and
E2's false exact equality.

**Scope and reproducibility.** E1 and E2 are checked from definitions
by [the new review auditor](../../tools/audit_v040_review.py), without importing
the original model or experiment modules. Its evidence includes a direct
first-step mass update for E1, the factorization establishing the corrected
criterion, and exhaustive sampling enumeration for E2. The current run is
[review-evidence.json](evidence/review-evidence.json). Neither correction is
used by the v0.4.0 proof chain; that chain has its own independently checked
prefix, entry, residence, W0 and tail arguments. This addendum is a correction
to the historical mathematical exposition, not a new release.
