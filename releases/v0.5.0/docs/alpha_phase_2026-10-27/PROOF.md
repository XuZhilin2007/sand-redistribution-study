# Sharp universal alpha interval: proof and exact certification

Post-release research, 27 October 2026 (Asia/Shanghai). The canonical model,
v0.4.0 scientific snapshot, prior reviews, and corrections are unchanged.
This is a computer-assisted mathematical proof, with explicit analytic
dependencies; it is not a claim of human peer review or proof-assistant checking.

Define

\[
 W_K=\{\alpha\in(0,1):\ \forall t\in\mathbb Z_{\ge0},\ U(p_t)\ge U(p_2)\},
 \qquad W=\bigcap_{K\ge6}W_K.
\]

All states, selectors, and losses below use the frozen definitions:

\[
 p_{0,i}=\frac{2K-2i+1}{K^2},\quad
 s_t=\min\arg\max_j\left(\sum_{i\le j}p_{t,i}-\frac jK\right),\quad
 U(p)=K\sum_i(p_i-1/K)^2.
\]

The respray weights are the bin integrals of the XO density with slopes
1/2, 3/2, 2, 0 and breaks 1/3, 2/3, 5/6. The update is
\(p'_i=(1-\alpha\mathbf1_{i\le s})p_i+\alpha F_sg_i\).
The activity threshold remains the existing one. Theorem A's invariant
\(D_m=b(1-b)\ge10/121\) guarantees activity forever, for all parameters here.

**1. The sharp result.** Let

\[
\begin{aligned}
 P_-(a)={}&42889525a^5-240679790a^4+555506236a^3\\
          &-613816360a^2+272381200a-37296000,\\
 P_+(a)={}&1164a^3-4361a^2+5068a-1136.
\end{aligned}
\]

Let \(a_-\) be the unique zero of \(P_-\) in \((6/25,1/4)\), and
\(a_+\) the unique zero of \(P_+\) in \((19/100,1/3)\). Then

\[
 \boxed{W=W_6\cap W_{10}=[a_-,a_+]},
 \qquad
 a_-=0.24615482985838558\ldots,\quad
 a_+=0.29164543617016553\ldots.
\]

The endpoints have algebraic degrees 5 and 3, respectively: the primitive
polynomials are irreducible modulo 29 and 11. For degree 5, excluding all
linear and monic quadratic factors suffices; for degree 3, excluding linear
factors suffices. The verifier performs these finite-field divisions exactly.
Root uniqueness is checked by Sturm counting and, independently, a strict
Bernstein derivative sign with opposite rational endpoint values. Exact
100-bisection rational brackets are in `exact_results.json`.

Two inclusions prove the result, without a connectedness assumption:

* Every parameter outside \([a_-,a_+]\) is excluded by K=10 or K=6.
* The entire closed interval belongs to every \(W_K\), with an exact finite
  complement and a uniform proof for K>=204.

Write

\[
 B=[123/500,7/24].
\]

A stronger intermediate fact is

\[
 B\subset W_K\quad\text{for every }K\ge7,\ K\ne10.
\]

Thus every other size admits a common interval extending strictly beyond
both universal boundaries. K=6 and K=10 determine the intersection exactly.

**2. The exact K=6 set and the global outer bound.** Define

\[
 R_6(a)=2328a^5-13378a^4+32735a^3-38636a^2+17183a-2080
\]

and let \(a_6\) be its unique zero in \((19/100,1/5)\). The exact result is

\[
 W_6=[a_6,a_+],\qquad a_6=0.19040773653717682\ldots.
\]

At time 1 the exact K=6 selector is 3 for \(a\le4/27\), 4 for
\(4/27<a\le1/3\), and 5 for \(a>1/3\). On
\([19/100,1/3]\), the word producing states 1,2,3 is (3,4,5), and

\[
 U_3-U_2=\frac{aR_6(a)}{2592},\qquad
 U_1-U_2=-\frac{aP_+(a)}{1296},\qquad U_0>U_2.
\]

On this entire interval, \(x_2=p_{2,5}>h\), and \(cq>h\), so the selector
is permanently m=5 from time 2. The remaining state lies on a geometric
line. More precisely, with \(\rho=1-a\), its loss is \(V(\rho^{t-2})\),
where \(V(z)=Az^2+Bz+C\) and \(A\ge0\).

The elementary integer-time criterion is exact: if
\(V(\rho)\ge V(1)\), then
\(B\le-A(1+\rho)\), so
\(V'(z)\le-A(1-\rho)\le0\) for \(0\le z\le\rho\).
Consequently every later integer state has loss at least that of time 3,
which is at least that of time 2. This uses the discrete comparison, in
accordance with the October 27 correction to the historical M0 criterion.

Sturm/derivative signs isolate the two loss crossings. Below 19/100, a
finite disjunctive proof gives \(U_3<U_2\). Its three left-open cells are

\[
 (0,4/27],\quad(4/27,913/5400],\quad(913/5400,19/100].
\]

On each cell the true time-1 selector is licensed. Every candidate time-2
selector is either excluded by another prefix or has a strictly negative
\(U_3-U_2\) polynomial. No selector continuity across a switch is used.
Between 19/100 and \(a_6\), the displayed polynomial is negative.
For \(a_+<a\le1/3\), \(U_1<U_2\). For the branch with time-1 selector 5,
\(U_1<U_2\) on \((1/3,2/3]\), and \(U_0<U_2\) on \([2/3,1)\).
These facts prove the entire K=6 description.

For K=10, the time-1 selector is 5 on \((0,4/45]\), 6 on
\((4/45,6/25]\), and 7 on \((6/25,1/3]\). A similar disjunctive proof
establishes \(U_3<U_2\) on

\[
 (0,4/45],\quad(4/45,37/225],\quad(37/225,6/25].
\]

On \((6/25,1/4]\) the first three selectors are (5,7,8), and

\[
 U_3-U_2=\frac{aP_-(a)}{45000000}.
\]

The polynomial has one simple zero, with negative sign below it. Hence
every \(a<a_-\) fails at K=10, time 3. This and the complete K=6 result
prove \(W\subset[a_-,a_+]\).

On the same K=10 interval, \(U_0,U_1>U_2\). A square-completion certificate
starting at time 4 gives all later comparisons, so

\[
 W_{10}\cap(6/25,1/4]=[a_-,1/4].
\]

In particular the lower endpoint is admitted, with exactly times 2 and 3
tied. At the upper endpoint, K=6 has exactly times 1 and 2 tied.

**3. The finite complement: all real parameters and all integer times.**
For K=7..203 except 10, the verifier covers B. K=10 is separately covered
on \([1/4,7/24]\). K=6 and the remaining K=10 interval use paragraph 2.

The exact affine upper envelope of \(D_{1,j}(a)\) gives rational closed
cells with a fixed time-1 selector. Every envelope junction is checked
against the actual minimum-index tie winner. Both limiting branches may
be proved on a closed cell, but the all-time argument is applied only to
the branch actually chosen at that parameter. No fictitious right branch
at a tie is substituted for the true trajectory.

For each cell the checker reconstructs both full updates. If time 1 was
not already a boundary selection, it strictly dominates every internal
time-2 prefix by a boundary prefix. The suffix is excluded by Theorem A.
Theorem A then gives permanent residence in \(\{\ell,m\}\).

It explicitly checks times 0 and 1. Remaining early times are handled by
a binary tree of boundary choices. The branch m is impossible on a cell
if \(x\le h\) throughout; the branch ell is impossible only if \(x>h\)
throughout. Each unpruned loss polynomial must be strictly above U2.
This over-approximation permits selector switches inside a cell without
requiring irrational switch locations to become rational cell endpoints.

An independent-of-word tail bound, described in paragraph 5, covers the
remaining times. The latest tail start actually needed is time 5.
There are 1,045 cells for the 196 ordinary sizes, plus one K=10 upper cell.
The producer checks signs by rational Bernstein subdivision. The witness
is then read and reconstructed using exact Sturm signs, with the recorded
tail start enforced and all endpoint tie winners checked again. Negative
controls reject a missing time-1 comparison, a true time-3 failure, a tail
started too early, and a missing K=32 smallest-index triple-tie branch.

**4. Uniform two-peak prefix enclosure for K>=204.** Put

\[
 h=1/K,\quad r=K\bmod6,\quad u_0=\tfrac12-(r\bmod2)h/2,\quad
 P(y)=2y-y^2,\quad \lambda=aP(u_0),\quad\rho=1-a.
\]

At grid points, the unscanned middle and upper pieces of D1 are

\[
 -y^2+(1+3\lambda/2)y-4\lambda/3,\qquad
 -y^2+(1+2\lambda)y-5\lambda/3.
\]

Their vertices are
\(y_M=1/2+3\lambda/4\), \(y_U=1/2+\lambda\).
On B, all scanned pieces are strictly increasing; a rounded middle peak
strictly beats the scanned endpoint, and a rounded upper peak strictly
beats the frozen suffix. Both nearest-grid vertices lie strictly inside
their proper kernel segments. Every global maximum is therefore a nearest
grid point of one of these two peaks, with lower-index tie selection.

Thus the actual time-1 selector satisfies one of the exact formulas

\[
 u_1=1/2+3\lambda/4+\eta h\quad(M),\qquad
 u_1=1/2+\lambda+\eta h\quad(U),\qquad |\eta|\le1/2.
\]

This is a two-branch enclosure, not a claim that both branches occur for
each parameter. The continuous peak difference is
\(\lambda(4-21\lambda)/48\), and grid loss is at most \(h^2/4\).
These bounds show M is the unique peak family for \(a\le1/4\) and U is
the unique peak family for \(a\ge13/50\). It suffices to cover

\[
 M:\ a\in[123/500,13/50],\qquad
 U:\ a\in[1/4,7/24].
\]

The overlap includes all possible changes of dominant peak. This avoids
asserting a single continuum trajectory at 16/63, where the existing
review correctly records different subsequential branches.

All constants are rational. With \(h_* =1/204\),
\(F_*=P(1/2-h_*/2)\), \(\lambda_-=123F_*/500\),
\(\lambda_+=7/32\), the checker verifies strict positivity of

\[
\begin{gathered}
 1/3-(4/3)(7/24)+\lambda_-/2,\quad
 (123/500)(3F_*/2-1),\\
 2/3-(1/2+3\lambda_+/4+h_*/2),\quad
 1/2+\lambda_--h_*/2-2/3,\quad
 5/6-(1/2+\lambda_++h_*/2),\\
 (3\lambda_-/4)^2-h_*^2/4,\quad
 (5/6-1/2-\lambda_+)^2-h_*^2/4,\\
 (3/16)(4-21(3/16))/48-h_*^2/4,\\
 [(13/50)F_*]\,[21(13/50)F_*-4]/48-h_*^2/4.
\end{gathered}
\]

The peak-gap function is monotone on the ranges used; the checker also
verifies \(\lambda_->2/21\). The strict central/suffix comparisons rule
out other peak families. Ties inside or between the two permitted families
remain covered by the closed eta enclosure.

**5. Exact moments, uniform entry, W0, and the infinite tail.** Set

\[
\begin{gathered}
 m=\lceil5K/6\rceil,\quad \ell=m-1,\quad
 b=5/6+(r/6)h,\quad L=b-h,\quad c=P(b),\quad q=(6-r)h/3,\\
 \nu=a[P(u_1)-\lambda(1-G(u_1))],\quad
 \mu=\lambda+\nu,\quad\beta=\rho\lambda+\nu.
\end{gathered}
\]

The true two-step density, whose bin integrals equal p2, is
\(\rho^2f_0+\beta k\) up to u0, \(\rho f_0+\beta k\) from u0 to u1,
and \(f_0+\mu k\) above u1, where \(f_0(y)=2-2y\).
For M, split at 1/3,u0,u1,2/3,L. For U, split at 1/3,u0,2/3,u1,L.
In particular

\[
 x_2=hX,\quad X=2-2b+h+\mu(6-r)/3,\quad
 \gamma=(c-hX)/(1-q),\quad e_i=p_{2,i}-\gamma g_i\ (i<m).
\]

Then \(\sum e_i=0\). Define
\(C_g=K\sum_{i<m}g_i^2\), \(Z=K\sum_{i<m}g_ie_i\),
\(E=K\sum_{i<m}e_i^2\).

For an affine density piece \(A f_0+Ck\) on [v,w], its average is
\(A(2-v-w)+Ck\); its square integral is
\((w-v)[\overline f^2+A^2(w-v)^2/3]\).
The exact cell covariance corrections are necessary, including their signs.
Let

\[
 t_1=(r\bmod3)/3,\quad t_2=(2r\bmod3)/3,\quad
 v_j=t_j(1-t_j),\quad \theta=(6-r)/3.
\]

Put \((A_2,C_2)=(1,\mu-\gamma)\) for M and
\((\rho,\beta-\gamma)\) for U. With
\(I_Z=\int_0^L k(f_2-\gamma k)\), \(I_E=\int_0^L(f_2-\gamma k)^2\),

\[
\begin{aligned}
 C_g={}&3/2-2\theta h-h(v_1+v_2/4),\\
 Z={}&I_Z+h^2(v_1\rho^2+v_2A_2/2)
       -h[v_1(\beta-\gamma)+v_2C_2/4],\\
 E={}&I_E-\frac{h^2}{3}[\rho^4u_0+\rho^2(u_1-u_0)+L-u_1]\\
      &-h[v_1(\beta-\gamma)^2+v_2C_2^2/4]
       +2h^2[v_1\rho^2(\beta-\gamma)+v_2A_2C_2/2].
\end{aligned}
\]

The code independently sums seven runs of bin-average arithmetic
progressions and compares coefficients in \(\mathbb Q[h,a,\eta,\gamma]\).
Both branch formulas and all residues agree identically. Including the
direct p2 prefix square norm and the complete p1 square norm gives 60
formal identities. Controls detect a reversed covariance and a missing
mixed term. Ninety-five full-vector cases bind the formulas to the actual
updates; 40 further replays use integer masses and kernel intersection
lengths. These replays are checks, not the source of the infinite quantifier.

W0 is checked by direct density norms with exact cell variances. For p1,
also let \(t_3=(5r\bmod6)/6\), \(v_3=t_3(1-t_3)\). Then

\[
\begin{aligned}
 U_1={}&\int_0^1 f_1^2-1
 -\frac{h^2}{3}[\rho^2u_0+1-u_0]
 -h\lambda^2(v_1+v_2/4+4v_3)\\
 &+2h^2\lambda(v_1\rho+v_2/2-2v_3).
\end{aligned}
\]

For U2, the same prefix moment formula is used with gamma=0, plus
\(hX^2\) and the frozen suffix norm
\(4(1-b)^3/3-(1-b)h^2/3\), then subtract 1.
Also \(U_0=(1-h^2)/3\).

For entry, the three lower-region D2 quadratics are

\[
\begin{aligned}
 -\rho^2y^2+(2\rho^2+\beta/2-1)y,\\
 -\rho^2y^2+(2\rho^2+3\beta/2-1)y-\beta/3,\\
 -\rho y^2+(2\rho+3\beta/2-1)y-4\rho\lambda/3-\nu/3.
\end{aligned}
\]

Coefficientwise upper bounds and \(\max_y(-Ay^2+By+C)=C+B^2/(4A)\)
bound them. Any unswept middle interval increases to a value no greater
than \(2/9-\mu/3\). On the upper region, the final derivative is bounded
below by \(2\mu-2/3\). For U, the intermediate swept upper segment also
has positive derivative, bounded by
\(2\rho(1-u_1)+2\beta-1\). The boundary value is at least
\(5/36-(2\mu_{\max}-1/3)h_{\max}\). The interval proof establishes a
strict gap above all lower-region coordinates; upper-region coordinates
are strictly increasing to ell. Thus all global maxima of D2 are boundary
indices. Theorem A supplies permanent residence thereafter.

After entry, for \(z=\rho^{t-2}\), both boundary updates preserve

\[
 p_{t,i}=\frac{c-x_t}{1-q}g_i+ze_i\quad(i<m).
\]

Both branches multiply the residual by rho; no fixed future word is needed.
The frozen suffix is unchanged. Expansion gives

\[
 U_t-U_2=A\left(x_t+\frac{B+2dz}{2A}\right)^2+H(z),\quad
 H_2=E-d^2/A\ge0
\]

by the Gram inequality on the true vectors. With

\[
 \widehat A=hA=1+hC_g/(1-q)^2,\quad
 B=-2cC_g/(1-q)^2,\quad d=-Z/(1-q),\quad f=cZ/(1-q),
\]

the exact tail expressions are

\[
\begin{aligned}
 H'(\rho)&=2(E-d^2h/\widehat A)\rho+2f-Bdh/\widehat A,\\
 H(\rho)&=E(\rho^2-1)+2f(\rho-1)-hX(\widehat A X+B+2d)
           -h(B+2d\rho)^2/(4\widehat A).
\end{aligned}
\]

These extend regularly to h=0. Positive interval denominators are checked.
Since H2>=0, a negative H'(rho) makes H decreasing on [0,rho]. A positive
H(rho) therefore covers every integer time t>=3. For the finite complement
the same argument uses \(w=\rho^{T-2}\), covering every t>=T.

**6. Continuous boxes and complete assembly.** Each branch's alpha interval
in paragraph 4 is initially divided into 32 closed rational intervals.
For each of six residues, divide h in [0,1/204] into four closed intervals
and eta into [-1/2,0] and [0,1/2]. This gives 3,072 initial boxes.

Natural exact Fraction interval arithmetic evaluates all formulas above.
An unresolved box is bisected along the largest of the alpha width, h
width, and eta width multiplied by h_max. Both children are checked; no
sampling or pruning removes an unproved region. Depth exhaustion is a
failure, never acceptance. The final run visits 20,554 boxes and has
11,813 terminal boxes, with no unresolved boxes. The binary covering
identity is 20,554 = 2(11,813) - 3,072.

Every terminal box proves

\[
 U_0-U_2>3/20,\qquad U_1-U_2>1/1000,\qquad
 H'(\rho)<-1/10,\qquad H(\rho)>1/200,
\]

as well as a boundary-entry margin greater than 3/100 and positive required
derivatives. Box variables include nonphysical rounding combinations; the
enclosures are sufficient bounds, and Gram is invoked only for true states.

These boxes cover **every K>=204 and every real alpha in B**. The finite
complement meets them at 203/204. Combining this with paragraph 2 gives
\([a_-,a_+]\subset W\), completing paragraph 1.

The sole global minimum is time 2 for every K throughout the open sharp
interval. At its lower endpoint K=10 has times 2 and 3 tied; at its upper
endpoint K=6 has times 1 and 2 tied. The strict certificates exclude other
equalities there. No late-time or K-to-infinity limit shrinks these bounds.

**7. Fixed-size discontinuities and a fourth-state obstruction.** Exact
certificates also give

\[
 [21/25,43/50]\subset W_{10},\qquad (13/40,1/3]\subset W_{13}.
\]

Since 1/4 is in W10 and 1/2 is not, W10 is disconnected. At K=13,
alpha=13/40 gives selectors (6,9,11,...) under the true tie rule and

\[
 U_1-U_2=-31191/22497280<0.
\]

Immediately to its right s1=10; the right branch is certified all-time
on the stated half-open interval. Its time-2 loss jumps by

\[
 U_2(13/40+)-U_2(13/40)=-155929/16224000.
\]

The left branch is a fixed polynomial near this tie, with strictly
negative U1-U2 at the tie, so an excluded left neighborhood exists.
Also alpha=81/250 is an exact failure between 1/4 and the right interval.
Thus W13 is disconnected and is not closed relative to (0,1).

At K=10, alpha=9/10 has word (5,8,9,8,...) and

\[
 U_0>U_2,\quad U_1>U_2,\quad U_3>U_2,\qquad
 U_4-U_2=-4697222139/20000000000000<0.
\]

It disproves a general reduction of fixed-K membership to just times
0,1,2,3. This is a finite transient obstruction; it is not evidence of an
asymptotic boundary for the universal interval.

**Reproduction and trust boundary.** Run

```text
python -B experiments/alpha_phase_2026_10_27/verify_all.py
```

The mathematical runtime uses the standard library only. It reconstructs
the proof, checks target coverage, root ordering and source/report hashes,
reruns the published baseline verification, and checks that its 37 frozen
scientific files and the prior review/correction files remain unchanged.
The exact interval, polynomial and source records are kept in the new
experiment directory. The optional plot uses an existing matplotlib
environment; it is not part of the proof.

Common analytic dependencies are Theorem A's invariants and residence,
elementary quadratic completion and Gram, the validity of exact rational
interval operations, and the stated grid alignment behind the formal sums.
The two sign methods and the bin-sum audit reduce implementation risk;
they do not constitute a formal proof assistant or independent human review.
No existing release, tag, certificate, review record, or Git history is changed.
