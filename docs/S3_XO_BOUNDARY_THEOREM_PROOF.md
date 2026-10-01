# S3 — Canonical XO Boundary-Reduction Theorem：证明本体

> Deposit date: 2026-09-30. **S3: CLOSED — PASS**.
>
> **[New Mathematical Result] — proof deposited.** This document deposits Astra's completed P1–P4 proof from the preceding S3 proof-attempt response. It does not initiate new research or extend the audited statement. [Registration](S3_XO_BOUNDARY_THEOREM_REGISTRATION.md) records T1–T5; [audit provenance](S3_XO_BOUNDARY_THEOREM_AUDIT.md) records the independent GPT-6.1 Sol verdict supplied by the Owner: **PASS — no substantive gap found**.

## 1. Model, notation, and exact scope

Fix an integer $`K\ge6`$. The state is a nonnegative mass vector $`p_t`$ with total mass one. Set

$$
F_{t,j}=\sum_{i\le j}p_{t,i},\qquad
D_{t,j}=F_{t,j}-j/K,\qquad d_t=\max_jD_{t,j}.
$$

Selection is $`s_t=\min\operatorname*{argmax}_jD_{t,j}`$. The stopping rule is $`d_t\le10^{-12}`$. On an active step, remove one quarter of each selected-prefix bin and redistribute the removed mass $`M_t=F_{t,s_t}/4`$:

$$
p_{t+1,i}=\begin{cases}
\frac34p_{t,i}+M_tg_i,&i\le s_t,\\
p_{t,i}+M_tg_i,&i>s_t,
\end{cases}
\qquad g_i=G(i/K)-G((i-1)/K).
$$

The frozen XO CDF is

$$
G(y)=\begin{cases}
y/2,&0\le y\le1/3,\\
3y/2-1/3,&1/3\le y\le2/3,\\
2y-2/3,&2/3\le y\le5/6,\\
1,&5/6\le y\le1.
\end{cases}
$$

Canonical initialization is the grid difference of $`2y-y^2`$:

$$
p_{0,i}=\frac{2K-2i+1}{K^2},\qquad D_{0,j}=(j/K)(1-j/K).
$$

All statements concern exact arithmetic, this initialization, this XO law, removal fraction $`\alpha=1/4`$, and smallest-index selection. They make no claim for other laws, initial states, removal fractions, finite-precision implementations, or 2D models.

| Symbol | Definition |
|---|---|
| $`h`$ | $`1/K`$, the uniform bin share |
| $`r`$ | $`K\bmod6`$, in $`\{0,1,2,3,4,5\}`$ |
| $`m`$ | $`\lceil5K/6\rceil`$ |
| $`\ell`$ | $`m-1`$, the left boundary index |
| $`b`$ | $`m/K`$ |
| $`c`$ | $`2b-b^2`$, the conserved mass $`F_{t,m}`$ |
| $`a`$ | $`b(1-b)=c-b`$, the conserved discrepancy $`D_{t,m}`$ |
| $`q`$ | $`g_m=(6-r)/(3K)`$, the redistribution share of bin $`m`$ |
| $`x_t`$ | $`p_{t,m}`$, the boundary coordinate |
| $`t_*`$ | $`\min\{t:s_t\in\{\ell,m\}\}`$, whose existence is proved below |

For the closed form of $`q`$, write $`K=6n+r`$, so $`m=5n+r`$ and $`b=5/6+rh/6`$. Then $`2/3\le(m-1)h<5/6`$, and

$$
q=1-G((m-1)h)=1-[2(b-h)-2/3]=(6-r)h/3.
$$

## 2. Existing invariant dependency — T1 [Existing Math]

The previously established XO support-endpoint invariant supplies, for every time,

$$
s_t\le m,\quad F_{t,m}=c,\quad D_{t,m}=a,
\quad D_{t,j}<a\ (j>m),\quad a\ge10/121>10^{-12}.
\tag{I}
$$

Here is the complete earlier proof, reproduced so that this release does not
depend on an unpublished review draft. Write $`G_j=G(j/K)`$. By the support
endpoint of XO, $`G_j=1`$ for every $`j\ge m`$. At time zero,
$`D_{0,j}=(j/K)(1-j/K)`$; since $`m/K\ge5/6>1/2`$, this profile is strictly
decreasing for $`j\ge m`$.

Inductively assume $`D_{t,j}=D_{0,j}`$ for every $`j\ge m`$. Then each
$`j>m`$ has $`D_{t,j}<D_{t,m}`$, so the selected prefix satisfies
$`s_t\le m`$. If $`j> s_t`$ and $`j\ge m`$, the cumulative update gives
$`F_{t+1,j}=F_{t,j}-M_t(1-G_j)=F_{t,j}`$. The only remaining case is
$`j=s_t=m`$, where
$`F_{t+1,m}=(1-\alpha)F_{t,m}+\alpha F_{t,m}G_m=F_{t,m}`$.
Thus all cumulative coordinates at and beyond $`m`$ are conserved, and the
induction closes. In particular,
$`F_{t,m}=F_{0,m}=2b-b^2=c`$ and
$`D_{t,m}=D_{0,m}=b(1-b)=a`$.

For the uniform bound, $`5/6\le b\le10/11`$: if $`K\ge11`$, use
$`\lceil5K/6\rceil\le5K/6+5/6`$; for $`K=6,\ldots,10`$, the values of $`b`$
are $`5/6,6/7,7/8,8/9,9/10`$. Since $`x(1-x)`$ decreases on
$`[1/2,1]`$, $`a=b(1-b)\ge10/121>10^{-12}`$. Therefore $`d_t\ge a`$
and the tolerance stopping rule never fires in the exact model. This is the
previously proved T1 dependency, reproduced without changing its status or
claiming it as a new S3 result.

We also use

$$
35/36\le c<1,\qquad h\le1/6.
\tag{B}
$$

The full cumulative update is

$$
D'_{j}=\begin{cases}
\frac34D_j+\frac14(F_sG_j-jh),&j\le s,\\
D_j-\frac14F_s(1-G_j),&j>s.
\end{cases}
\tag{U}
$$

## 3. P1 — finite boundary entry (T2)

For each interior index $`j\le m-2`$, put $`y=jh`$ and $`L_j=cG(y)-y`$. The decisive bound is

$$
L_j<a.
\tag{1}
$$

If $`y\le2/3`$, then $`G(y)\le y`$, so $`L_j\le(c-1)y\le0<a`$. Otherwise $`2/3<y<5/6`$. Put $`z=b-y=(m-j)h\ge2h`$. Since $`b-5/6=rh/6`$,

$$
\begin{aligned}
a-L_j
&=c(1-G(y))-(b-y)\\
&=(2c-1)z-crh/3\\
&\ge h[c(4-r/3)-2]\\
&\ge h[(35/36)(7/3)-2]=29h/108>0.
\end{aligned}
\tag{2}
$$

For a swept interior coordinate, $`F_s\le c`$ and (U) imply

$$
D'_{j}\le\tfrac34D_j+\tfrac14L_j.
\tag{3}
$$

For an unswept coordinate, (U) subtracts $`F_s(1-G_j)/4>0`$. Therefore, once $`D_j<a`$, it remains strictly below $`a`$. While $`D_j\ge a`$, it decreases by at least

$$
\tfrac14\min\{a-L_j,\ a(1-G_j)\}>0,
\tag{4}
$$

because $`F_s=sh+d\ge a`$. Each interior coordinate must cross below $`a`$ in finite time; there are only finitely many such coordinates. Consequently the stronger region

$$
D_j<a\quad\text{for every }j\le m-2
\tag{R}
$$

is reached in finite time and is forward-invariant. With (I), only $`\ell`$ and $`m`$ can then maximize discrepancy.

For the explicit bound, (2) and the first case give $`a-L_j\ge29h/108`$ for all interior indices. Also $`1-G_j\ge2h`$: for $`jh\le2/3`$, use $`1-G_j\ge1/3\ge2h`$; on the last linear segment, use $`1-G_j=2z-rh/3\ge7h/3`$. Thus the decrease in (4) is at least

$$
\frac14\min\{29h/108,\ 2ah\}\ge\frac{5}{121K}.
$$

Initially $`D_{0,j}\le1/4`$ and $`1/4-a\le81/484`$. Hence (R) holds by a time

$$
T\le\left\lfloor\frac{81K}{20}\right\rfloor+1,
\qquad t_*\le T.
\tag{5}
$$

This conservative bound proves entry, not the exact short entry times observed in S2.

## 4. P2 — persistence from first boundary selection (T3, K ≥ 7)

Entry into (R) may occur after the first boundary selection. To cover that first selection, consider the weaker condition

$$
D_j<\max(D_\ell,D_m)\quad(j\le m-2),
\tag{6}
$$

together with (I). If selection is $`m`$, then $`d=a`$, so (R) already holds and persists by §3.

If selection is $`\ell`$, write $`x=p_m`$. Then $`x\le h`$ and

$$
w=F_\ell=c-x\ge c-h\ge29/36.
$$

Both $`\ell`$ and every interior coordinate are swept. Define $`H(y)=wG(y)-y`$. Equation (U) gives

$$
D'_\ell-D'_j=\tfrac34(D_\ell-D_j)+\tfrac14[H(\ell h)-H(jh)].
\tag{7}
$$

On the XO segments through $`\ell h<5/6`$, the slopes of $`H`$ are $`w/2-1<0`$, $`3w/2-1>0`$, and $`2w-1>0`$. It decreases up to $`1/3`$ and increases thereafter. Thus it suffices to compare the grid endpoints:

$$
\begin{aligned}
H(\ell h)-H(h)
&=w(1-q-h/2)-b+2h\\
&\ge Z_K:=(c-h)(1-q-h/2)-b+2h.
\end{aligned}
\tag{8}
$$

For $`K=6n+r`$ and $`u=n-1`$, exact expansion yields

$$
\begin{aligned}
6K^3Z_K={}&180u^3+(231+136r)u^2\\
&+(12+152r+30r^2)u\\
&+2r^3+19r^2+31r-39.
\end{aligned}
\tag{9}
$$

For $`r=1,\ldots,5`$, the constant terms are respectively $`13,115,279,517,841`$, and all other coefficients are positive. For $`r=0`$ with $`K\ge7`$, $`u\ge1`$, so the expression is at least $`180+231+12-39=384>0`$. Hence $`Z_K>0`$ for every $`K\ge7`$.

It follows that $`H(\ell h)>H(jh)`$ for every interior index. By (7), $`D'_\ell>D'_j`$. Condition (6) is therefore forward-invariant for $`K\ge7`$ on the inherited invariant set.

At the first selection of $`\ell`$ or $`m`$, an interior tie at the maximum is impossible: the smaller interior index would have been selected. Thus (6) holds at $`t_*`$, proving persistence from first entry for $`K\ge7`$.

## 5. K = 6 exceptional bridge (completion of P2/T3)

The endpoint bound (8) is negative at $`K=6`$; the proof does not apply it there. Instead the canonical trajectory has an exact two-update bridge. Here $`m=5`$, $`\ell=4`$, $`a=5/36`$,

$$
p_0=(11,9,7,5,3,1)/36,\qquad
g=(1/12,1/12,1/4,1/4,1/3,0),\qquad s_0=3.
$$

Direct rational updates give:

| State time | $`D_1`$ | $`D_2`$ | $`D_3`$ | $`D_4`$ | $`D_5=a`$ | Selection |
|---|---|---|---|---|---|---|
| 1 | $`5/64`$ | $`11/96`$ | $`9/64`$ | $`23/144`$ | $`5/36`$ | $`4=\ell`$ |
| 2 | $`59/1728`$ | $`1/27`$ | $`115/1728`$ | $`157/1728`$ | $`5/36`$ | $`5=m`$ |

Time 1 is the first boundary entry and satisfies (6). At time 2, all interior coordinates satisfy (R), which persists by §3. This completes P2 for the canonical trajectory for every $`K\ge6`$.

## 6. P3 — boundary-coordinate scalar closure (T4)

For every $`t\ge t_*`$, the maximizing set is contained in $`\{\ell,m\}`$, and

$$
D_{t,\ell}=a+h-x_t.
$$

Smallest-index tie handling therefore gives

$$
s_t=\begin{cases}m,&x_t>h,\\\ell,&x_t\le h,\end{cases}
\qquad d_t=a+(h-x_t)^+.
\tag{10}
$$

Sweeping through $`m`$ removes a quarter of $`x_t`$ and redistributes total mass $`c/4`$. Sweeping through $`\ell`$ leaves bin $`m`$ unswept and redistributes total mass $`(c-x_t)/4`$. Hence

$$
x_{t+1}=\begin{cases}
\frac34x_t+\frac14cq,&x_t>h,\\
x_t+\frac14(c-x_t)q,&x_t\le h.
\end{cases}
\tag{11}
$$

The moved mass is $`c/4`$ on the upper branch and $`(c-x_t)/4`$ on the lower branch. Equations (10)–(11) determine the boundary coordinate, selection, discrepancy, boundary ties, moved mass, and discrepancy increments. They do not reconstruct the full mass vector.

Threshold equality must be distinguished from strict deficit:

| Scalar situation | Discrepancy change |
|---|---|
| $`x_t>h`$, $`x_{t+1}<h`$ | Strict rebound |
| $`x_t>h`$, $`x_{t+1}\ge h`$ | Exact tie |
| $`x_t<h`$ | Strict contraction |
| $`x_t=h`$ | Exact tie; selection is $`\ell`$ |

At equality, $`x_{t+1}=h+(c-h)q/4>h`$ but $`d_t=d_{t+1}=a`$. Thus the S2 shortcut “branch B strictly contracts” needs $`x_t<h`$ in a general statement. The S2 census recorded no such boundary-coordinate equality cases; its historical checks remain [VCR].

## 7. P4 — residue-class dichotomy (T5)

Within (11), the upper branch has fixed point $`L=cq`$. While it remains active,

$$
x_{t+n}=L+(3/4)^n(x_t-L).
\tag{12}
$$

The lower branch has fixed point $`c>h`$, with

$$
c-x_{t+1}=(1-q/4)(c-x_t).
\tag{13}
$$

If $`cq\ge h`$, every upper-branch state stays strictly above $`h`$. A lower-branch run must cross above $`h`$ in finite time, since (13) tends to $`c>h`$. Thereafter selection remains $`m`$ forever.

If $`cq<h`$, an upper-branch run cannot persist forever, since (12) tends below $`h`$. A lower-branch run cannot persist forever either, since (13) tends above $`h`$. Thus both $`x_t>h`$ and $`x_t\le h`$, and correspondingly both selected indices, occur infinitely often. A downward crossing here includes landing exactly at $`h`$, which selects $`\ell`$. No assertion that threshold equality is impossible is needed.

For canonical XO,

$$
\frac{cq}{h}=c\frac{6-r}{3}.
$$

If $`r\in\{0,1,2\}`$, this is at least $`(35/36)(4/3)=35/27>1`$. If $`r\in\{3,4,5\}`$, it is at most $`c<1`$. Therefore, for every integer $`K\ge6`$,

$$
\begin{aligned}
K\bmod6\in\{0,1,2\}&\Longrightarrow
\exists T\ \forall t\ge T:\ s_t=m,\\
K\bmod6\in\{3,4,5\}&\Longrightarrow
s_t=\ell\text{ and }s_t=m\text{ each infinitely often}.
\end{aligned}
\tag{14}
$$

This is an infinite-time branch dichotomy, not a periodicity result. The experimental labels LOCK/CYCLE retain their frozen finite-horizon meaning at $`H=10K`$. This proof does not additionally establish that the frozen classifier matches (14) at that cutoff for every untested $`K`$.

## 8. Audit provenance, evidence hierarchy, and closure

- **Proof source:** Astra's completed “Stage 4 S3 — Bounded XO Boundary-Reduction Proof Attempt” response in this chat, before this repository deposit round. The present document transcribes its argument and scope.
- **Original verification pass:** the proof response recorded an alternative tail-mass derivation of the interior bound, exact coefficient checks of the endpoint polynomial, exact $`K=6`$ arithmetic, and a separate threshold-tie check. These are distinct from the independent audit.
- **[Independent Audit]:** GPT-6.1 Sol, **PASS — no substantive gap found**, as supplied by the Owner in the 2026-09-30 deposit request. The [audit record](S3_XO_BOUNDARY_THEOREM_AUDIT.md) preserves source, scope, and the limit of the available transcript; this deposit does not claim to rerun that audit.
- **[Existing Math]:** T1, the original XO invariant, remains the dependency in §2.
- **[New Mathematical Result]:** T2–T5, supported by the deposited P1–P4 proof in §§3–7.
- **[VCR]:** S1/S2 exact computational evidence remains discovery and validation support, not the basis of universal quantification.

**S3: CLOSED — PASS.** No P1–P4 blocker remains within the audited canonical scope. Exact first-entry times, periodicity, full-state asymptotics, limiting dwell statistics, other initial states/laws/removal fractions, and 2D extensions are not settled by this theorem. No additional research is initiated by the deposit.
