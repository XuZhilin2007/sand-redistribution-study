# Literature ledger and exact theory interfaces for Sand v0.4.0

Search and reading dates: 26–27 October 2026, Asia/Shanghai. This is a
structured literature review with selective reading of the most relevant
theorems. It is not an exhaustive mathematical novelty search. A reference
below is credited only at the access level stated. No conclusion relies on
AI review verdicts, search-result absence, citation counts, or third-party
paper summaries.

The relationship classifications used here are: **direct consequence** when
a cited theorem applies after a checked translation; **specialization** when
Sand instantiates an established object or elementary method; **technical
extension** when Sand supplies additional explicit or uniform information;
**potential independent contribution** when no implication has been established;
and **unresolved** when an interface or assumption remains unverified. These
describe this review's evidence, rather than priority in the entire literature.

**The exact interface that matters most.** After boundary residence, Sand has

\[
 A(x)=(1-\alpha)x+\alpha cq\quad(x>h),\qquad
 B(x)=(1-\alpha q)x+\alpha cq\quad(x\le h),
\]

with h=1/K, c=2b−b², b=⌈5K/6⌉/K and q=(6−K mod 6)/(3K). In a CYCLE class
cq<h<c. Put \(u_-=A(h)\), \(u_+=B(h)\). The interval
\((u_-,u_+]\) is invariant: each increasing branch moves toward its fixed
point on the opposite side of h, and its endpoint images stay in this interval.
Each orbit in (0,c] enters it after finitely many threshold crossings. Under
the orientation reversing change \(z=(u_+-x)/(u_+-u_-)\), the map is

\[
 f(z)=\begin{cases}
 \lambda z+\delta,&0\le z<\eta,\\
 \mu(\lambda z+\delta-1),&\eta\le z<1,
 \end{cases}\qquad
 \lambda=1-\alpha,\quad
 \mu=\frac{1-\alpha q}{1-\alpha},\quad
 \eta=\frac{q(c-h)}{(1-q)h},\quad
 \delta=1-\lambda\eta.
\]

These are Sand's own algebraic identities, independently rechecked here. The
interval orientation and assignment of x=h to B matter: z=η maps to 0.
For every CYCLE K and 0<α<1, 0<η<1, 0<λ<1, λμ<1 and 1−λ<δ<1.
They meet the standing hypotheses of Laurent–Nogueira's three-parameter
family [L1]. For algebraic α the three parameters are algebraic. Rational
rotation number then follows from [L1, Theorem 2]. Actual periodic-orbit
attraction additionally requires its non-exceptional endpoint case. The
explicit words and exclusion of the exceptional endpoint at canonical α=1/4
remain work supplied by Sand's certificates. This interface does not determine
the early loss U₂ or compare it with other times.

**L1. Michel Laurent and Arnaldo Nogueira, _Dynamics of 2-interval piecewise
affine maps and Hecke-Mahler series_.** Journal of Modern Dynamics 17 (2021),
33–63, [DOI 10.3934/jmd.2021002](https://www.aimsciences.org/article/doi/10.3934/jmd.2021002).
Access: publisher metadata and [arXiv:1907.08655v1 full text](https://arxiv.org/html/1907.08655v1),
Definition 1, Theorems 1–3 and §8. Relationship: **direct consequence** for
algebraic-parameter rationality after the dictionary above; conditional
specialization for periodic attraction. Theorem 3 distinguishes a rational
plateau's right endpoint, where the original map need not have a periodic
orbit. The presence of two unequal slopes is already covered here. Sand's
all-K transient minimum is not a conclusion of these theorems.

**L2. José Pedro Gaivão, Michel Laurent and Arnaldo Nogueira, _Rotation number
of 2-interval piecewise affine maps_.** Aequationes Mathematicae 99 (2025),
511–530; first published online 22 April 2024,
[DOI 10.1007/s00010-024-01064-2](https://link.springer.com/article/10.1007/s00010-024-01064-2).
Access: open publisher full text, §§2–6, Theorem 5.2, Corollary 5.3 and
Theorem 5.5; also [arXiv:2212.10944](https://arxiv.org/abs/2212.10944).
The five-parameter family includes the three-parameter case with the paper's
b=1, c=0; these names are not Sand's b,c. Relationship: **known extension**
of L1's framework, not a missing generalization to be claimed by Sand. The
fixed-slope variation in its rotation theorem does not by itself prove a
monotonic phase diagram along Sand's α path, which changes both slopes.

**L3. Michel Laurent and Arnaldo Nogueira, _Rotation number of contracted
rotations_.** Journal of Modern Dynamics 12 (2018), 175–191,
[DOI 10.3934/jmd.2018007](https://www.aimsciences.org/article/doi/10.3934/jmd.2018007).
Access: publisher metadata and [arXiv:1704.05130v3 full text](https://arxiv.org/html/1704.05130v3),
introduction and main statements. The arXiv title contains “interval”. Its
map \(\{\lambda x+\delta\}\) has equal slopes. Sand's slopes differ because
0<q<1 and α>0. Relationship: **related antecedent**, with L1 providing the
appropriate direct theorem; the equal-slope result alone is insufficient.

**L4. Svante Janson and Anders Öberg, _A piecewise contractive dynamical
system and Phragmén's election method_.** Bulletin de la Société Mathématique
de France 147(3) (2019), 395–441,
[DOI 10.24033/bsmf.2787](https://numdam.org/articles/10.24033/bsmf.2787/).
Access: publisher record and [arXiv:1709.06398v1 full text](https://arxiv.org/html/1709.06398v1),
§§1–2 and classification statements. This treats the two choices of value at
the discontinuity of an equal-slope contracted rotation, with rotation and
symbolic descriptions. Relationship: **related specialization**, especially
for distinguishing invariant sets, limit sets and endpoint conventions. Its
election application and universal-cycle results do not establish Sand's
unequal-slope loss minimum.

**L5. Arnaldo Nogueira and Benito Pires, _Dynamics of piecewise contractions
of the interval_.** Ergodic Theory and Dynamical Systems 35(7) (2015),
2198–2215, [DOI 10.1017/etds.2014.16](https://doi.org/10.1017/etds.2014.16).
Access: publisher metadata and [arXiv:1206.5676v1 full text](https://arxiv.org/html/1206.5676v1),
introduction and main statements. Injectivity and contraction on n pieces
give an upper bound of n periodic orbits. The result does not guarantee
existence or convergence for every parameter. Sand's invariant interval U,
after normalization, is the relevant injective domain; its whole scalar
domain (0,c] is not injective. Relationship: **conditional specialization**;
L1 is sharper for the present two-branch family.

**L6. Arnaldo Nogueira, Benito Pires and Rafael A. Rosales, _Asymptotically
periodic piecewise contractions of the interval_.** Nonlinearity 27(7) (2014),
1603–1610, [DOI 10.1088/0951-7715/27/7/1603](https://doi.org/10.1088/0951-7715/27/7/1603).
Access: [arXiv:1310.5784v1 full text](https://arxiv.org/html/1310.5784v1),
Theorems 1.1, 2.1 and Example 4.1. The assumptions include fixed smooth
contractions with disjoint images and almost every partition parameter.
Relationship: **related generic theorem**, not an all-K proof for Sand's
countable, correlated parameter set. A null exceptional set can contain every
point of a specified countable set. Its example also distinguishes eventual
symbolic periodicity from existence of a periodic orbit.

**L7. Eleonora Catsigeras, Pierre Guiraud, Arnaud Meyroneinc and Edgardo
Ugalde, _On the asymptotic properties of piecewise contracting maps_.**
Dynamical Systems 31(2) (2016), 107–135,
[DOI 10.1080/14689367.2015.1068274](https://doi.org/10.1080/14689367.2015.1068274).
Access: publisher metadata and
[authors' 31 March 2014 manuscript](https://www.fing.edu.uy/~eleonora/Recopilacion/Archivos/PiecewiseContractions31March2014.pdf),
Theorem 3.1 and surrounding definitions. A positive distance between the
attractor and the discontinuity set ensures a finite union of periodic
orbits. Relationship: **conditional specialization** once that separation is
proved; it cannot supply the separation itself. Crossing h repeatedly is
different from an attractor containing h. This is useful for perturbation
research, not a theorem about U₂.

**L8. Sakshi Jain and Carlangelo Liverani, _Piecewise contractions_.**
Ergodic Theory and Dynamical Systems 45(5) (2025), 1503–1540,
[DOI 10.1017/etds.2024.78](https://doi.org/10.1017/etds.2024.78).
Access: [publisher full text](https://www.cambridge.org/core/journals/ergodic-theory-and-dynamical-systems/article/piecewise-contractions/24054E955E5DDBDA27C85F25B4268C74),
§§2–3, including its definitions of Markov and piecewise injective contractions.
The openness and stability results use specific contraction, partition and
topology hypotheses. Relationship: **unresolved for perturbations of Sand**.
The full Sand state has conserved suffix directions, so it is not automatically
a strict contraction on every piece. Applying a reduced theorem for each
fixed K would still not give a K-independent stability radius.

**L9. Albert Granados, Lluís Alsedà and Maciej Krupa, _The Period Adding and
Incrementing Bifurcations: From Rotation Theory to Applications_.** SIAM
Review 59(2) (2017), 225–292,
[DOI 10.1137/140996598](https://epubs.siam.org/doi/10.1137/140996598).
Access: publisher metadata and [arXiv:1407.1895v2 full text](https://arxiv.org/html/1407.1895v2),
§2 and Theorem 2.6. This provides established Farey/rotation/symbolic
mechanisms under branch and parameter-path assumptions. Relationship:
**established mechanism, incomplete direct interface** for Sand's K or α
paths. The required monotone offset variation has not been shown for those
paths; observed period tables alone do not verify that hypothesis.

**L10. Mario di Bernardo, Christopher J. Budd, Alan R. Champneys and Piotr
Kowalczyk, _Piecewise-smooth Dynamical Systems: Theory and Applications_.**
Applied Mathematical Sciences 163, Springer (2008),
[DOI 10.1007/978-1-84628-708-4](https://link.springer.com/book/10.1007/978-1-84628-708-4).
Access: publisher preview and contents only. Relationship: **background
framework**, not a verified implication. The discontinuity of Sand's feedback
map makes this terminology appropriate; no specific theorem from the book
was used to certify Sand or its originality in this review.

**L11. Igal Sason, _On Data-Processing and Majorization Inequalities for
f-Divergences with Applications_.** Entropy 21(10) (2019), article 1022,
[DOI 10.3390/e21101022](https://doi.org/10.3390/e21101022).
Access: [arXiv:1911.02436v1 full text](https://arxiv.org/html/1911.02436v1),
§I.A, §II.A and the chi-square specialization. Sand's U is exactly Pearson
chi-square divergence from the uniform law. A fixed-selector update is a
stochastic linear operator and contracts L1 differences; data processing
compares both laws after applying that same operator. Sand's uniform target
is generally not invariant under it. Relationship: **known metric and operator
tools**, without a monotonic-U or two-update conclusion. The state-dependent
choice of operator is an additional issue.

**L12. Mikael Johansson and Anders Rantzer, _Computation of piecewise
quadratic Lyapunov functions for hybrid systems_.** IEEE Transactions on
Automatic Control 43(4) (1998), 555–559,
[DOI 10.1109/9.664157](https://doi.org/10.1109/9.664157).
Access: [published paper text](https://people.ece.ubc.ca/~moishi/eece571m/articles/JohannsonRantzer98.pdf),
§§II–IV. It formulates piecewise quadratic stability certificates for
continuous-time affine vector fields on polyhedral cells. Sand is discrete
and its loss need not decrease. Relationship: **methodological precedent**,
not a directly applicable stability theorem or a source of t=2. Quadratic
expansion, Gram positivity and completing a square are standard tools.

**L13. Stephen Prajna and Ali Jadbabaie, _Safety Verification of Hybrid
Systems Using Barrier Certificates_.** HSCC 2004, LNCS 2993, 477–492,
[DOI 10.1007/978-3-540-24743-2_32](https://doi.org/10.1007/978-3-540-24743-2_32).
Access: [authors' full manuscript](https://web.mit.edu/~jadbabai/www/papers/hscc04_2.pdf),
definition of safety and Theorem 2, including discrete resets. Sand's
infinite-time comparison can be interpreted as avoiding the set U<U₂ after
augmenting the state with its decaying residual variable. Relationship:
**established verification viewpoint**. This is an interpretation of Sand's
certificate; a literal barrier satisfying that paper's conditions has not
been synthesized here. It supplies no model-dependent constants or prefix.

**L14. Mohamed Amin Ben Sassi, Sriram Sankaranarayanan, Xin Chen and Erika
Ábrahám, _Linear relaxations of polynomial positivity for polynomial
Lyapunov function synthesis_.** IMA Journal of Mathematical Control and
Information 33(3) (2016), 723–756,
[DOI 10.1093/imamci/dnv003](https://doi.org/10.1093/imamci/dnv003).
Access: author-institution metadata and [arXiv:1407.2952v1 full text](https://arxiv.org/html/1407.2952v1),
§3.1, Lemma 3.1 and Corollary 3.1. The nonnegative Bernstein basis forms a
partition of unity, bounding a polynomial by its coefficients. Relationship:
**specialization of a standard positivity method** in Sand's finite checker.
Neither Bernstein bounds nor subdivision should be claimed as new. Sand
supplies the relevant exact polynomial obligations and coverage.

**L15. Lionel Levine and Yuval Peres, _Strong Spherical Asymptotics for
Rotor-Router Aggregation and the Divisible Sandpile_.** Potential Analysis
30(1) (2009), 1–27,
[DOI 10.1007/s11118-008-9104-6](https://doi.org/10.1007/s11118-008-9104-6).
Access: [arXiv:0704.0688v5 full text](https://arxiv.org/html/0704.0688v5),
definition of the divisible sandpile and Proposition 1.2. Their toppling
retains capacity one at a lattice site and sends excess to neighbors; the
order-independent limiting odometer is central. Sand removes a fraction
from a globally selected prefix and uses a fixed respray kernel. Relationship:
**distinct deterministic redistribution model**. Shared sand terminology does
not transfer Abelian stabilization or potential-theoretic shape theorems.

**L16. Akshay Kashyap, Tamer Başar and R. Srikant, _Quantized consensus_.**
Automatica 43(7) (2007), 1192–1203,
[DOI 10.1016/j.automatica.2007.01.002](https://www.sciencedirect.com/science/article/pii/S0005109807001045).
Access: publisher abstract/introduction and
[authors' paper](https://tamerbasar.csl.illinois.edu/KashyapBasarSrikantAutoVol42Iss7_2007.pdf)
search-visible text; no complete proof audit. It studies integer node loads,
connected graphs and randomized averaging. Relationship: **related objective,
different assumptions**. Sand's real masses, biased respray and prescribed
prefix feedback do not meet that setup. Convergence of an averaging algorithm
is not an early minimum of a nonmonotone prescribed orbit.

**L17. Oliver Jenkinson, _Ergodic Optimization_.** Discrete and Continuous
Dynamical Systems 15(1) (2006), 197–224,
[DOI 10.3934/dcds.2006.15.197](https://www.aimsciences.org/article/doi/10.3934/dcds.2006.15.197).
Access: publisher abstract/introduction only. Its objective is an extremal
ergodic average over orbits or invariant measures. Relationship: **different
optimization problem**. Sand compares instantaneous losses at integer times
along one fixed orbit. This source does not imply that minimum or its
attainment at an early time.

**Search method and limits.** The model was reconstructed before searching.
The following query families and citation trails were actually used:

| Family | Representative queries / trail | Purpose |
|---|---|---|
| Interval dynamics | `contracted rotation map lambda delta devil staircase Laurent Nogueira`; exact titles of L1 and L2 | Find the closest family and read its standing hypotheses and endpoint cases |
| General contraction theory | `piecewise contractions asymptotically periodic maps interval Nogueira Pires 2015`; L1 references → L5/L6; L6 references and newer L8 context | Distinguish generic periodicity from a statement at every canonical K |
| Symbolic structure | `Granados Alsedà Krupa Period Adding`; L1 → Janson–Öberg | Check whether word patterns and frequency locking are established mechanisms |
| Energy and certification | `Csiszar Information-type measures 1967 pdf`; `Sason f-Divergence data processing 2019`; `Computation of piecewise quadratic Lyapunov functions`; `Linear relaxations of polynomial positivity`; `Safety verification ... barrier certificates` | Check metric inheritance, contraction premises and proof-method precedent |
| Orbit minimization | `finite-time orbit minimization piecewise affine`; `finite time minimum piecewise quadratic dynamical systems stopping`; `ergodic optimization Birkhoff averages Jenkinson 2006` | Separate the actual transient comparison from stability, control synthesis and average optimization |
| Redistribution | `redistribution maximal discrepancy deterministic dynamics`; `divisible sandpile Levine Peres 2009`; `quantized consensus Kashyap 2007` | Test alternative formulations beyond the project's terminology |

Gambaudo–Tresser, Brémont, Bugeaud–Conze and switched-server/polling-system
references were encountered through the primary bibliographies. They are
**leads, not content-verified sources of a Sand theorem** in this assessment.
The historical Csiszár citation was searched, but the original 1967 paper
was not retrieved in a sufficiently usable form; L11 supplies the accessible
primary mathematical treatment used here. Some publisher full texts were
restricted, and several accessible arXiv texts are earlier versions of the
published papers. The ledger states the versions actually read.

The connected search service was used first; one transport failure was
followed by the built-in web interface. Publisher and author/arXiv material
were used for mathematical claims. No subscription MathSciNet search,
exhaustive citation census, or systematic survey of every language and venue
was performed. The search did not establish an existing theorem implying
Sand's full uniform transient result; it also did not establish priority or
exhaustive novelty. The most useful next literature check is a specialist's
assessment of the explicit reduction and the uniform finite-time certificate.
