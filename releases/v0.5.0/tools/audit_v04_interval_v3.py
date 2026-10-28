"""Independent, exact audit of v04 interval witnesses, without the generator.

Only the standard library is used.  No production experiment or Sturm kernel is
imported.  Polynomial signs are certified by nonnegative rational Bernstein
coefficients, with exact subdivision and fail-closed termination.  Boundary
residence and absence of stopping are analytic dependencies on theorem A.

The default scope is the COMMON interval [1999/8000, 2001/8000].  --scope
recorded instead checks each witness's own, usually wider, recorded interval.
This script never writes the generator's checkpoint or changes its source.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from math import comb
from pathlib import Path
import subprocess
import time


ROOT = Path(__file__).resolve().parents[1]
TARGETS = tuple(range(6, 204)) + (300, 600, 999, 1000, 1001)
ALPHA0 = Q(1, 4)
COMMON = (Q(1999, 8000), Q(2001, 8000))
ENGINE = "v04-interval-v3"
AUDITOR = "v04-independent-bernstein-v1"


class AuditFailure(Exception):
    pass


def require(condition, message):
    if not condition:
        raise AuditFailure(message)


def trim(p):
    p = tuple(Q(x) for x in p)
    n = len(p)
    while n and not p[n - 1]:
        n -= 1
    return p[:n]


def add(p, q):
    return trim((p[i] if i < len(p) else 0) +
                (q[i] if i < len(q) else 0)
                for i in range(max(len(p), len(q))))


def scale(p, k):
    return trim(x * k for x in p)


def sub(p, q):
    return add(p, scale(q, -1))


def mul(p, q):
    if not p or not q:
        return ()
    result = [Q(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x:
            for j, y in enumerate(q):
                if y:
                    result[i + j] += x * y
    return trim(result)


def value(p, x):
    result = Q(0)
    for coefficient in reversed(p):
        result = result * x + coefficient
    return result


def total(polynomials):
    result = ()
    for p in polynomials:
        result = add(result, p)
    return result


def power(p, n):
    result = (Q(1),)
    for _ in range(n):
        result = mul(result, p)
    return result


@lru_cache(maxsize=None)
def bernstein_matrix(n):
    return tuple(tuple(Q(comb(k, j), comb(n, j)) for j in range(k + 1))
                 for k in range(n + 1))


def bernstein_coefficients(p, lo, hi):
    """p(lo + (hi-lo)u) = sum b[k] C(n,k) u^k (1-u)^(n-k)."""
    if not p:
        return (Q(0),)
    n = len(p) - 1
    lp = [lo ** i for i in range(n + 1)]
    wp = [(hi - lo) ** i for i in range(n + 1)]
    shifted = [sum((p[i] * comb(i, j) * lp[i - j]
                    for i in range(j, n + 1)), Q(0)) * wp[j]
               for j in range(n + 1)]
    return tuple(sum((shifted[j] * factor for j, factor in enumerate(row)), Q(0))
                 for row in bernstein_matrix(n))


def split_bernstein(b):
    """Exact de Casteljau subdivision at the midpoint."""
    row = tuple(b)
    left, right = [row[0]], [row[-1]]
    while len(row) > 1:
        row = tuple((row[i] + row[i + 1]) / 2 for i in range(len(row) - 1))
        left.append(row[0])
        right.append(row[-1])
    return tuple(left), tuple(reversed(right))


def certify_sign(p, lo, hi, strict=False, left_open=False, max_depth=12):
    """Return (proved, leaves, deepest); inability to prove is never acceptance.

    Nonnegative coefficients imply nonnegativity on the closed interval.
    If at least one coefficient is positive, every interior value is positive.
    Checking the included endpoints therefore also certifies strict positivity.
    """
    def visit(b, open_left, depth):
        if strict:
            if b[-1] <= 0 or (not open_left and b[0] <= 0):
                return False, 0, depth
        elif b[0] < 0 or b[-1] < 0:
            return False, 0, depth
        if all(x >= 0 for x in b):
            return True, 1, depth
        if depth >= max_depth:
            return False, 0, depth
        first, second = split_bernstein(b)
        r1 = visit(first, open_left, depth + 1)
        if not r1[0]:
            return r1
        r2 = visit(second, False, depth + 1)
        return r2[0], r1[1] + r2[1], max(r1[2], r2[2])

    require(lo < hi, "invalid sign interval")
    return visit(bernstein_coefficients(trim(p), lo, hi), left_open, 0)


def xo_cdf(y):
    if y <= Q(1, 3):
        return y / 2
    if y <= Q(2, 3):
        return Q(3, 2) * y - Q(1, 3)
    if y <= Q(5, 6):
        return 2 * y - Q(2, 3)
    return Q(1)


def update_full(state, selected, weights):
    transferred = mul((Q(0), Q(1)), total(state[:selected]))
    return [add(mul((Q(1), Q(-1)), p) if i < selected else p,
                scale(transferred, weights[i]))
            for i, p in enumerate(state)]


def full_loss(state, k):
    return scale(total(mul(sub(p, (Q(1, k),)), sub(p, (Q(1, k),)))
                       for p in state), k)


class PieceAudit:
    def __init__(self, k, lo, hi, left_open, word, horizon):
        self.k, self.lo, self.hi = k, lo, hi
        self.left_open = left_open
        self.word, self.horizon = word, horizon
        self.obligations = 0
        self.leaves = 0
        self.max_depth = 0
        self.max_degree = 0
        self.digest = hashlib.sha256()

    def sign(self, p, strict, label):
        p = trim(p)
        ok, leaves, depth = certify_sign(p, self.lo, self.hi,
                                        strict, self.left_open)
        require(ok, f"K={self.k}, {label}: sign not certified")
        self.obligations += 1
        self.leaves += leaves
        self.max_depth = max(self.max_depth, depth)
        self.max_degree = max(self.max_degree, len(p) - 1)
        payload = [label, strict, str(self.lo), str(self.hi), self.left_open,
                   [str(x) for x in p]]
        self.digest.update(json.dumps(payload, separators=(",", ":")).encode())
        self.digest.update(b"\n")

    def license_global(self, state, selected, label):
        prefix, discrepancies = (), []
        for j, p in enumerate(state, 1):
            prefix = add(prefix, p)
            discrepancies.append(sub(prefix, (Q(j, self.k),)))
        chosen = discrepancies[selected - 1]
        for j, other in enumerate(discrepancies, 1):
            if j != selected:
                self.sign(sub(chosen, other), j < selected, f"{label}:s>{j}")

    def run(self):
        k, word, t1 = self.k, self.word, self.horizon
        require(type(t1) is int and 3 <= t1 <= 500, "invalid T1")
        require(isinstance(word, list) and len(word) == t1, "incomplete word")
        require(all(type(s) is int and 1 <= s <= k for s in word), "invalid selector")
        require(0 < self.lo < self.hi < 1, "alpha outside (0,1)")
        m = (5 * k + 5) // 6
        h, b = Q(1, k), Q(m, k)
        c = 2 * b - b * b
        weights = [xo_cdf(Q(i, k)) - xo_cdf(Q(i - 1, k))
                   for i in range(1, k + 1)]
        q = weights[m - 1]
        require(q == Q(6 - k % 6, 3 * k), "boundary weight identity")
        require(sum(weights) == 1 and min(weights) >= 0, "invalid XO weights")
        require(all(g == 0 for g in weights[m:]), "nonzero suffix weight")
        require(c - b > Q(1, 10 ** 12), "stopping barrier missing")
        p0 = [(Q(2 * k - 2 * i + 1, k * k),) for i in range(1, k + 1)]
        require(total(p0[:m]) == (c,), "canonical prefix identity")
        states = [p0]
        for t in (0, 1):
            require(word[t] <= m, "selector leaves invariant prefix")
            self.license_global(states[t], word[t], f"global t={t}")
            states.append(update_full(states[t], word[t], weights))
        for state in states:
            require(total(state) == (Q(1),), "mass identity")
            require(total(state[:m]) == (c,), "prefix mass identity")
            require(state[m:] == p0[m:], "frozen suffix identity")
        t0 = 1 if k == 6 else 2
        require(word[t0] in (m - 1, m), "no boundary entry at anchor")
        if t0 == 2:
            self.license_global(states[2], word[2], "global t=2 entry")
        u0, u1, u2 = [full_loss(state, k) for state in states]
        self.sign(sub(u0, u2), True, "prefix U0-U2")
        self.sign(sub(u1, u2), True, "prefix U1-U2")

        # Vector expansion p = base + x * direction + z * error.
        # At the boundary coordinate direction=1 and error=base=0;
        # above the boundary the canonical suffix is fixed.
        x0 = states[t0][m - 1]
        gamma = scale(sub((c,), x0), 1 / (1 - q))
        error = [sub(states[t0][i], scale(gamma, weights[i])) for i in range(m - 1)]
        require(not total(error), "error vector must have zero sum")
        base = [c * g / (1 - q) for g in weights[:m - 1]] + [Q(0)]
        base += [p[0] for p in p0[m:]]
        direction = [-g / (1 - q) for g in weights[:m - 1]] + [Q(1)]
        direction += [Q(0)] * (k - m)
        error += [()] * (k - m + 1)
        require(sum(base) == 1 and sum(direction) == 0, "affine mass identity")
        # Compare coefficients of formal independent variables x and z after
        # the full update, for BOTH boundary branches.  This checks the vector
        # dynamics identity rather than inferring it from a finite trajectory.
        rho = (Q(1), Q(-1))
        scalar_shift = (Q(0), c * q)
        for selected in (m - 1, m):
            delta = int(selected == m - 1)  # transferred mass alpha*(c-delta*x)
            require(sum(base[:selected]) == c and
                    sum(direction[:selected]) == -delta and not total(error[:selected]),
                    "formal selected prefix mass c-delta*x")
            scalar_factor = (Q(1), -q if delta else Q(-1))
            for i in range(k):
                shrink = rho if i < selected else (Q(1),)
                next_base = add(scale(shrink, base[i]), (Q(0), c * weights[i]))
                next_direction = add(scale(shrink, direction[i]),
                                     (Q(0), -delta * weights[i]))
                require(next_base == add((base[i],), scale(scalar_shift, direction[i])),
                        "two-branch identity: constant coefficient")
                require(next_direction == scale(scalar_factor, direction[i]),
                        "two-branch identity: x coefficient")
                require(mul(shrink, error[i]) == mul(rho, error[i]),
                        "two-branch identity: z coefficient")
        offset = [x - h for x in base]
        a = k * sum(x * x for x in direction)
        bb = 2 * k * sum(x * y for x, y in zip(offset, direction))
        cc = k * sum(x * x for x in offset)
        dd = scale(total(scale(e, x) for x, e in zip(direction, error)), k)
        ff = scale(total(scale(e, x) for x, e in zip(offset, error)), k)
        ee = scale(total(mul(e, e) for e in error), k)
        require(a > 0, "loss quadratic coefficient A must be positive")

        def reduced_loss(x, z):
            return total((scale(mul(x, x), a), scale(x, bb), (cc,),
                          scale(mul(z, add(mul(dd, x), ff)), 2),
                          mul(ee, mul(z, z))))

        # Theorem A supplies permanent boundary residence.  The following
        # vector identity holds for either boundary branch; scalar replay
        # licenses the precise selector, including its equality convention.
        x, z = x0, (Q(1),)
        for t in range(t0, t1 + 1):
            ut = reduced_loss(x, z)
            if t <= 2:
                require(ut == full_loss(states[t], k), f"full/reduced U{t} identity")
            if t >= 3:
                self.sign(sub(ut, u2), True, f"finite U{t}-U2")
            if t == t1:
                break
            selected = word[t]
            require(selected in (m - 1, m), f"nonboundary word at t={t}")
            if selected == m:
                self.sign(sub(x, (h,)), True, f"boundary t={t}:x>h")
                x = add(mul((Q(1), Q(-1)), x), (Q(0), c * q))
            else:
                self.sign(sub((h,), x), False, f"boundary t={t}:x<=h")
                x = add(mul((Q(1), -q), x), (Q(0), c * q))
            z = mul(z, (Q(1), Q(-1)))

        # A(U-U2) = (A*x + B/2 + d*z)^2 + S2*z^2 + AH1*z + AH0.
        # S2=A*E-d^2 >= 0 by Cauchy-Schwarz on direction,error.
        # Their boundary components (1,0) also give the stronger Gram form.
        s2 = sub(scale(ee, a), mul(dd, dd))
        error_norm2 = scale(ee, Q(1, k))
        pref_norm2 = sum(x * x for x in direction[:m - 1])
        inner = scale(dd, Q(1, k))
        gram = scale(sub(scale(error_norm2, 1 + pref_norm2), mul(inner, inner)), k * k)
        require(s2 == gram, "independent vector Gram identity")
        ah1 = sub(scale(ff, 2 * a), scale(dd, bb))
        ah0 = sub(scale(sub((cc,), u2), a), (bb * bb / 4,))
        w = power((Q(1), Q(-1)), t1 + 1 - t0)
        self.sign(scale(add(scale(mul(s2, w), 2), ah1), -1), False,
                  "tail derivative at w <=0")
        self.sign(total((mul(s2, mul(w, w)), mul(ah1, w), ah0)), True,
                  "tail completed-square bound at w >0")
        return {"K": k, "lo": str(self.lo), "hi": str(self.hi),
                "mode": "leftopen" if self.left_open else "closed", "t0": t0,
                "T1": t1, "tail_starts": t1 + 1, "status": "VERIFIED",
                "obligations": self.obligations, "bernstein_leaves": self.leaves,
                "max_subdivision_depth": self.max_depth, "max_degree": self.max_degree,
                "obligations_sha256": self.digest.hexdigest()}


def validate_rows(rows):
    require(all(isinstance(r, dict) and type(r.get("K")) is int for r in rows),
            "invalid record/K type")
    require(Counter(r["K"] for r in rows) == Counter(TARGETS),
            "target set must be covered exactly once")
    for r in rows:
        require(r.get("engine") == ENGINE and r.get("class") == "PASS",
                f"K={r['K']}: wrong engine/class")
        require(r.get("mode") in ("closed", "tie-split"), "unsupported mode")
        lo, hi = Q(r["lo"]), Q(r["hi"])
        require(0 < lo < ALPHA0 < hi < 1, "invalid recorded bounds")
        require(lo <= COMMON[0] < COMMON[1] <= hi, "common interval not covered")
        if r["mode"] == "tie-split":
            require(isinstance(r.get("T1"), list) and len(r["T1"]) == 2,
                    "split T1 missing")
            require(isinstance(r.get("word_head"), dict) and
                    set(r["word_head"]) == {"left", "right"}, "split words missing")


def audit_row(row, scope):
    lo, hi = (Q(row["lo"]), Q(row["hi"])) if scope == "recorded" else COMMON
    k = row["K"]
    if row["mode"] == "closed":
        return [PieceAudit(k, lo, hi, False, row["word_head"], row["T1"]).run()]
    words = row["word_head"]
    return [PieceAudit(k, lo, ALPHA0, False, words["left"], row["T1"][0]).run(),
            PieceAudit(k, ALPHA0, hi, True, words["right"], row["T1"][1]).run()]


def self_test():
    tests = []
    def check(label, condition):
        require(condition, f"auditor self-test failed: {label}")
        tests.append(label)
    def sign(p, lo, hi, strict=False, left_open=False):
        return certify_sign(p, lo, hi, strict, left_open)[0]
    check("positive constant", sign((Q(1),), Q(0), Q(1), True))
    check("zero weak", sign((), Q(0), Q(1)))
    check("zero strict rejected", not sign((), Q(0), Q(1), True))
    check("negative constant rejected", not sign((Q(-1),), Q(0), Q(1)))
    linear = (-ALPHA0, Q(1))
    check("left-open strict endpoint root", sign(linear, ALPHA0, Q(1, 2), True, True))
    check("closed strict endpoint root rejected", not sign(linear, ALPHA0, Q(1, 2), True))
    check("included right root strict rejected", not sign(scale(linear, -1), Q(0), ALPHA0, True, True))
    check("interior tangent accepted weakly", sign(mul(linear, linear), Q(0), Q(1, 2)))
    sq = mul(linear, linear)
    old_counterexample = mul(sq, sub(sq, (Q(1, 1600),)))
    check("review I1 negative region rejected", not sign(old_counterexample, Q(1, 5), Q(3, 10)))
    eps = Q(1, 10 ** 40)
    narrow = mul(linear, (-ALPHA0 - eps, Q(1)))
    check("close-root exact midpoint negative", value(narrow, ALPHA0 + eps / 2) < 0)
    check("close-root weak closed rejected", not sign(narrow, ALPHA0, COMMON[1]))
    check("close-root weak left-open rejected", not sign(narrow, ALPHA0, COMMON[1], False, True))
    # Check the basis identity independently at several rational points.
    p = (Q(1, 7), Q(-2, 3), Q(5, 11), Q(9, 13))
    lo, hi = Q(-2, 5), Q(7, 8)
    coeffs = bernstein_coefficients(p, lo, hi)
    for u in (Q(0), Q(1, 7), Q(2, 3), Q(1)):
        expanded = sum((b * comb(3, i) * u ** i * (1 - u) ** (3 - i)
                        for i, b in enumerate(coeffs)), Q(0))
        check(f"Bernstein identity at u={u}", expanded == value(p, lo + (hi - lo) * u))
    left, right = split_bernstein(coeffs)
    check("subdivision left basis identity", left == bernstein_coefficients(p, lo, (lo + hi) / 2))
    check("subdivision right basis identity", right == bernstein_coefficients(p, (lo + hi) / 2, hi))
    k6word = [3, 4] + [5] * 13
    check("valid K6 whole-time witness accepted",
          PieceAudit(6, *COMMON, False, k6word, 15).run()["status"] == "VERIFIED")
    def rejected_piece(lo, hi, word, horizon, expected):
        try:
            PieceAudit(6, lo, hi, False, word, horizon).run()
        except AuditFailure as exc:
            return expected in str(exc)
        return False
    check("false K6 wide whole-time claim rejected by U1 prefix",
          rejected_piece(Q(1, 5), Q(3, 10), k6word[:5], 5, "prefix U1-U2"))
    check("tampered initial selector rejected",
          rejected_piece(*COMMON, [4] + k6word[1:], 15, "global t=0"))
    check("tampered final selector rejected",
          rejected_piece(*COMMON, k6word[:-1] + [4], 15, "boundary t=14"))
    check("truncated word rejected",
          rejected_piece(*COMMON, k6word[:-1], 15, "incomplete word"))
    return tests


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path,
                        default=ROOT / "experiments/v04_w_interval/interval_certs_v3.jsonl")
    parser.add_argument("--scope", choices=("common", "recorded"), default="common")
    parser.add_argument("--report", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    started = time.monotonic()
    controls = self_test()
    print(f"Auditor controls: {len(controls)} PASS", flush=True)
    if args.self_test:
        return 0
    sources = [Path(__file__).resolve(), args.input.resolve(),
               ROOT / "experiments/v04_w_interval/interval_certs.py",
               ROOT / "experiments/v04_w_interval/symkernel.py",
               ROOT / "docs/XO_ALPHA_EXTENSION_PROOF.md",
               ROOT / "docs/S3_XO_BOUNDARY_THEOREM_PROOF.md"]
    source_hashes = {str(p.relative_to(ROOT)): file_hash(p) for p in sources}
    git_head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                       text=True).strip()
    input_bytes = args.input.read_bytes()
    require(hashlib.sha256(input_bytes).hexdigest() ==
            source_hashes[str(args.input.resolve().relative_to(ROOT))], "input changed while reading")
    rows = [json.loads(line) for line in input_bytes.decode("utf-8").splitlines()
            if line.strip()]
    validate_rows(rows)
    pieces, failures = [], []
    for count, row in enumerate(rows, 1):
        try:
            pieces.extend(audit_row(row, args.scope))
        except (AuditFailure, KeyError, ValueError, TypeError, ArithmeticError) as exc:
            failures.append({"K": row["K"], "problem": str(exc)})
            print(f"FAIL K={row['K']}: {exc}", flush=True)
        if count % 10 == 0 or count == len(rows):
            print(f"Audited {count}/{len(rows)} K; failures={len(failures)}; "
                  f"elapsed={time.monotonic() - started:.1f}s", flush=True)
    require(source_hashes == {str(p.relative_to(ROOT)): file_hash(p) for p in sources},
            "source/input/dependency changed during audit")
    require(git_head == subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                               text=True).strip(), "git HEAD changed during audit")
    report = {"auditor": AUDITOR, "scope": args.scope,
              "status": "VERIFIED" if not failures else "NOT_VERIFIED",
              "git_head": git_head, "source_sha256": source_hashes,
              "targets": list(TARGETS), "target_count": len(TARGETS),
              "verified_K": sorted({p["K"] for p in pieces} - {f["K"] for f in failures}),
              "piece_count": len(pieces), "common_interval": [str(x) for x in COMMON],
              "analytic_dependencies": ["theorem A: canonical no-stopping/prefix invariants",
                                        "theorem A: permanent boundary residence/scalar map",
                                        "two-branch vector identity p=base+x*direction+z*error",
                                        "Cauchy-Schwarz: A*E-d^2>=0",
                                        "completed-square monotone tail bound",
                                        "Bernstein basis nonnegativity/subdivision identity"],
              "controls": controls, "pieces": pieces, "failures": failures,
              "elapsed_seconds": round(time.monotonic() - started, 3)}
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"{report['status']}: {len(report['verified_K'])}/{len(TARGETS)} K, "
          f"{len(pieces)} pieces, {sum(p['obligations'] for p in pieces)} exact sign obligations",
          flush=True)
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
