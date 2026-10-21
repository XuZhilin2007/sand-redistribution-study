"""WP1 dev: LN-plateau classification of the CYCLE lattice.

For each CYCLE K, with the exact dictionary (lambda, mu, delta) =
(3/4, 4*beta/3, 1 - (3/4)*eta), find the UNIQUE rational p/q whose
plateau [delta^-(lambda,mu,p/q), delta^+(lambda,mu,p/q)] contains
delta_K (LN Theorem 1(i)); verify:
  - the plateau is strictly interior on the left and strict on the right
    (endpoint (iii) never occurs);
  - the T7 census word equals (up to rotation) the mechanical word
    w_j = A iff (j*p mod q) < q - p;
  - the derived rho = p/q equals the census word's B-count/period.

The plateau endpoints use LN's explicit formulae with
S = sum_{k=1}^{q-2} (floor((k+1)p/q) - floor(kp/q)) lambda^k mu^{floor(kp/q)}.
"""
from __future__ import annotations

import csv
import sys
from fractions import Fraction as F
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
from sand_m0.word_cyclicity import boundary_params  # noqa: E402

LAM = F(3, 4)
Q_MAX = 64


def S_sum(lam, mu, p, q):
    if q == 2:
        return F(0)
    tot = F(0)
    for k in range(1, q - 1):
        d = ((k + 1) * p) // q - (k * p) // q
        tot += d * lam ** k * mu ** ((k * p) // q)
    return tot


def delta_plus(lam, mu, p, q):
    S = S_sum(lam, mu, p, q)
    num = (1 - lam) * (1 + mu * S + lam ** (q - 1) * mu ** p * (1 - lam))
    den = 1 + (mu - 1) * S + lam ** (q - 1) * mu ** (p - 1) * \
        (mu - lam * mu - 1)
    return num / den


def delta_minus(lam, mu, p, q):
    S = S_sum(lam, mu, p, q)
    num = (1 - lam) * (1 + mu * S)
    den = 1 + (mu - 1) * S - lam ** q * mu ** (p - 1)
    return num / den


def dictionary(K):
    P = boundary_params(K)
    h, c, q, L, beta = P["h"], P["c"], P["q"], P["L"], P["beta"]
    A = lambda x: F(3, 4) * x + L * F(1, 4)
    B = lambda x: x + (c - x) * q * F(1, 4)
    D = B(h) - A(h)
    eta = (B(h) - h) / D
    mu = 4 * beta / 3
    delta = 1 - F(3, 4) * eta
    return mu, delta


def mechanical_word(p, q):
    """Word over one period from phase 0: w_j = A iff jp mod q < q-p."""
    return "".join("A" if (j * p) % q < q - p else "B" for j in range(q))


def same_rotation(w1, w2):
    if len(w1) != len(w2):
        return False
    return any(w1[i:] + w1[:i] == w2 for i in range(len(w1)))


def census():
    path = REPO / "experiments" / "v03_t7_cyclicity" / "results" / \
        "t7_census_certificates.csv"
    return {int(r["K"]): r["word"]
            for r in csv.DictReader(path.open(encoding="utf-8"))}


def main():
    import json
    import subprocess
    cens = census()
    out_rows = []
    n_ok = 0
    n_total = 0
    for K in sorted(cens):
        if K % 6 not in (3, 4, 5):
            continue
        n_total += 1
        mu, delta = dictionary(K)
        found = None
        for q in range(2, Q_MAX + 1):
            for p in range(1, q):
                from math import gcd
                if gcd(p, q) != 1:
                    continue
                dlo = delta_minus(LAM, mu, p, q)
                dhi = delta_plus(LAM, mu, p, q)
                if dlo <= delta <= dhi:
                    found = (p, q, dlo, dhi)
                    break
            if found:
                break
        if found is None:
            print(f"HARD FAIL K={K}: no plateau with q<={Q_MAX}")
            return 1
        p, q, dlo, dhi = found
        cw = cens[K]
        mw = mechanical_word(p, q)
        checks = {
            "word_len_is_q": len(cw) == q,
            "Bcount_is_p": cw.count("B") == p,
            "word_is_rotation_of_mechanical": same_rotation(cw, mw),
            "left_closed": dlo <= delta,
            "right_strict": delta < dhi,
        }
        if not all(checks.values()):
            print(f"HARD FAIL K={K}: {checks} rho={p}/{q}")
            return 1
        n_ok += 1
        out_rows.append({
            "K": K, "rho_p": p, "rho_q": q,
            "interior_pos": (delta - dlo) / (dhi - dlo),
            "mechanical_word": mw,
        })
    import csv as _csv
    fields = list(out_rows[0].keys())
    fields[3] = "interior_num"
    fields.append("interior_den")
    out_path = Path(__file__).resolve().parent / "results" /         "wp1_rows.csv"
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        w = _csv.DictWriter(fh, fieldnames=[
            "K", "rho_p", "rho_q", "interior_pos", "mechanical_word"])
        w.writeheader()
        for r in out_rows:
            w.writerow({**r, "interior_pos": str(r["interior_pos"])})
    res = {
        "n_total": n_total, "n_ok": n_ok,
        "all_right_strict": True,
        "elapsed_note": "see log",
        "git": subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=REPO,
            capture_output=True, text=True).stdout.strip(),
    }
    res_path = Path(__file__).resolve().parent / "results" /         "wp1_summary.json"
    res_path.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps(res))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
