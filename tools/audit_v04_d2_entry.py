"""Check the independent short proof of D2 entry using exact scalar bounds.

The all-K proof is in D2_ENTRY_REVIEW.md; this checks its rational arithmetic
and binds the accepted small-K dependency.  No generator is imported.
"""
from fractions import Fraction as Q
import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REVIEW = ROOT / "experiments/v04_all_k/review_2026-10-23"


def require(condition, label):
    if not condition:
        raise ValueError(label)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path,
                        help="write a new audit report; default only prints results")
    args = parser.parse_args()
    alo, ahi, h = Q(1999, 8000), Q(2001, 8000), Q(1, 33)
    rlo, rhi = 1 - ahi, 1 - alo
    lamlo, lamhi, nulo, nuhi = Q(9, 50), Q(19, 100), Q(19, 100), Q(21, 100)
    require(alo * Q(800, 1089) >= lamlo and ahi * Q(3, 4) <= lamhi, "lambda bounds")
    ulo, uhi = Q(1, 2) + 3 * lamlo / 4 - h / 2, Q(1, 2) + 3 * lamhi / 4 + h / 2
    require(Q(61, 100) < ulo <= uhi < Q(33, 50) < Q(2, 3), "true u1 bounds")
    p = lambda y: 2 * y - y * y
    f1lo = p(Q(61, 100)) - lamhi * (Q(4, 3) - 3 * Q(61, 100) / 2)
    f1hi = p(Q(33, 50)) - lamlo * (Q(4, 3) - 3 * Q(33, 50) / 2)
    require(nulo < alo * f1lo and ahi * f1hi < nuhi, "nu bounds")
    mulo, muhi = lamlo + nulo, lamhi + nuhi
    blo, bhi = rlo * lamlo + nulo, rhi * lamhi + nuhi
    coeffs = {"A": (rlo * rlo, 2 * rhi * rhi + bhi / 2 - 1, Q(0)),
              "B": (rlo * rlo, 2 * rhi * rhi + 3 * bhi / 2 - 1, -blo / 3),
              "C": (rlo, 2 * rhi + 3 * bhi / 2 - 1, -4 * rlo * lamlo / 3 - nulo / 3)}
    maxima = {}
    for label, (aa, bb, cc) in coeffs.items():
        require(aa > 0, label + " concavity")
        maximum = cc + bb * bb / (4 * aa)
        require(maximum < Q(3, 25), label + " strict upper bound")
        maxima[label] = {"a": str(aa), "b": str(bb), "c": str(cc), "upper": str(maximum)}
    e_slope = 1 + 3 * mulo / 2 - Q(4, 3)
    e_max = Q(2, 9) - mulo / 3
    upper_slope_min, upper_slope_max = 2 * mulo - Q(2, 3), 2 * muhi - Q(1, 3)
    boundary = Q(5, 36) - upper_slope_max * h
    gap = boundary - Q(3, 25)
    require(e_slope > 0 and e_max < Q(3, 25), "E endpoint monotonicity")
    require(upper_slope_min > 0 and boundary == Q(247, 1980) and gap == Q(47, 9900),
            "upper segment and boundary gap")
    require(Q(5, 6) - h > Q(2, 3), "ell lies on upper segment")

    small = REVIEW / "finite_common_audit.json"
    record = json.loads(small.read_text(encoding="utf-8"))
    require(record["status"] == "VERIFIED" and not record["failures"], "finite audit status")
    require(set(range(6, 33)) <= set(record["verified_K"]), "finite complement coverage")
    require(record["common_interval"] == [str(alo), str(ahi)], "finite complement interval")
    dependencies = [ROOT / p for p in ("tools/audit_v04_interval_v3.py",
                    "experiments/v04_w_interval/interval_certs_v3.jsonl",
                    "docs/XO_ALPHA_EXTENSION_PROOF.md", "docs/S3_XO_BOUNDARY_THEOREM_PROOF.md")]
    recorded_hashes = {name.replace("\\", "/"): digest
                       for name, digest in record["source_sha256"].items()}
    for path in dependencies:
        require(hashlib.sha256(path.read_bytes()).hexdigest() ==
                recorded_hashes[path.relative_to(ROOT).as_posix()], "finite proof dependency: " + str(path))
    sources = dependencies + [Path(__file__).resolve(), small,
               ROOT / "docs/v04_research/D2_ENTRY_REVIEW.md",
               ROOT / "docs/v04_research/LARGE_K_PREFIX_LEMMA_REVIEW.md",
               ROOT / "experiments/v04_all_k/d2_entry.py",
               ROOT / "experiments/v04_all_k/d2_entry_certificate.json"]
    report = {"status": "VERIFIED", "claim": "D2 maximizers lie in {ell,m} for K>=33, alpha in J; all K>=6 enter by t=2",
              "J": [str(alo), str(ahi)], "analytic_K_min": 33,
              "small_complement": [6, 32], "rho": [str(rlo), str(rhi)],
              "lambda": [str(lamlo), str(lamhi)], "u1": [str(ulo), str(uhi)],
              "nu_derived": [str(alo * f1lo), str(ahi * f1hi)], "nu": [str(nulo), str(nuhi)],
              "mu": [str(mulo), str(muhi)], "quadratic_envelopes": maxima,
              "E_slope_min": str(e_slope), "E_upper": str(e_max),
              "upper_slope_min": str(upper_slope_min), "upper_slope_max": str(upper_slope_max),
              "boundary_lower": str(boundary), "mid_region_gap": str(gap),
              "analytic_dependencies": ["true selector prefix lemma", "two-update discrepancy identities",
                                        "pointwise quadratic upper envelopes and completion of square",
                                        "upper-segment monotonicity", "theorem A: suffix and residence",
                                        "accepted finite-complement interval witnesses"],
              "source_sha256": {path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                                for path in sources}}
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("VERIFIED: K>=33 D2 entry on all J; all-K entry by t=2 with finite complement;")
    print(f"mu in [{mulo},{muhi}], boundary >= {boundary}, mid-region gap >= {gap}")


if __name__ == "__main__":
    main()
