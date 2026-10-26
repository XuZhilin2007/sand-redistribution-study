# -*- coding: utf-8 -*-
"""v04-A scan analysis: per-K critical alpha thresholds (float evidence),
alpha*argmin organisation on the left side, and the common-window envelope.

Reads scan_phase.json (float, exploratory). Prints:
  - per-K lower threshold (largest probed alpha < 1/4 where K fails)
  - per-K upper threshold (smallest probed alpha > 1/4 where K fails)
  - envelope brackets of the common two-round window
  - alpha*argmin stability on the left side
"""
from __future__ import annotations

import json
from pathlib import Path

rows = json.loads(Path(__file__).with_name("scan_phase.json").read_text())
left = [r for r in rows if r["alpha"] < 0.25]
right = [r for r in rows if r["alpha"] > 0.25]

# --- per-K thresholds -------------------------------------------------
ks = sorted({r["K"] for r in rows})
lo_thr, hi_thr = {}, {}
for K in ks:
    fails_lo = [r["alpha"] for r in left if r["K"] == K and r["viol"]]
    fails_hi = [r["alpha"] for r in right if r["K"] == K and r["viol"]]
    lo_thr[K] = max(fails_lo) if fails_lo else None   # K fails at <= this
    hi_thr[K] = min(fails_hi) if fails_hi else None

left_failers = sorted(K for K in ks if lo_thr[K] is not None)
right_failers = sorted(K for K in ks if hi_thr[K] is not None)
print(f"K range probed: {min(ks)}..{max(ks)}  n={len(ks)}")
print(f"\nLEFT side (alpha<1/4): {len(left_failers)} K ever fail")
print("  K : last failing alpha   (float)")
for K in left_failers:
    print(f"  {K:>5} : {lo_thr[K]:.6f}")
print(f"\nRIGHT side (alpha>1/4): {len(right_failers)} K ever fail; "
      f"smallest 12 by threshold:")
for K in sorted(right_failers, key=lambda K: hi_thr[K])[:12]:
    print(f"  {K:>5} : {hi_thr[K]:.6f}")

# --- common window envelope -------------------------------------------
inner_alphas = sorted({r["alpha"] for r in rows})
clean = [a for a in inner_alphas
         if not any(r["viol"] or r["flag"] for r in rows if r["alpha"] == a)]
below = [a for a in clean if a < 0.25]
above = [a for a in clean if a > 0.25]
lo_env = max(below)   # largest clean alpha below 1/4
hi_env = min(above)   # smallest clean alpha above 1/4
lo_fail_near = min((lo_thr[K] for K in left_failers), default=None)
hi_fail_near = min((hi_thr[K] for K in right_failers), default=None)
print(f"\nCOMMON WINDOW (probed K set, float):")
print(f"  bracket: ({max(lo_fail_near, lo_env):.6f} .. {lo_env:.6f}]"
      f" fail-free up to fail at {hi_fail_near:.6f}")
print(f"  i.e. clean at least [{lo_env}, {hi_env}] around 1/4")

# --- alpha * argmin on the left ---------------------------------------
print("\nalpha * argmin (left side, failing rows only, K<=30):")
hdr = "K   " + "  ".join(f"{a:.3f}" for a in sorted({r['alpha'] for r in left})[:9])
print(hdr)
for K in ks:
    if K > 30:
        break
    cells = []
    for a in sorted({r["alpha"] for r in left})[:9]:
        r = next(r for r in left if r["K"] == K and r["alpha"] == a)
        cells.append(f"{a * r['argmin']:5.2f}" if r["viol"] else "    .")
    print(f"{K:>3} " + " ".join(cells))
print("\nalphas:", [round(a, 4) for a in sorted({r['alpha'] for r in left})[:9]])
