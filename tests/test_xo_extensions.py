"""Fast exact regression for A/T6 and the declared scope boundaries.

Collected by tests/count_tests.py. The full matrices stay in the separate
experiments/xo_extensions/validate.py command. No archived results are read.
"""
from fractions import Fraction as F
import importlib.util
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments/xo_extensions"))
from exact_engine import frozen_classifier, integer_orbit, boundary_step


def load(relative):
    spec = importlib.util.spec_from_file_location(Path(relative).stem, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check(label, condition):
    if not condition:
        print("[FAIL]", label)
        raise AssertionError(label)
    print("[PASS]", label)


def main():
    if not __debug__:
        raise RuntimeError("Assertions are required; do not use -O")
    frozen = frozen_classifier()
    A = load("experiments/xo_alpha_extension/verify_alpha.py")
    T6 = load("experiments/xo_t6_classifier/verify_t6.py")
    D = load("experiments/xo_supporting_results/mine_descent.py")
    for K in range(6, 13):
        for alpha in (F(1, 8), F(1, 4), F(3, 4)):
            row, choices = integer_orbit(K, alpha, 10 * K, frozen)
            other = A.run_alpha(K, alpha, 10 * K)
            check(f"K={K} alpha={alpha}: full-mass orbit vs Fraction entry, confinement and frozen label",
                  row["t_star"] == other["t_star"] and row["class"] == other["finite_classifier"]["class"]
                  and row["R_m"] == other["finite_classifier"]["R_m"])
    for K in range(6, 18):
        row = T6.scalar_classifier(K)
        value = frozen(row["s_list"], (5 * K + 5) // 6, K)
        check(f"K={K}: frozen H=10K classifier vs residue dichotomy",
              value["class"] == ("LOCK" if K % 6 < 3 else "CYCLE") and value["R_m"] == row["R_m"])
    for K in (6, 9, 10):
        h = F(1, K)
        s, d, next_x = boundary_step(K, F(1, 4), h)
        m, _, _, c, a, _, q = A.anchor(K)
        check(f"K={K}: threshold equality chooses lower; discrepancy stays tied",
              s == m - 1 and d == a and next_x == h + (c - h) * q / 4 and next_x > h)
    for alpha, expected in ((F(4, 27), 3), (F(4, 27) + F(1, 1000000), 4),
                            (F(1, 3), 4), (F(1, 3) + F(1, 1000000), 5)):
        _, choices = integer_orbit(6, alpha, 3, frozen)
        check(f"K=6 alpha={alpha}: smallest-index choice at initial branch boundary", choices[1] == expected)
    for K in (9, 10, 11, 15, 16, 17):
        row = T6.scalar_classifier(K)
        dwell = T6.dwell_summary(row)
        check(f"K={K}: initial and subsequent dwell bounds remain separate",
              dwell["initial_dwell"] <= 3 * K and dwell["maxA_after_B"] <= dwell["boundA_after_B"]
              and dwell["maxB_after_A"] <= dwell["boundB_after_A"])
    m = 8
    check("LOCK precedence over two-branch post-set", frozen([m-1] + [m] * 9, m, 9)["class"] == "LOCK")
    counter, _ = integer_orbit(9, F(1, 8), 90, frozen)
    check("outside T6: alpha=1/8 labels LOCK despite two branches",
          counter["class"] == "LOCK" and counter["post_set"] == [7, 8] and counter["R_m"] == 10)
    slow, _ = integer_orbit(10, F(1, 1000), 100, frozen)
    check("outside T6: slow alpha yields OTHER before entry", slow["class"] == "OTHER" and slow["t_star"] is None)
    _, choices = integer_orbit(50, F(1, 4), 3, frozen)
    check("s0/s1/s2 indexing: K=50 gives 25/32/41", choices == D.first_selections(50) == [25, 32, 41])
    for K in (6, 7, 8, 9, 10, 11, 20, 50, 100):
        check(f"K={K}: first-update candidate formula vs every bin", D.formula_argmax(K) == D.argmax_D1(K))


if __name__ == "__main__":
    main()
