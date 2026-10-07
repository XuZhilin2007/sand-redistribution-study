"""Validate the XO supplements without an extension-lab checkout.

Run at the repository root, for example:
  python -B experiments/xo_extensions/validate.py --scope full --out-dir ../xo-v020-run

Only a new or empty output directory is accepted. The frozen S2 source is
loaded from this repository, not from another checkout. The checked configs
describe the fixed verifiers and independently configure the full-state run.
"""
from __future__ import annotations

import argparse
import contextlib
from fractions import Fraction as F
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import platform
import sys
import time

from exact_engine import ROOT, FROZEN_S2, frozen_classifier, integer_orbit, proof_identities


def load(relative):
    path = ROOT / relative
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_config(relative):
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def save(directory, name, value):
    (directory / name).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n",
                                encoding="utf-8")


def checked_configs():
    alpha = read_config("experiments/xo_alpha_extension/config.json")
    assert alpha["K_values"] == list(range(6, 54)) + [100, 200, 400]
    assert alpha["alpha_values"] == ["1/8", "1/4", "3/8", "1/2", "3/4"]
    assert alpha["K6_extra_alpha"] == ["1/1000", "1/100", "4/27", "6/31", "1/3", "99/100", "999/1000"]
    assert alpha["requested_H_multiplier"] == 10 and alpha["large_K_full_state_cap"] == 300
    assert alpha["float_proxy_steps"] == 200000 and alpha["K6_extra_H"] == 400
    t6 = read_config("experiments/xo_t6_classifier/config.json")
    assert t6["K_values"] == list(range(6, 301)) and t6["alpha"] == "1/4" and t6["H_multiplier"] == 10
    aux = read_config("experiments/xo_supporting_results/config.json")
    assert aux["F"] == {"seed": 20261006, "canonical_samples": 40, "perturbed_samples": 20,
                        "sample_H": 150, "xo_sanity_K": 12, "xo_sanity_H": 2000,
                        "alpha": "1/4", "tolerance": "1/1000000000000"}
    assert aux["D"] == {"K_min": 6, "K_max": 4000, "second_step_K_max": 60, "alpha": "1/4"}
    assert aux["S12"]["K_values"] == [50, 200, 8, 14, 26]
    assert aux["S12"]["epsilon_values"] == ["1/100", "1/20", "1/10", "1/5"]
    assert aux["S12"]["steps"] == 200000 and aux["S12"]["arithmetic"] == "scalar float proxy after exact full-state setup"
    assert aux["S12"]["seed_rule"] == "1000*K + int(100*epsilon)"
    assert aux["S6"] == {"K_values": [20, 50, 100], "H": 1000, "seed": 42, "alpha": "1/4",
                         "kinds": ["canonical", "uniform", "near_uniform", "suffix_transfer_stop_control"]}
    assert aux["C1"] == {"R": 12, "T": 5, "alpha": "1/4", "seed": 20261006, "max_active_steps": 5,
                         "tolerance": "1/1000000000000"}
    config = read_config("experiments/xo_extensions/config.json")
    assert config["schema_version"] == 1 and config["model"] == "canonical XO"
    return {"A": alpha, "T6": t6, "auxiliary": aux, "independent": config}


def source_hashes():
    paths = {FROZEN_S2}
    for name in ("xo_alpha_extension", "xo_t6_classifier", "xo_supporting_results", "xo_extensions"):
        directory = ROOT / "experiments" / name
        paths.update(p.relative_to(ROOT).as_posix() for p in directory.glob("*.py"))
        paths.add(f"experiments/{name}/config.json")
    paths.update(p.relative_to(ROOT).as_posix() for p in (ROOT / "docs").glob("XO_*.md")
                 if p.name != "XO_EXTENSION_VALIDATION.md")
    paths.add("tests/test_xo_extensions.py")
    return {rel: hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() for rel in sorted(paths)}


def run_verifier(label, relative, output):
    print("VERIFIER START", label, flush=True)
    log = io.StringIO()
    try:
        with contextlib.redirect_stdout(log):
            value = load(relative).main()
    finally:
        (output / (label + ".log")).write_text(log.getvalue(), encoding="utf-8")
    assert value["status"] == "PASS", label
    save(output, label + ".json", value)
    print("VERIFIER PASS", label, flush=True)
    return value


def independent_checks(config, scalar, probes, descent, frozen):
    trajectories, stress = [], []
    assert config["alpha"] == "1/4" and config["H_multiplier"] == 10
    for K in config["full_mass_vector_K"]:
        print("INTEGER FULL-STATE START", K, flush=True)
        row, choices = integer_orbit(K, F(config["alpha"]), config["H_multiplier"] * K, frozen)
        other = scalar.scalar_classifier(K)
        assert choices == other["s_list"]
        assert row["class"] == other["class"] == ("LOCK" if K % 6 < 3 else "CYCLE")
        assert row["R_m"] == other["R_m"]
        trajectories.append(row)
    for value in config["K6_alpha_boundary_cases"]:
        row, _ = integer_orbit(6, F(value), config["K6_H"], frozen)
        assert row["t_star"] is not None and row["confinement_time"] is not None
        stress.append(row)
    counter, _ = integer_orbit(9, F(1, 8), 90, frozen)
    assert counter["class"] == "LOCK" and counter["R_m"] == 10 and counter["post_set"] == [7, 8]
    slow, _ = integer_orbit(10, F(1, 1000), 100, frozen)
    assert slow["class"] == "OTHER" and slow["t_star"] is None
    for K in config["first_selection_K"]:
        _, choices = integer_orbit(K, F(1, 4), 3, frozen)
        assert choices == descent.first_selections(K)
    forced = []
    pattern = config["forced_lower_pattern"] * config["forced_pattern_repetitions"]
    for K in config["forced_full_state_K"]:
        state = probes.full_state_to_upper(K)
        p = state["p"]
        for force in pattern:
            discrepancies = probes.discrepancies(p)
            assert max(discrepancies) == state["a"]
            assert discrepancies.index(max(discrepancies)) + 1 == state["m"]
            s = state["m"] - 1 if force else state["m"]
            old = [F(0)]
            for mass in p:
                old.append(old[-1] + mass)
            moved = old[s] / 4
            G = [F(0)]
            for mass in state["g"]:
                G.append(G[-1] + mass)
            new = [F(3, 4) * prefix + moved * G[j] if j <= s
                   else prefix - moved * (1 - G[j]) for j, prefix in enumerate(old)]
            p_new = [new[j] - new[j - 1] for j in range(1, K + 1)]
            assert p_new == probes.update(p, s, state["g"])
            p = p_new
        forced.append({"K": K, "exact_observed_states": len(pattern), "status": "PASS"})
    return {"status": "PASS", "engine": "integer mass vector, shared denominator, no scalar advance",
            "full_trajectories": trajectories, "K6_alpha_boundary_cases": stress,
            "finite_window_counterexample": counter, "slow_alpha_other": slow,
            "forced_rule_exact_checks": forced, "algebra_checks": proof_identities()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scope", choices=("core", "auxiliary", "full"), default="full")
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    if not __debug__:
        parser.error("Assertions are required; do not run with -O.")
    output = args.out_dir.resolve()
    if output == ROOT or output.is_relative_to(ROOT / ".git"):
        parser.error("Choose a separate validation output directory.")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error("Output directory must be new or empty; archived evidence is never overwritten.")
    output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    configs = checked_configs()
    hashes = source_hashes()
    results = {}
    if args.scope in ("core", "full"):
        results["A"] = run_verifier("A", "experiments/xo_alpha_extension/verify_alpha.py", output)
        assert len(results["A"]["matrix"]) == 255 and len(results["A"]["K6_extra"]) == 7
        for row in results["A"]["matrix"]:
            assert row["requested_H"] == 10 * row["K"]
            assert row["verified_H"] == (300 if row["K"] > 60 else 10 * row["K"])
        assert sum(row["steps"] for row in results["A"]["finite_float"]) == 51000000
        results["T6"] = run_verifier("T6", "experiments/xo_t6_classifier/verify_t6.py", output)
        assert len(results["T6"]["classifications"]) == 295 and len(results["T6"]["closed_form_checks"]) == 15
        frozen = frozen_classifier()
        scalar = load("experiments/xo_t6_classifier/verify_t6.py")
        frozen_rows = []
        for K in configs["T6"]["K_values"]:
            row = scalar.scalar_classifier(K)
            value = frozen(row["s_list"], (5 * K + 5) // 6, K)
            assert value["class"] == row["class"] and value["R_m"] == row["R_m"] and value["t_star"] == row["t_star"]
            frozen_rows.append({"K": K, "H": row["H"], **value})
        save(output, "frozen_classifier.json", {"status": "PASS", "source": FROZEN_S2, "rows": frozen_rows})
        independent = independent_checks(configs["independent"], scalar,
                       load("experiments/xo_supporting_results/probes.py"),
                       load("experiments/xo_supporting_results/mine_descent.py"), frozen)
        save(output, "independent.json", independent)
        results["independent"] = {"full_trajectories": len(independent["full_trajectories"]),
                                 "K6_alpha_boundary_cases": len(independent["K6_alpha_boundary_cases"]),
                                 "frozen_classifier_cases": len(frozen_rows)}
    if args.scope in ("auxiliary", "full"):
        for label, name in (("F", "verify_flattail"), ("D", "mine_descent"), ("probes", "probes"), ("C1", "check_c1")):
            results[label] = run_verifier(label, "experiments/xo_supporting_results/" + name + ".py", output)
        assert results["F"]["statistics"]["canonical_runs"] == 40
        assert results["F"]["statistics"]["perturbed_runs"] == 20
        assert results["D"]["first_step_checked"] == 3995 and len(results["D"]["second_step_records"]) == 55
        assert len(results["probes"]["restricted_noise"]) == 20 and len(results["probes"]["finite_initial_states"]) == 12
        assert all(row["arithmetic"] == configs["auxiliary"]["S12"]["arithmetic"] for row in results["probes"]["restricted_noise"])
        assert results["C1"]["verified_active_steps"] == 5
    assert source_hashes() == hashes, "Sources changed during validation"
    report = {"status": "PASS", "version": "0.2.0", "scope": args.scope, "python": platform.python_version(),
              "arithmetic": {"A": "Fraction full vectors + separately labelled float scalar proxies",
                             "T6": "Fraction full-state entry + exact scalar closure",
                             "independent": "integer full mass vectors with common denominator",
                             "S12": "scalar float proxy after exact full-state setup",
                             "F_D_S6_C1": "exact rational/integer"},
              "source_sha256": hashes, "frozen_configs": configs, "results": results,
              "elapsed_seconds": round(time.perf_counter() - started, 3),
              "scope_note": "Finite checks support the deposited proofs; no periodicity, arbitrary-alpha finite-window guarantee, or novelty claim."}
    save(output, "manifest.json", report)
    print("VALIDATION PASS", args.scope, "elapsed", report["elapsed_seconds"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
