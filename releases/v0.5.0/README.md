# Sand v0.5.0 release candidate

Candidate identifier: **v0.5.0-rc.1**. Prepared 28 October 2026, Asia/Shanghai.
This is an unpublished review package. No Git tag or public release has been
created or approved by this package.

The unchanged canonical XO theorem is

\[
 \bigcap_{K\ge6}W_K=W_6\cap W_{10}=[\alpha_-,\alpha_+],
 \quad\alpha_-=0.24615482985838558\ldots,
 \quad\alpha_+=0.29164543617016553\ldots.
\]

The exact polynomials, endpoint cases, fixed-K results and analytic
dependencies remain in the original [proof](docs/alpha_phase_2026-10-27/PROOF.md)
and [assessment](docs/alpha_phase_2026-10-27/RESULTS.md). The
[October 28 audit](docs/alpha_phase_audit_2026-10-28/AUDIT_CN.md) found no
mathematical blocker and identified the reproduction/provenance issues this
candidate addresses. The complete fixed-K phase diagram remains outside
the claim.

**Use the new reproduction entry point.** The scripts, reports and command
examples in the historical alpha-phase directory are preserved verbatim.
Their old `verify_all.py` command writes in place. Do not execute that
historical entry point against preserved evidence. Use this command from
the extracted candidate directory, choosing a directory that does not exist:

```text
python -B -S tools/reproduce_v050.py --out-dir ../sand-v050-run-001
```

For verification against the externally recorded candidate identity, also
pass `--expected-manifest-sha256` with the value from `CANDIDATE_RECEIPT.json`
distributed alongside the ZIP. The receipt is external to avoid a self-hash
cycle. Verify the ZIP's SHA-256 against that receipt before extraction.
The package verifies the manifest's entire file list and every file hash
before creating output or launching a stage. The externally supplied digest
pins which manifest is expected; without that argument the check establishes
internal consistency, not authenticity of a substituted package.

The command requires a **new output directory outside the candidate**.
Existing output paths, paths inside the candidate, and preflight hash or
scope failures are rejected. Every stage is given an explicit report path;
the finite auditor is given the newly generated finite witness as its input.
All original reports are checked before and after reproduction. The frozen
package itself receives no output files.

Python's standard library is sufficient. This candidate was tested with
Python 3.14.6 on Windows, using `-B -S`, an empty child-process PATH, no Git
metadata and no network access. The existing sources use modern Python;
other interpreter versions have not been validated here. No plotting package
is required to reproduce the mathematics. PNG/SVG figures and the historical
optional matplotlib script are included as reference artifacts.

The run performs the five original alpha proof stages, verifies the
published v0.4.0 verification ZIP separately, and replays an unchanged copy
of the October 28 audit reconstruction. This replay is additional
reproduction evidence, not another independent mathematical review.

Output layout:

| Location | Meaning |
|---|---|
| `preflight.json` | Candidate identity and original-evidence checks performed before execution |
| `raw/*.json` | Exact original producer outputs, newly generated into this run |
| `reports/*.json` | Consolidated reports with explicit stage scope; the stale exact-stage scope is corrected here only |
| `metadata_corrections.json` | The single scope correction, old/new values, and a check that no scientific field changed |
| `baseline_extract/` | Separate extraction of the byte-preserved published v0.4.0 ZIP |
| `independent_checks.py/.json` | Unchanged copied audit program and its matching regenerated output |
| `logs/` | Actual stdout/stderr of every stage |
| `verification.json` | Complete reproduction outcome, stage exits, hashes, scope and limitations |

`raw/exact_results.json` deliberately retains its original producer's stale
scope string so the reference implementation remains byte-identical.
It is labelled raw historical-producer output. The authoritative new scope
is in `V050_CLAIMS.json` and `reports/exact_results.json`: the exact stage
provides low-dimensional results and exclusions; the new universal inner
inclusion comes from `finite_cover` and `uniform_probe`. The v0.4.0 baseline
supports only its original narrow alpha interval.

Regenerated mathematical content must equal the pre-frozen audited
reference content. Only top-level elapsed time, and the finite auditor's
witness hash after verifying it binds the actual regenerated witness, are
excluded from that semantic comparison. Source hashes, exact fractions,
root brackets, scopes in raw evidence, obligations, cells, coverage,
terminal-box digest, bounds and claims are compared. Any scientific
discrepancy stops the run; it cannot be accepted as a metadata correction.

The published baseline is kept as its original ZIP and publication receipt.
It is reverified in its own extracted tree, with its own 78-file public
manifest. The 37-file private research freeze is separately preserved under
`references/research-v0.4.0/`; it is not substituted for the public baseline
or silently adapted to the new package README.

The [provenance guide](docs/V050_PROVENANCE.md) explains file roles,
independence limits and content identity. `V050_MANIFEST.json` is the frozen
inventory of the actual package, including this README, both proof documents,
the figures, original evidence, audit records, new tooling and scope ledger.
It rejects missing, modified or additional files. The manifest itself is
the sole file excluded from its own inventory; the external receipt pins it
and the complete ZIP.

**Interpretation of success.** A successful run reproduces the audited
calculation and its documented analytic assembly in this exact candidate.
It does not independently prove the theorem, formalize Theorem A, establish
novelty, imply human peer review, or authorize publication. Review and explicit
owner approval remain necessary before tagging or publishing v0.5.0.
