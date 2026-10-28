# Sand v0.5.0-rc.1 — release-readiness assessment

Prepared 28 October 2026, Asia/Shanghai. **Recommendation: the exact candidate
is ready for owner review and release approval. The three identified
reproduction/provenance/scope blockers are resolved.** This recommendation
does not authorize a Git tag or public release. Neither was created.

This work consolidates the already-audited alpha-phase theorem. No new
mathematical research, model change, proof amendment or scientific output
change was introduced. The October 28 mathematical audit remains an
AI-assisted audit, not independent human peer review or formalization.

**Candidate identity and contents.**

| Item | Exact candidate |
|---|---|
| Identifier/status | `v0.5.0-rc.1`, unpublished release candidate |
| Source ZIP | [sand-v0.5.0-rc.1-verification-source.zip](sand-v0.5.0-rc.1-verification-source.zip) |
| ZIP size | 2,612,306 bytes |
| ZIP SHA-256 | `3e49c412a5df4982ad14035b4d870e44a2eb489ec4349af2c3f0fa697de83057` |
| Frozen inventory | [package/V050_MANIFEST.json](package/V050_MANIFEST.json) |
| Manifest SHA-256 | `5dcbc85e889f3c54cc1577dbe68c1e5f0c63cdb586fc8739499890c046ad76b0` |
| Complete extracted package | 114 files, 5,337,935 bytes; 113 entries plus the manifest itself |
| Detached identity receipt | [CANDIDATE_RECEIPT.json](CANDIDATE_RECEIPT.json) |
| Authoritative current scope ledger | [package/V050_CLAIMS.json](package/V050_CLAIMS.json) |

The inventory enumerates every actual package file, including the new README
and engineering tools; both original alpha proof documents and both figures;
all original alpha sources, evidence and exploratory records; the October 28
audit and specifically linked evidence; earlier reviews and corrections;
the four mathematical library/kernel dependencies; and licensing material.
Each entry records SHA-256, size, role, original path/hash and whether its
workspace origin was tracked, modified or untracked.

The byte-preserved published v0.4.0 verification ZIP and publication receipt
are under `references/published-v0.4.0/`. The separate original 37-file private
research freeze and manifest are under `references/research-v0.4.0/`. Its
original README is preserved there and is not confused with the new RC
README. The complete public baseline has its own 78-file manifest.

The ZIP is a selected verification bundle, not the entire private research
repository. Unrelated private work, Git metadata and the large duplicate
source tree of the October 28 isolation run are omitted. The audit's linked
reproduction reports are included. The full original repository and audit
tree remain unchanged in the workspace.

**Disposition of the audit findings.**

| Finding | Resolution | Evidence |
|---|---|---|
| Original aggregate overwrites existing evidence, and its `--report` does not redirect child stages | New `tools/reproduce_v050.py` requires a new external output directory, supplies every child `--report`, supplies the regenerated finite witness via `--input`, and rejects reuse or in-package paths | 29 [guard controls](GUARD_TESTS.json); isolated full run; originals unchanged |
| Hashes generated after reproduction do not establish preservation of old evidence; new proof documents were not frozen | Pre-existing complete manifest is checked before any output/process and after the run; old report/source/witness bindings and pre-frozen mathematical payload digests are checked independently of new hashes | Manifest; [preflight](isolated_run/preflight.json); [isolated verification](ISOLATED_VERIFICATION.json) |
| Exact-stage scope wrongly credits v0.4.0 with the wider universal inclusion | New frozen scope ledger attributes the new inclusion to `finite_cover` and `uniform_probe`; only new consolidated reports receive the corrected scope and per-stage responsibilities | [metadata correction record](isolated_run/metadata_corrections.json); [new exact report](isolated_run/reports/exact_results.json) |
| Independence could be overstated | Guides, ledger and aggregate distinguish second sign algorithms, shared trajectory algebra, analytic dependencies, audit-program replay and independent human/formal proof | [provenance guide](package/docs/V050_PROVENANCE.md); aggregate qualification |

Original `exact_results.json` and the original producer's newly regenerated
raw output intentionally retain the old scope sentence as historical/raw
evidence. They are not the authoritative new release scope. The correction
record states its old and new values; restoring that value and removing the
new `proof_scope` metadata reproduces the raw report exactly. No other output
field is changed during consolidation.

**Reproduction procedure.** Verify the ZIP against the detached receipt,
extract it to a fresh location, then run from its root:

```text
python -B -S tools/reproduce_v050.py --out-dir ../sand-v050-run-001 --expected-manifest-sha256 5dcbc85e889f3c54cc1577dbe68c1e5f0c63cdb586fc8739499890c046ad76b0
```

Use a different, nonexistent output directory for each run. The historical
alpha `verify_all.py` and its in-place command examples remain preserved
reference material; they are not the v0.5.0 reproduction entry point.
The new entry point cannot overwrite their saved reports.

The five mathematical reports are compared in full against the pre-frozen
audited originals except for top-level `elapsed_seconds` and the finite
auditor's dependent witness digest. That digest is separately required to
match the actual regenerated witness; the original witness binding is
checked too. Scientific fields, source hashes, exact fractions, obligations,
cells, tail starts, terminal-box digest, exact bounds and claims are not
excluded. An unexpected mathematical discrepancy stops reproduction.

**Actual isolated verification.** The exact ZIP identified above was extracted
and verified as a 114-file package before execution. Python 3.14.6 ran with
`-B -S`; child processes had an empty PATH, no site packages, no package Git
metadata, and used no network. The run completed in 258.487 seconds:

| Stage | Exit | Actual result |
|---|---|---|
| `verify_phase` | 0 | Boundary polynomials, root brackets, exact exclusions and fixed-K results match the audited original |
| `audit_phase` | 0 | 60 formal identities, 95 exact moment cases and 40 integer cases match |
| `finite_cover` | 0 | 196 ordinary K, 1,045 cells and the additional K=10 upper cell match |
| `audit_cover` | 0 | 1,046 cells and four mathematical rejection controls match |
| `uniform_probe` | 0 | 20,554 visited / 11,813 terminal boxes; all exact bounds and the terminal-box digest match |
| Published v0.4.0 verification | 0 | 6/6 checks passed; all 78 public frozen files unchanged |
| Copied October 28 audit program | 0 | Entire regenerated audit JSON equals the saved audit result |

The [full new run record](isolated_run/verification.json), raw reports,
consolidated reports, logs and correction ledger are retained separately
from original evidence. A separate [summary](ISOLATED_VERIFICATION.json)
pins the candidate identity and records the checks. Repacking the extracted
package with its included builder, without Git on PATH, produced a
**byte-identical source ZIP**, recorded in [REPACK_RECEIPT.json](REPACK_RECEIPT.json).
That determinism was tested in this Python/zlib environment; it is not a
claim about all compression-library versions.

**Preservation and actual diff.** All 1,676 pre-existing files in the
preparation inventory retained their bytes, including all 699 files pinned
by the October 28 audit, historical reviews/corrections, the alpha proof and
evidence, tracked research files and v0.4.0 release records. Git HEAD remains
`2c337bb4a52a7e79430c6430413893dd79c6e9e4`. The two pre-existing maintained-doc
modifications remain unchanged and are not this task's edits.

The engineering additions are four new Python modules and two new release
guides. New generated materials are the explicit candidate scope/inventory,
ZIP/receipts and separate reproduction/readiness evidence. All 108 files
copied from workspace sources into the bundle match their respective source
bytes, with new versus historical roles distinguished in the inventory.
There are zero modifications to original scientific algorithms, proofs or
reports. The [actual additions patch](CONSOLIDATION_DIFF.patch),
[diff review](DIFF_REVIEW.json) and [workspace preservation check](WORKSPACE_AFTER.json)
make this distinction inspectable.

**Remaining gates and recommendation.** There is no unresolved blocker among
the identified safe-reproduction, complete frozen provenance and scope
metadata issues. The exact archive is ready for review and explicit approval.
Successful execution reproduces the audited computation and its documented
analytic assembly; it is not an independent proof of the theorem.

The current HEAD is the older baseline and does not contain the untracked
alpha theorem or this RC. **Do not tag that HEAD as v0.5.0.** After the owner
reviews and approves the exact candidate, the release commit/tree must contain
the approved bytes and be checked against this inventory before any tag or
publication. If any reviewed byte changes, the current hashes and isolated
acceptance no longer identify that new candidate: prepare a new candidate
identity and repeat the relevant freeze/verification steps.

No tag, public release, upload, or Git-history mutation has been performed.
