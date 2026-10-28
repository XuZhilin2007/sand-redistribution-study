# Sand v0.5.0 release and reproduction

v0.5.0 publishes the owner-approved `v0.5.0-rc.1` candidate without changing
its scientific claims, proof, algorithms, original evidence or frozen bytes.
The complete approved 114-file package is committed at
[`releases/v0.5.0/`](releases/v0.5.0/). Its RC identifier/status fields are
retained as historical provenance; the new
[release binding](release_records/v0.5.0/RELEASE_BINDING.json) records the
authorization and exact path mapping for this publication.

The claimed canonical-XO universal effective set is exactly the closed
interval with endpoints approximately 0.24615482985838558 and
0.29164543617016553, determined by K=10 and K=6. The
[approved proof](releases/v0.5.0/docs/alpha_phase_2026-10-27/PROOF.md),
[October 28 audit](releases/v0.5.0/docs/alpha_phase_audit_2026-10-28/AUDIT_CN.md),
and [scope ledger](releases/v0.5.0/V050_CLAIMS.json) state the exact theorem,
qualifications, and stage responsibilities. This release does not claim a
complete fixed-K phase diagram or independent human/formal verification.

From the Git checkout, use a new output directory outside the frozen package:

```text
python -B -S releases/v0.5.0/tools/reproduce_v050.py --out-dir ../sand-v050-run-001 --expected-manifest-sha256 5dcbc85e889f3c54cc1577dbe68c1e5f0c63cdb586fc8739499890c046ad76b0
```

From the extracted verification ZIP, use `tools/reproduce_v050.py` instead.
The standard-library command redirects all generated evidence, refuses
existing or in-package output paths, verifies the original evidence before
and after execution, and compares regenerated mathematical content against
the audited references. Use the new entry point rather than the preserved
historical alpha `verify_all.py`, whose default reports write in place.

The immutable approved source ZIP is published with its original filename:
`sand-v0.5.0-rc.1-verification-source.zip`, 2,612,306 bytes, SHA-256
`3e49c412a5df4982ad14035b4d870e44a2eb489ec4349af2c3f0fa697de83057`.
The [original candidate receipt](release_records/v0.5.0/CANDIDATE_RECEIPT.json)
and [frozen manifest](releases/v0.5.0/V050_MANIFEST.json) retain their approved
bytes. The archive is uploaded directly, without repacking or renaming it.

The release-tree verification and publication records are under
[`release_records/v0.5.0/`](release_records/v0.5.0/). Any copied candidate
readiness document there is a byte-preserved historical snapshot; its
original relative links refer to the preparation record context. This file
is the current Git-tree reproduction guide.

All pre-existing public v0.4.0 tree files and Git history are preserved.
The root README and original root verification entry remain the v0.4.0
baseline deliberately; v0.5.0 is self-contained in its versioned directory.
Its published-v0.4.0 ZIP and private 37-file research freeze remain distinct
references inside that directory. Successful reproduction is execution
evidence for the audited computer-assisted proof, not another independent
proof or review.
