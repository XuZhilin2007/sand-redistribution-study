Publishes the approved `v0.5.0-rc.1` verification package byte-for-byte.

- Canonical XO universal effective set: `[alpha_-, alpha_+]`, approximately `[0.24615482985838558, 0.29164543617016553]`, with K=10 and K=6 determining the boundaries. Exact polynomials, endpoint cases and proof are included.
- Includes the frozen proof/evidence, October 28 audit, corrected stage-scope ledger, safe reproduction tooling and release provenance. Existing v0.4.0 files and history are preserved.
- Reproduction checks the complete freeze before/after execution, writes only to a new external directory, and compares mathematical outputs with the audited references. This is reproducibility evidence, not new independent human or formal proof.

The Git-tree entry point is `releases/v0.5.0/tools/reproduce_v050.py`; see [the release/reproduction guide](https://github.com/XuZhilin2007/sand-redistribution-study/blob/v0.5.0/V050_RELEASE.md). The original root README/verifier remain the preserved v0.4.0 baseline.

The source asset deliberately retains its approved filename `sand-v0.5.0-rc.1-verification-source.zip` and immutable RC metadata. Its SHA-256 is `3e49c412a5df4982ad14035b4d870e44a2eb489ec4349af2c3f0fa697de83057` (2,612,306 bytes). It was uploaded without rebuilding. The attached receipt and manifest identify the exact approved content.

Full fixed-K phase diagrams and a general finite membership horizon remain outside this release's claims.
