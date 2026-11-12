# Comparing a rerun with this baseline

`RUN_INTEGRITY.verified.json` holds the hash of every numerical file in the
verified run, and `results.sha256` holds the hashes of the published results.
A rerun should match both, except for the files below, which record something
that changes from run to run. A mismatch in any other file is a real change.

Compare a run made in the same folder name as the baseline (`pipeline/rebuild`);
several outputs record the absolute path of the run directory.

| Files | Why they differ | Count |
| --- | --- | --- |
| `status/*.json` | each records the stage's run time in `seconds` | 63 |
| `legacy/results_v236/protocol/protocol_receipt.json` | records the UTC time of the run | 1 |
| `legacy/results_v237/protocol/lock_receipt.json` | records the UTC time of the run | 1 |
| `current/outputs/ICM_groundlike_preliminary/cell_scores.csv` | the baseline file has its columns in a random order (REVIEW.md F39); the code now sorts them, so the first rerun differs once. Compare it by column name | 1 |

For the first three rows, compare every field except the time. For
`cell_scores.csv`, sort the columns of both files by name and compare the
values. A rerun made after the F39 fix produces a stable file; replace the
baseline hash when the baseline is next refreshed.
