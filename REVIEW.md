# Release Review

Tracks the review of the received ICM v281 reproduction package. The package is
ready for release when every box below is ticked. Findings are numbered (F1,
F2, …) and referenced from the checklist; record the decision and commit when
closing one.

**Per-file workflow.** Confirm the file matches `PACKAGE_SHA256.json`, commit it
as received, commit any `ruff format` changes on their own, then commit each
logical change separately. Do not change numerical behaviour (thread counts,
seeds, pinned versions) without explicit sign-off.

## 1. Environment

- [x] Python environment managed with uv (`pyproject.toml`, `uv.lock`, `.python-version`)
- [x] Skip dotfiles (`.DS_Store`) when walking pipeline folders
- [x] Show and record run paths relative to the repository root
- [x] ruff configuration (PEP 8, Google docstrings, received code excluded)
- [ ] Pin R exactly: `R = "= 4.4.2"` in `rproj.toml`, then `rig proj lock`
- [ ] Commit `.gitignore` rig block and `.rvenvlib/`
- [ ] Commit `README.md` (Python and R sections)
- [ ] Drop `black` from the dev dependency group (`ruff format` replaces it)
- [ ] Make the pipeline use the pinned R rather than the first `Rscript` on `PATH` (F6)

## 2. Pipeline launchers (`pipeline/*.py`)

Status per file: received → formatted → reviewed → findings resolved.

- [ ] `run_all.py` — reviewed; resolve F1–F3, F7, F8
- [ ] `run_stage.py` — reviewed; resolve F4–F6, F8
- [ ] `materialize.py` — resolve F13
- [ ] `run_integrity.py` — resolve F10, F11
- [ ] `stage_guard.py`
- [ ] `prepare_primary_inputs.py`
- [ ] `prepare_display_tables.py`
- [ ] `prepare_render_inputs.py`
- [ ] `run_figures.py` — includes F17
- [ ] `document_provenance.py`
- [ ] `verify_inputs.py`
- [ ] `verify_code.py`
- [ ] `validate_workflow.py`
- [ ] Configuration: `stages.json`, `root_map.json`, `input_manifest.json`,
      `display_input_map.json`, `render_jobs.json`, `render_requests.json`,
      `workspace_directories.json`, `environment.json`
- [ ] Decide fate of legacy `pipeline/requirements.txt` and `requirements-lock.txt`

## 3. Notebook layer (`code/`)

- [ ] `extend_analysis.py` — source patching first (F9)
- [ ] `reproduce.py` — includes F15
- [ ] `check_primary_predictions.py`
- [ ] `check_expression_and_text.py`
- [ ] `check_correspondence.py`
- [ ] `check_genomewide.py`
- [ ] `check_integrity_guard.py`
- [ ] `audit_result_annotations.py`
- [ ] `collect_figures.py`
- [ ] `collect_source_tables.py`
- [ ] `rebuild_supplement.py`
- [ ] `supplement_draw_job.py`

## 4. Entry points

- [ ] `execute_notebook.py`
- [ ] `verify_package.py`
- [ ] `Reproduce_ICM_v281.ipynb` — update cell 474fd814 (still cites
      `requirements-lock.txt`); review the resume heuristic (F15)

## 5. Numerical code (`pipeline/code`, `extensions`, `renderers`)

About 88,500 lines in 504 scripts; 47 stages each name one entry script. Review
by evidence rather than line by line.

- [ ] Map each stage's entry script and everything it imports or calls
- [ ] Identify scripts never executed; confirm against
      `validation/EXCLUDED_UNUSED_FRAGMENTS.json`; decide whether to ship them
- [ ] Trace each main-text figure and the molecular table to its producers
      (`FIGURE_MAP.csv`, `stages.json`) and review those producers in full
- [ ] Review `pipeline/extensions/` (16 scripts)
- [ ] Review `pipeline/renderers/` (24 scripts)
- [ ] Replace hardcoded `/tmp` and `/private/tmp` `MPLCONFIGDIR` values (F16)
- [ ] Remove source-text patching dependencies, then format and lift the ruff exclusion

## 6. Documentation

- [ ] Merge `README.txt` into `README.md`; remove `README.txt`
- [ ] `ENVIRONMENT.json` — claims "pip check passed"; align with uv
- [ ] Check `COMPREHENSIVE_AUDIT.txt`, `RESULT_HARDCODING_AUDIT.txt`,
      `PREPROCESSING.txt` against the code as released
- [ ] Check `DATA_PROVENANCE.csv`, `DATA_LICENSES.txt`, `REFERENCES.txt`
- [ ] Decide what of `validation/` and `reference/` ships

## 7. Release

- [ ] Decide how `data/` is distributed (currently ignored by git)
- [ ] Ignore generated folders: `pipeline/rebuild*/`, `pipeline/cache/`,
      `pipeline/audit/`, `pipeline/reports/`, `pipeline/figures/`,
      `pipeline/project/`, `pipeline/input_aliases.json`
- [ ] Delete `pipeline/rebuild.bak/` and editor backups (`*~`)
- [ ] Replace `PACKAGE_SHA256.json` and `pipeline/SHA256_CODE.txt` with checksums
      of the reviewed release
- [ ] Fresh full run from a clean clone; all checks pass
- [ ] Notebook outputs contain no machine-specific paths (`grep /Users/`)
- [ ] Tag the release

## Findings

| ID | Where | Finding | Proposed action | Status |
| -- | ----- | ------- | --------------- | ------ |
| F1 | `run_all.py` | `pipeline/reports/` and `pipeline/audit/` are created only by `code/reproduce.py:59`; running `run_all.py` directly, as its docstring documents, fails after all 47 stages | Create the folders in `run_all.py` or in each writer | Open |
| F2 | `run_all.py:77` | `verify_inputs.py` runs after all analyses, too late to prevent a bad run (`materialize.py` already checks input hashes) | Move before `materialize`, or drop as redundant | Open |
| F3 | `run_all.py:52-69` | Hand-coded file copies and `prepare_display_tables.py` before two stages are dependencies missing from `stages.json` | Document; consider moving into `stages.json` | Open |
| F4 | `run_stage.py:41-44` | Stage output goes only to its log; a failure surfaces as a bare `CalledProcessError` | Print the log tail on failure | Open |
| F5 | launchers | `assert` used for runtime checks (about 32), skipped under `python -O` | Raise explicit exceptions | Open |
| F6 | `run_stage.py:35` | `Rscript` taken from `PATH`; R version not enforced | Check `R.version` or resolve R via rig | Open |
| F7 | `run_all.py:10` | Positional `sys.argv` parsing; `--resume` before the path is misread as the directory | Use `argparse` | Open |
| F8 | launchers | Single-letter globals, multi-imports, mid-file import, no `main()`, missing docstrings, long lines | Restyle | Open |
| F9 | `code/extend_analysis.py:101` | Patch `.replace("W=P.parents[2]", ...)` never matches: the target already reads `W = P.parents[1]` (formatted, see `validation/black_pipeline.log`) | Remove the dead patch | Open |
| F10 | `run_integrity.py` | Receipt does not hash `code/`, so edits there pass a resume check | Add `code/` to the code snapshot | Open |
| F11 | `run_integrity.py` | Receipt recorded twice (`run_all.py` and `code/reproduce.py`); the second overwrites the first | Record once | Open |
| F12 | `run_integrity.py:77` | Fingerprints every file in `data/source/atlas/` by glob rather than the declared inputs | Fingerprint the named CSVs only | Open |
| F13 | `materialize.py` | Unused `import shutil`; `assert` for hash checks (see F5) | Clean up | Open |
| F14 | many scripts | Thread variables set in about 89 places | Keep (single-threaded by design); consider one shared helper later | Open |
| F15 | `code/reproduce.py`, notebook | Resume chosen by `(REBUILD / "status").is_dir()`, so an interrupted run is resumed and then rejected | Base resume on the receipt's existence | Open |
| F16 | `pipeline/code`, `extensions`, `renderers` | Hardcoded `/tmp` and `/private/tmp` `MPLCONFIGDIR` | Point at the run's cache folder | Open |
| F17 | `run_figures.py` | Author-only `--discover` branch hardcodes another user's absolute paths | Remove before release | Open |
