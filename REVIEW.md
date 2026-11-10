# Release Review

Tracks the review of the received ICM v281 reproduction package. The package is
ready for release when every box below is ticked. Findings are numbered (F1,
F2, …) and referenced from the checklist; record the decision and commit when
closing one.

**Order.** Review in execution order: sections 3–6 follow the call chain of a
notebook run, so each file is reviewed after everything that produces its
inputs. Jump ahead only on evidence: a finding that affects the file in hand, a
failed or changed full run, or a release decision that removes files from
scope. Section 2 comes first. The independent checks (section 7) recompute
from source data and can be reviewed at any time.

**Per-file workflow.** Confirm the file matches `PACKAGE_SHA256.json` and commit
it as received; review it and record findings; commit any `ruff format` changes
on their own; then commit each approved fix separately, and test restructures
against the previous version. Do not change numerical behaviour (thread counts,
seeds, pinned versions) without explicit sign-off.

**Cross-file couplings.** Before changing a file, check the couplings table
below for anything that patches, hashes or names it, and update the table when
a change adds or removes a coupling.

## 1. Environment

- [x] Python environment managed with uv (`pyproject.toml`, `uv.lock`, `.python-version`)
- [x] Skip dotfiles (`.DS_Store`) when walking pipeline folders
- [x] Show and record run paths relative to the repository root
- [x] ruff configuration (PEP 8, Google docstrings, received code excluded)
- [x] Pin R exactly: `R = "= 4.4.2"` in `rproj.toml`, then `rig proj lock`
      (bfa99ee, whose message says 4.2.2 by mistake; lock confirmed current)
- [x] Commit `.gitignore` rig block; keep `.rvenvlib/` ignored (F19)
- [x] Commit `README.md` (Python and R sections) (510ebc9; line 37 changes with the
      `black` item below, the R-on-`PATH` section with F6)
- [ ] Drop `black` from the dev dependency group (`ruff format` replaces it)
- [ ] Make the pipeline use the pinned R rather than the first `Rscript` on `PATH` (F6)

## 2. Safeguards (before further changes)

- [ ] Keep the verified run's receipt (`pipeline/rebuild/RUN_INTEGRITY.json`,
      1,944 output hashes) and `results/` hashes as the regression baseline
- [ ] Inventory cross-file couplings (table below)
- [ ] Make source patches fail loudly when their target text is missing (F18)

## 3. Setup

- [ ] `run_all.py` — reviewed; F1, F7, F8 fixed; F2, F3 open
- [ ] `materialize.py` — resolve F13
- [ ] `prepare_primary_inputs.py`
- [ ] Configuration read during setup: `root_map.json`, `input_manifest.json`,
      `workspace_directories.json`, `stages.json`, `environment.json`

## 4. Analysis stages

- [ ] `run_stage.py` — reviewed; F4, F5, F8 fixed; F6 open
- [ ] `stage_guard.py`
- [ ] Map each stage's entry script and everything it imports or calls
      (47 stages, one entry script each; about 88,500 lines in
      `pipeline/code`, `extensions`, `renderers`)
- [ ] Identify scripts never executed; confirm against
      `validation/EXCLUDED_UNUSED_FRAGMENTS.json`; decide whether to ship them
- [ ] Trace each main-text figure and the molecular table to its producers
      (`FIGURE_MAP.csv`, `stages.json`) and review those producers in full
- [ ] Replace hardcoded `/tmp` and `/private/tmp` `MPLCONFIGDIR` values (F16)

## 5. Figures

- [ ] `prepare_display_tables.py`
- [ ] `prepare_render_inputs.py`
- [ ] `run_figures.py` — includes F17
- [ ] Configuration: `display_input_map.json`, `render_requests.json`,
      `render_jobs.json`
- [ ] Review `pipeline/renderers/` (24 scripts)

## 6. Run wrap-up and notebook layer

- [ ] `document_provenance.py`
- [ ] `verify_inputs.py`
- [ ] `run_integrity.py` — resolve F10, F11, F12
- [ ] `code/reproduce.py` — includes F15
- [ ] `code/extend_analysis.py` — source patching (F9, F18)
- [ ] Review `pipeline/extensions/` (16 scripts, run by `extend_analysis.py`)
- [ ] `code/collect_figures.py`
- [ ] `code/collect_source_tables.py`
- [ ] `code/rebuild_supplement.py`
- [ ] `code/supplement_draw_job.py`
- [ ] Remove source-text patching dependencies, format the numerical code in one
      commit, lift the ruff exclusion, and confirm a fresh run matches the baseline

## 7. Checks and entry points

- [ ] `code/check_primary_predictions.py`
- [ ] `code/check_expression_and_text.py`
- [ ] `code/check_correspondence.py`
- [ ] `code/check_genomewide.py`
- [ ] `code/check_integrity_guard.py`
- [ ] `code/audit_result_annotations.py`
- [ ] `execute_notebook.py`
- [ ] `verify_package.py`
- [ ] `pipeline/verify_code.py`
- [ ] `pipeline/validate_workflow.py`
- [ ] `Reproduce_ICM_v281.ipynb` — update cell 474fd814 (still cites
      `requirements-lock.txt`); review the resume heuristic (F15)
- [ ] Decide fate of legacy `pipeline/requirements.txt` and `requirements-lock.txt`

## 8. Documentation

- [ ] Merge `README.txt` into `README.md`; remove `README.txt`; point
      `readme` in `pyproject.toml` at `README.md`
- [ ] `ENVIRONMENT.json` — claims "pip check passed"; align with uv
- [ ] Check `COMPREHENSIVE_AUDIT.txt`, `RESULT_HARDCODING_AUDIT.txt`,
      `PREPROCESSING.txt` against the code as released
- [ ] Check `DATA_PROVENANCE.csv`, `DATA_LICENSES.txt`, `REFERENCES.txt`
- [ ] Decide what of `validation/` and `reference/` ships

## 9. Release

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

## Couplings

Places where one file depends on another's exact text, name or bytes. A change
to the right-hand file can silently break the left-hand one.

| Dependent | Depends on | Kind | Fails loudly? |
| --------- | ---------- | ---- | ------------- |
| `materialize.py` | historical absolute paths in `pipeline/code/**` (`root_map.json`) | text rewrite (`re.sub`) | No: an unmatched path is left unchanged |
| `stage_guard.py` | `root_map.json` roots | path checks at run time | To review |
| `code/extend_analysis.py:82` | paths in `pipeline/extensions/**` (`root_map.json`) | text rewrite (`re.sub`) | No |
| `code/extend_analysis.py:101` | `W=P.parents[2]` in `figure1_extension/analyze.py` | text patch (`str.replace`) | No: already a no-op (F9) |
| `run_figures.py` | renderer source, `pd.read_csv` and `np.load` | `exec` with patched I/O | To review |
| `verify_code.py` | 513 files in `SHA256_CODE.txt` | byte hashes | Yes |
| `run_integrity.py` | all `pipeline/` code and run outputs | byte hashes | Yes, on resume |
| `run_all.py`, `run_stage.py` | script paths in `stages.json` | file names | Yes |

## Findings

| ID | Where | Finding | Proposed action | Status |
| -- | ----- | ------- | --------------- | ------ |
| F1 | `run_all.py` | `pipeline/reports/` and `pipeline/audit/` are created only by `code/reproduce.py:59`; running `run_all.py` directly, as its docstring documents, fails after all 47 stages | Create the folders in `run_all.py` or in each writer | Fixed in `run_all.py` (4681a80) |
| F2 | `run_all.py:77` | `verify_inputs.py` runs after all analyses, too late to prevent a bad run (`materialize.py` already checks input hashes) | Move before `materialize`, or drop as redundant | Open |
| F3 | `run_all.py:52-69` | Hand-coded file copies and `prepare_display_tables.py` before two stages are dependencies missing from `stages.json` | Document; consider moving into `stages.json` | Open |
| F4 | `run_stage.py:41-44` | Stage output goes only to its log; a failure surfaces as a bare `CalledProcessError` | Print the log tail on failure | Fixed (2705207) |
| F5 | launchers | `assert` used for runtime checks (about 32), skipped under `python -O` | Raise explicit exceptions | `run_stage.py` fixed (025b166); other launchers open |
| F6 | `run_stage.py:35` | `Rscript` taken from `PATH`; R version not enforced | Check `R.version` or resolve R via rig | Open |
| F7 | `run_all.py:10` | Positional `sys.argv` parsing; `--resume` before the path is misread as the directory | Use `argparse` | Fixed (f8666f4) |
| F8 | launchers | Single-letter globals, multi-imports, mid-file import, no `main()`, missing docstrings, long lines | Restyle | `run_all.py` (f718f78), `run_stage.py` (4058727) done; other launchers open |
| F9 | `code/extend_analysis.py:101` | Patch `.replace("W=P.parents[2]", ...)` never matches: the target already reads `W = P.parents[1]` (formatted, see `validation/black_pipeline.log`) | Remove the dead patch | Open |
| F10 | `run_integrity.py` | Receipt does not hash `code/`, so edits there pass a resume check | Add `code/` to the code snapshot | Open |
| F11 | `run_integrity.py` | Receipt recorded twice (`run_all.py` and `code/reproduce.py`); the second overwrites the first | Record once | Open |
| F12 | `run_integrity.py:77` | Fingerprints every file in `data/source/atlas/` by glob rather than the declared inputs | Fingerprint the named CSVs only | Open |
| F13 | `materialize.py` | Unused `import shutil`; `assert` for hash checks (see F5) | Clean up | Open |
| F14 | many scripts | Thread variables set in about 89 places | Keep (single-threaded by design); consider one shared helper later | Open |
| F15 | `code/reproduce.py`, notebook | Resume chosen by `(REBUILD / "status").is_dir()`, so an interrupted run is resumed and then rejected | Base resume on the receipt's existence | Open |
| F16 | `pipeline/code`, `extensions`, `renderers` | Hardcoded `/tmp` and `/private/tmp` `MPLCONFIGDIR` | Point at the run's cache folder | Open |
| F17 | `run_figures.py` | Author-only `--discover` branch hardcodes another user's absolute paths | Remove before release | Open |
| F18 | source patches | `str.replace` and `re.sub` on source text do nothing when the target text is absent, so a later edit can disable a patch without any error | Route every patch through a helper that checks the expected number of matches | Open |
| F19 | `.gitignore`, `.rvenvlib/` | `/.rvenvlib/` is ignored inside rig's managed `# rig rvenv` block, so `rig proj lock`/`sync` (which `README.md` tells users to run) deletes the line and `.rvenvlib/` shows as untracked; the block's comment and `.rvenvlib/.gitignore` both say the shim is committed; a fresh clone warns "package 'rvenv' … was not found" until `rig proj sync` | Either commit `.rvenvlib/` as rig intends, or keep it ignored with the rule moved outside rig's markers and the comment corrected | Kept ignored; rule moved after rig's block, which survives `rig proj lock`/`sync` (bd77b97) |
| F20 | `run_integrity.py:114-126` | The R fingerprint runs `Rscript` in the caller's working directory, while the stage runs in `<run>/current`; launched from the repository root it reads `.Renviron` and the project library, so it could record a different `Matrix` than the stage loads (identical today: `rproj.lock` has no packages) | Run the probe from the same directory as the stage | Open |
| F21 | `run_integrity.py:106` | Fingerprints `python-igraph`, which is not installed (the distribution is `igraph` 1.0.0); the verified receipt records `"python-igraph": "not installed"`, so a changed igraph (Leiden clustering) passes a resume check | Fingerprint `igraph` | Open |
| F22 | working tree | 70 files differ from `PACKAGE_SHA256.json` and root `requirements.txt`, `requirements-lock.txt` are missing. `README.txt` and `Reproduce_ICM_v281.ipynb` were edited before review (uv migration); most others (`results/`, `pipeline/figure_notes/`, `data/analytic/`, some `validation/` CSVs) were rewritten by the verified run. The original archive (`~/Downloads/ICM_v281_collaborator_comprehensive.zip`) matches the manifest for all 1,196 files | Take each "as received" commit from the archive, not the working tree, then commit the local change separately; decide whether run-rewritten outputs ship as received or as regenerated | Open |
