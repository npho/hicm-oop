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
- [x] Drop `black` from the dev dependency group (`ruff format` replaces it)
      (1e5dd69; leftover constraint pins removed in 1e328aa)
- [x] Make the pipeline use the pinned R rather than the first `Rscript` on `PATH` (F6)

## 2. Safeguards (before further changes)

- [x] Keep the verified run's receipt (`pipeline/rebuild/RUN_INTEGRITY.json`,
      1,944 output hashes) and `results/` hashes as the regression baseline
      (`validation/baseline/`; full run archived outside the repo; F24)
- [x] Inventory cross-file couplings (table below; four rows added, F25–F27 found;
      the `run_figures.py` cell stays "To review" until §5)
- [x] Make source patches fail loudly when their target text is missing (F18)
      (launcher rewrites done; patching inside legacy scripts to be settled in §4)

## 3. Setup

- [x] `run_all.py` — reviewed; F1, F2, F7, F8 fixed; F3 closed as documented; F29 open (§5)
- [x] `materialize.py` — reviewed; F13 fixed, F5 and F8 done for it, F30 (a–c) fixed
- [x] `prepare_primary_inputs.py` — reviewed; F31 (a, b, d) fixed, (c) open (§6)
- [x] Configuration read during setup: `root_map.json`, `input_manifest.json`,
      `workspace_directories.json`, `stages.json`, `environment.json` — reviewed;
      F34 closed, F32 deferred to §8, F33 deferred to §9

## 4. Analysis stages

- [x] `run_stage.py` — reviewed; F4, F5, F6, F8 fixed
- [x] `stage_guard.py` — reviewed; F27, F35, F36, F37 fixed; F23 open (R stage runs unguarded)
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
- [ ] `ENVIRONMENT.json` — claims "pip check passed"; align with uv; it lists
      `black` and its dependencies (`click`, `mypy-extensions`, `pathspec`,
      `tomli`), which the lock no longer has
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
| `materialize.py` | historical absolute paths in `pipeline/code/**` (`root_map.json`) | text rewrite (`source_patch.rewrite_roots`) | Partly: each of the four roots must match a file, else it raises; an unmapped path is printed, not fatal (F28) |
| `stage_guard.py` | `root_map.json` roots | path checks at run time, on resolved paths | Yes: raises on a read under a historical root |
| `code/extend_analysis.py` | paths in `pipeline/extensions/**` (`root_map.json`) | text rewrite (`source_patch.rewrite_roots`) | Partly: `archive`, `legacy` and `cesfw` must each match a file (`current` never does), else it raises; unmapped paths are printed |
| `materialize.py`, `code/extend_analysis.py` | `source_patch.py` | import | Yes |
| `run_figures.py` | renderer source, `pd.read_csv` and `np.load` | `exec` with patched I/O | To review |
| `verify_code.py` | 513 files in `SHA256_CODE.txt` | byte hashes | Yes |
| `run_integrity.py` | all `pipeline/` code and run outputs | byte hashes | Yes, on resume |
| `run_all.py`, `run_stage.py` | script paths in `stages.json` | file names | Yes |
| `run_integrity.py` | `run_stage.r_versions` and `R_VERSION_PROBE` | import; receipt records the probe's output | Yes: a changed probe changes `runtime` and fails resume |
| `run_stage.py` `R_VERSION`, `MATRIX_VERSION` | `rproj.toml`, `ENVIRONMENT.json` | duplicated version numbers | Yes, at the R stage, if the installed R differs; No if only the files drift apart |
| `materialize.py`, `code/extend_analysis.py`, `stage_guard.py` | `root_map.json` (four historical roots) | one file, three consumers: two text rewrites and the run-time access block | Only `stage_guard.py` raises; the rewrites leave unmatched paths unchanged (F18) |
| `code/*.py` (11 files: `reproduce.py`, `extend_analysis.py`, the `check_*` scripts, collectors, `rebuild_supplement.py`, `audit_result_annotations.py`), `input_aliases.json` | the folder name `pipeline/rebuild` | hardcoded path (F25) | Mostly: a missing path raises; stale files in shared folders (`project/`, `figures/`, `reports/`) do not |
| `code/check_*`, collectors, `audit_result_annotations.py`; `run_figures.py` via `display_input_map.json` and `input_aliases.json` | exact output file names and folders written by the stages | hardcoded relative paths | `pd.read_csv` raises on a missing file, not checked per reader; `prepare_render_inputs.py` handling of a missing input: to review |
| `code/extend_analysis.py` | `run_all.py`'s receipt, and `status/*.json` | `verify` at start; skips a script whose status is PASS with the same script hash | Receipt: yes. Status: no, a stale PASS file skips the script silently. It also writes into the run after the receipt was recorded (F11) |

## Findings

| ID | Where | Finding | Proposed action | Status |
| -- | ----- | ------- | --------------- | ------ |
| F1 | `run_all.py` | `pipeline/reports/` and `pipeline/audit/` are created only by `code/reproduce.py:59`; running `run_all.py` directly, as its docstring documents, fails after all 47 stages | Create the folders in `run_all.py` or in each writer | Fixed in `run_all.py` (4681a80) |
| F2 | `run_all.py:77` | `verify_inputs.py` runs after all analyses, too late to prevent a bad run (`materialize.py` already checks input hashes) | Move before `materialize`, or drop as redundant | Fixed: runs once at the start (about 4 s); the end-of-run call is dropped, since the receipt records the input hashes |
| F3 | `run_all.py:52-69` | Hand-coded file copies and `prepare_display_tables.py` before two stages are dependencies missing from `stages.json` | Document; consider moving into `stages.json` | Closed as documented: both copies raise `FileNotFoundError` when a source is missing, so a missing dependency fails loudly; explained in the `prepare_stage` docstring; left in code. The double call is F29 |
| F4 | `run_stage.py:41-44` | Stage output goes only to its log; a failure surfaces as a bare `CalledProcessError` | Print the log tail on failure | Fixed (2705207) |
| F5 | launchers | `assert` used for runtime checks (about 32), skipped under `python -O` | Raise explicit exceptions | `run_stage.py` fixed (025b166), `materialize.py` and `prepare_primary_inputs.py` fixed; other launchers open |
| F6 | `run_stage.py:50` | `Rscript` taken from `PATH`; R version not enforced. The stage also inherits the caller's environment and, with no startup files in `<run>/current`, reads `~/.Rprofile`, `~/.Renviron` and the user library, which precedes R's bundled library even under `--vanilla`: a user-installed `Matrix` would silently replace 1.7.1. A wrong R is caught only by a later resume check, never on a fresh run | Check `R.version` or resolve R via rig; isolate R from user startup files and libraries | Fixed: `--vanilla`, bundled library only, version check before the stage; outputs byte-identical (baf1876) |
| F7 | `run_all.py:10` | Positional `sys.argv` parsing; `--resume` before the path is misread as the directory | Use `argparse` | Fixed (f8666f4) |
| F8 | launchers | Single-letter globals, multi-imports, mid-file import, no `main()`, missing docstrings, long lines | Restyle | `run_all.py` (f718f78), `run_stage.py` (4058727), `materialize.py`, `prepare_primary_inputs.py` done; other launchers open |
| F9 | `code/extend_analysis.py:101` | Patch `.replace("W=P.parents[2]", ...)` never matches: the target already reads `W = P.parents[1]` (formatted, see `validation/black_pipeline.log`) | Remove the dead patch | Fixed: patch removed; extension copies byte-identical |
| F10 | `run_integrity.py` | Receipt does not hash `code/`, so edits there pass a resume check | Add `code/` to the code snapshot | Open |
| F11 | `run_integrity.py` | Receipt recorded twice (`run_all.py` and `code/reproduce.py`); the second overwrites the first | Record once | Open |
| F12 | `run_integrity.py:77` | Fingerprints every file in `data/source/atlas/` by glob rather than the declared inputs | Fingerprint the named CSVs only | Open |
| F13 | `materialize.py` | Unused `import shutil`; `assert` for hash checks (see F5) | Clean up | Fixed: unused import removed, imports split, `assert`s replaced by explicit errors; workspace byte-identical (528 entries) |
| F14 | many scripts | Thread variables set in about 89 places | Keep (single-threaded by design); consider one shared helper later | Open |
| F15 | `code/reproduce.py`, notebook | Resume chosen by `(REBUILD / "status").is_dir()`, so an interrupted run is resumed and then rejected | Base resume on the receipt's existence | Open |
| F16 | `pipeline/code`, `extensions`, `renderers` | Hardcoded `/tmp` and `/private/tmp` `MPLCONFIGDIR` | Point at the run's cache folder | Open |
| F17 | `run_figures.py` | Author-only `--discover` branch hardcodes another user's absolute paths | Remove before release | Open |
| F18 | source patches | `str.replace` and `re.sub` on source text do nothing when the target text is absent, so a later edit can disable a patch without any error | Route every patch through a helper that checks the expected number of matches | Fixed for the launchers: `pipeline/source_patch.py` counts matches per root and raises if a required root matches nothing; materialized workspace (528 entries) and extension copies (16 files) byte-identical to before. Note: `re.sub` replaces every match, so a "root still present" check cannot fail and was not used. Scripts that patch other scripts' text (for example `legacy/scripts_v237/04_figures.py`) are not stage entry scripts; confirm in the §4 mapping |
| F19 | `.gitignore`, `.rvenvlib/` | `/.rvenvlib/` is ignored inside rig's managed `# rig rvenv` block, so `rig proj lock`/`sync` (which `README.md` tells users to run) deletes the line and `.rvenvlib/` shows as untracked; the block's comment and `.rvenvlib/.gitignore` both say the shim is committed; a fresh clone warns "package 'rvenv' … was not found" until `rig proj sync` | Either commit `.rvenvlib/` as rig intends, or keep it ignored with the rule moved outside rig's markers and the comment corrected | Kept ignored; rule moved after rig's block, which survives `rig proj lock`/`sync` (bd77b97) |
| F20 | `run_integrity.py:114-126` | The R fingerprint runs `Rscript` in the caller's working directory, while the stage runs in `<run>/current`; launched from the repository root it reads `.Renviron` and the project library, so it could record a different `Matrix` than the stage loads (identical today: `rproj.lock` has no packages) | Run the probe from the same directory as the stage | Fixed: shares `run_stage.r_versions`; recorded value unchanged (8baf231) |
| F21 | `run_integrity.py:106` | Fingerprints `python-igraph`, which is not installed (the distribution is `igraph` 1.0.0); the verified receipt records `"python-igraph": "not installed"`, so a changed igraph (Leiden clustering) passes a resume check | Fingerprint `igraph` | Open |
| F22 | working tree | 70 files differ from `PACKAGE_SHA256.json` and root `requirements.txt`, `requirements-lock.txt` are missing. `README.txt` and `Reproduce_ICM_v281.ipynb` were edited before review (uv migration); most others (`results/`, `pipeline/figure_notes/`, `data/analytic/`, some `validation/` CSVs) were rewritten by the verified run. The original archive (`~/Downloads/ICM_v281_collaborator_comprehensive.zip`) matches the manifest for all 1,196 files | Take each "as received" commit from the archive, not the working tree, then commit the local change separately; decide whether run-rewritten outputs ship as received or as regenerated | Open |
| F23 | `run_stage.py:50`, `stage_guard.py` | The R stage runs without `stage_guard.py`, so its reads are neither restricted to the run nor recorded; today the run's copy of `export_cultures.R` has every path rewritten into the run directory | Review with `stage_guard.py`; at least check the materialized R script for historical roots | Open |
| F24 | `pipeline/rebuild/RUN_INTEGRITY.json`, `results/` | The verified run's receipt is intact (1,944 of 1,944 outputs and 49 of 49 inputs match disk) but is the only copy, in an untracked folder that a fresh run needs moved aside; it holds hashes only, no content to diff against. It does not cover `results/` (145 files, rewritten by the run, F22), so there is no baseline for them. Its `code` section matches 540 of 546 files; the 6 that differ are the launchers edited in review (`run_all.py`, `run_stage.py`, `run_integrity.py`, `materialize.py`, `run_figures.py`, `prepare_render_inputs.py`), so the receipt cannot pass a resume check any more, by design. It records `python-igraph: not installed` (F21) | Keep a copy of the receipt and a `results/` hash list outside `pipeline/rebuild/`; optionally archive the run's outputs; compare a fresh run against these | Fixed: `validation/baseline/RUN_INTEGRITY.verified.json` (byte-identical copy) and `results.sha256` (145 files, all verified against disk); `~/hicm-oop-baseline/` holds both plus `rebuild-verified-run.tar.gz` (1.75 GB, the whole of `pipeline/rebuild/` minus `.DS_Store`) |
| F25 | `pipeline/rebuild/`, `run_all.py`, `code/*.py`, `input_aliases.json` | A run needs an empty directory so stale files cannot leak into a "fresh from declared inputs" run and the receipt covers one run only. Reruns are therefore not idempotent: the name `rebuild` is hardcoded in 11 files outside the pipeline (`code/` checks and collectors, `reproduce.py`), `input_aliases.json` records paths into it, and `pipeline/project/`, `figures/`, `cache/`, `reports/`, `audit/` are shared by all runs and never cleared, so a second run in a new folder would overwrite them | Define the run directory once and allow a timestamped name (needs the 11 scripts reviewed first); decide whether the shared folders move inside the run | Deferred by decision: no change for now; revisit in §6 |
| F26 | `run_figures.py` (job loop) | When a renderer's expected output was not recorded by the patched `savefig`, the fallback takes any file of that name under `pipeline/project/` modified within a second of the job's start. `project/` is shared across runs (F25), so a stale file could be picked up instead of raising | Review in §5; drop the fallback or restrict it to the job's own folder | Open (found in the couplings inventory) |
| F27 | `stage_guard.py` | The user-home read guard only triggers for paths starting `/Users/` or `/home/`; a package located elsewhere (for example under `/srv` or `/mnt`) would get no home-directory protection, and the check does not apply on Windows. The historical-root check compares path strings | Review in §4 together with F23; build the guard from the real home directory instead of two prefixes | Fixed (0fcb6f4): the guard also blocks the real `Path.home()`. The Windows case remains: not a target of this package |
| F28 | `pipeline/code/current/outputs/ICM_shared_precursor_validation/run.py:40`, `ICM_pathway_module_exploration/analyze.py:17` | Historical locations outside `root_map.json` are left unrewritten: `/Users/adelosan/.cache/gseapy/…` and, in 12 legacy build scripts, `/Users/adelosan/miniconda3/…` (pandoc, conda Pythons, bamCoverage). `materialize.py` now prints each one (14 lines per run). In the verified run it is a dead fallback (the stage read the in-run `Hallmark_2020.gmt`; its log shows `gmt exists True`), and `stage_guard.py` would block it loudly. `analyze.py` is not a stage | Leave until the §4 script mapping; then drop the fallback or declare the path | Open (found reviewing F18) |
| F29 | `run_all.py` (`prepare_stage` and `main`) | `prepare_display_tables.py` runs twice: inside `prepare_stage` before `figure3_display` and again after all stages. `figure3_display` is the last stage and the script does not mention `v264` or `figure3`, so the second run probably rewrites identical files, but this is unproven (153 lines, received numerical code) | Review the script in §5; compare the bytes it writes after each call on a scratch copy of a run; drop the second call if identical | Open |
| F30 | `materialize.py` | (a) Hash checks use `read_bytes()`, loading each input whole (largest 2.8 GB, `hema_adult_bmref.h5ad`) while `verify_inputs.py` reads 8 MB blocks; (b) a declared input marked unavailable is skipped silently (none today: 54 rows, 0 unavailable); (c) a run folder outside the package is detected only by the final `print`, after all the work, as a bare `ValueError`; (d) `read_text()`/`write_text()` use the platform's default encoding | (a) hash in blocks; (b) raise; (c) check at the start; (d) state the encoding (optional) | (a), (b), (c) fixed, workspace byte-identical; (d) not applied (no change on this machine) |
| F31 | `prepare_primary_inputs.py`; `code/reproduce.py` `prepare_flat_files` | (a) Data checks are `assert`s (F5), and the script lacks the non-zero-total check per cell that `prepare_flat_files` has; (b) the cast to `float32` is exact only because the counts are integers below 2^24 (true today: max 3,685,601), which nothing checks; (c) the same ICM/preIm-Epi/Hypo subset is built twice: here for the pipeline (`expression.npz`, `metadata.csv`) and in the notebook for `data/analytic/primary_*`, which no code reads; the metadata files are byte-identical (`b60d6039…`); (d) style: unused `hashlib`, multi-import, ambiguous name `O`, long line, no `main()` (F8). Tested: two runs 3 s apart and the verified run give the same hashes for `expression.npz` (`12370459…`) and `metadata.csv`, so the output is byte-reproducible | (a), (b) explicit errors plus the two checks (outputs unchanged); (d) restyle; (c) decide when `reproduce.py` is reviewed (§6): share one function or keep both | (a), (b), (d) fixed: `expression.npz` and `metadata.csv` byte-identical to the verified run; (c) open until `reproduce.py` is reviewed (§6) |
| F32 | `pipeline/environment.json` | No code reads it. It is a 19-package subset of the root `ENVIRONMENT.json` (about 95 packages, plus R, Matrix and a scope note), so two records of the tested environment exist and can drift. The names differ in case (`PyMuPDF`, `Pillow` against `pymupdf`, `pillow`). `pyproject.toml` and `uv.lock` are the real pins. | Decide with the `ENVIRONMENT.json` item in §8 whether to drop it or generate it from `uv.lock`. No change now. | open (§8) |
| F33 | `input_manifest.json` | All 54 rows carry `local_source` and `resolved_source`: the original author's home paths. No code reads them. They name a person and a machine in a release package, and they are the source of most of the unmapped-path lines `materialize.py` prints. | Release decision: strip the two fields (the file is hashed in the code receipt, not the numerical one). No change now. | open (§9) |
| F34 | `workspace_directories.json`, `stages.json` | Checked and sound: 387 unique sorted directories that cover every manifest destination; 47 unique stages whose dependencies all resolve, with no cycle, and whose scripts all exist under `pipeline/code`. Two small points: four entries are `__pycache__` folders left by the author's runs (empty directories, harmless), and the stored stage order is not topological (for example `regulatory_prediction` precedes `regulatory_protocol`), so the run order comes from `run_stages` taking the first eligible stage in file order. That order is deterministic, so keep it. | None needed. Optionally drop the four `__pycache__` lines. | closed, no change |
| F35 | `stage_guard.py` | Style (F8) and a small cost: single-letter globals `B`, `R`, `seen`, multi-import line, no docstrings or `main()`, `sys.argv` read at import, and the allowed-root list is rebuilt with `resolve()` calls on every `open` event | Restyle with the behaviour unchanged; compute the allowed roots once | Fixed (defe606); behaviour identical in a probe of ten cases under the old and new file |
| F36 | `stage_guard.py:39-44` | The historical-root test compared the unresolved path as a string, so a symlink in the run pointing into a historical root was not blocked (the home check on line 35 caught it only while the root sits under `/Users/` or `/home/`). The `PIPELINE` and `sys.prefix` clauses are exemptions for a package or interpreter located inside a historical root: harmless, so kept | Test the resolved path against resolved roots | Fixed (0fcb6f4); probe: a symlink into a historical root is now blocked, every other case unchanged |
| F37 | `stage_guard.py` | The access log records only Python `open` events: not h5py or other C-library reads, not subprocesses (R), not `listdir` or `stat`. The module docstring says it records "project file accesses". In the verified run 41 of 54 declared inputs appear in the logs; the other 13 include `hema_adult_bmref.h5ad`, `troph_donor_p13.h5ad`, `adata.h5ad` and `M2021BlastoidIMs.rds`, which are read by C libraries or R, or not at all (the log cannot tell which). Network blocking likewise covers only this process | Say so in the docstring and here; do not use the logs to prove an input was unused | Fixed by documentation (defe606); limits stay |
