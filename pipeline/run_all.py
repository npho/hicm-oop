"""Rebuild every analysis and figure from the declared processed inputs.

Usage: python run_all.py RUN_DIRECTORY [--resume]

The run directory must sit beneath this package and be new or empty, unless
--resume is given for a completed run whose integrity receipt still matches.
Archived results never substitute for a failed stage. Every completed stage
leaves a log and a status file in the run directory.
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path

from run_integrity import record as record_run_integrity
from run_integrity import verify as verify_run_integrity

PIPELINE = Path(__file__).resolve().parent

# Tables the figure3_display stage reads from results_v263/source_data, given
# relative to current/outputs. Each copy is named after its source path with
# "/" replaced by "__".
FIGURE3_DISPLAY_INPUTS = [
    "results_v239/source_data/precursor_prediction.csv",
    "results_v239/source_data/heldout_usage_contrasts.csv",
    "results_v239/source_data/precursor_robustness.csv",
    "ICM_state_alternatives/summary.csv",
    "ICM_state_alternatives/transfer_summary.csv",
    "ICM_claim_audit/summary.csv",
    "ICM_claim_audit/exact_bracketing/summary.csv",
]


def run(script, *args):
    """Run a pipeline script with this interpreter, raising if it fails."""
    subprocess.run(
        [sys.executable, str(PIPELINE / script), *map(str, args)], check=True
    )


def copy(src, dst):
    """Copy a file with its metadata, creating the destination folder."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def passed_stages(run_dir):
    """Return the names of stages whose status file records a pass."""
    return {
        path.stem
        for path in (run_dir / "status").glob("*.json")
        if json.loads(path.read_text()).get("status") == "PASS"
    }


def prepare_stage(run_dir, name):
    """Stage inputs that stages.json does not declare (see REVIEW.md F3)."""
    outputs = run_dir / "current/outputs"
    if name == "remodeling_calibration":
        copy(
            run_dir / "legacy/results_v233/remodeling/gene_effects.csv",
            outputs / "results_v239/source_data/remodeling_gene_effects.csv",
        )
    if name == "figure3_display":
        run("prepare_display_tables.py", run_dir)
        destination = outputs / "results_v263/source_data"
        for source in FIGURE3_DISPLAY_INPUTS:
            copy(outputs / source, destination / source.replace("/", "__"))


def run_stages(run_dir, stages, done):
    """Run every stage not in done, each once its dependencies have passed.

    Adds each completed stage's name to done.
    """
    pending = [stage for stage in stages if stage["name"] not in done]
    while pending:
        eligible = next(
            (stage for stage in pending if set(stage["depends"]) <= done), None
        )
        if eligible is None:
            raise RuntimeError("Unresolved stage dependencies: " + repr(pending))
        name = eligible["name"]
        prepare_stage(run_dir, name)
        run("run_stage.py", run_dir, name)
        done.add(name)
        pending.remove(eligible)


def main():
    """Materialize the workspace, run all stages and figures, and record a receipt."""
    run_dir = Path(sys.argv[1]).resolve()
    if PIPELINE not in run_dir.parents:
        raise SystemExit("Create the run directory beneath the package directory.")
    resume = "--resume" in sys.argv

    if resume:
        verify_run_integrity(PIPELINE, run_dir)
    if run_dir.exists() and any(run_dir.iterdir()) and not resume:
        raise SystemExit(
            "Use an empty run directory or explicitly --resume this generated run."
        )

    if not resume:
        run("materialize.py", run_dir)
        run("prepare_primary_inputs.py", run_dir)
    stages = json.loads((PIPELINE / "stages.json").read_text())
    done = passed_stages(run_dir) if resume else set()
    run_stages(run_dir, stages, done)

    run("prepare_display_tables.py", run_dir)
    run("prepare_render_inputs.py", run_dir)
    run("run_figures.py", "all")
    run("document_provenance.py", run_dir)
    run("verify_inputs.py")
    last_run = {
        "analysis_stages": sum(stage["name"] in done for stage in stages),
        "run_directory": str(run_dir.relative_to(PIPELINE.parent)),
        "status": "PASS",
        "scope": (
            "Fresh analyses from declared processed inputs plus all figure drawing jobs"
        ),
    }
    (PIPELINE / "reports/LAST_RUN.json").write_text(json.dumps(last_run, indent=2))
    record_run_integrity(PIPELINE, run_dir)
    print("All declared analysis and rendering stages passed.")


if __name__ == "__main__":
    main()
