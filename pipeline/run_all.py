"""Clean processed-input rebuild. Never substitutes archived results.
Usage: python run_all.py /absolute/path/to/new_run_directory
A new directory is required. Every completed stage has a log and status file.
"""

from pathlib import Path
import sys, subprocess, json, shutil

B = Path(__file__).resolve().parent
R = Path(sys.argv[1]).resolve()
if B not in R.parents:
    raise SystemExit("Create the run directory beneath the package directory.")
resume = "--resume" in sys.argv
from run_integrity import verify as verify_run_integrity, record as record_run_integrity

if resume:
    verify_run_integrity(B, R)
if R.exists() and any(R.iterdir()) and not resume:
    raise SystemExit(
        "Use an empty run directory or explicitly --resume this generated run."
    )


def run(script, *args):
    subprocess.run([sys.executable, str(B / script), *map(str, args)], check=True)


def copy(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


if not resume:
    run("materialize.py", R)
    run("prepare_primary_inputs.py", R)
stages = json.loads((B / "stages.json").read_text())
done = (
    {
        p.stem
        for p in (R / "status").glob("*.json")
        if json.loads(p.read_text()).get("status") == "PASS"
    }
    if resume
    else set()
)
pending = [s for s in stages if s["name"] not in done]
while pending:
    eligible = next((s for s in pending if set(s["depends"]) <= done), None)
    if eligible is None:
        raise RuntimeError("Unresolved stage dependencies: " + repr(pending))
    name = eligible["name"]
    if name == "remodeling_calibration":
        copy(
            R / "legacy/results_v233/remodeling/gene_effects.csv",
            R / "current/outputs/results_v239/source_data/remodeling_gene_effects.csv",
        )
    if name == "figure3_display":
        run("prepare_display_tables.py", R)
        dest = R / "current/outputs/results_v263/source_data"
        for out, src in {
            "results_v239__source_data__precursor_prediction.csv": "current/outputs/results_v239/source_data/precursor_prediction.csv",
            "results_v239__source_data__heldout_usage_contrasts.csv": "current/outputs/results_v239/source_data/heldout_usage_contrasts.csv",
            "results_v239__source_data__precursor_robustness.csv": "current/outputs/results_v239/source_data/precursor_robustness.csv",
            "ICM_state_alternatives__summary.csv": "current/outputs/ICM_state_alternatives/summary.csv",
            "ICM_state_alternatives__transfer_summary.csv": "current/outputs/ICM_state_alternatives/transfer_summary.csv",
            "ICM_claim_audit__summary.csv": "current/outputs/ICM_claim_audit/summary.csv",
            "ICM_claim_audit__exact_bracketing__summary.csv": "current/outputs/ICM_claim_audit/exact_bracketing/summary.csv",
        }.items():
            copy(R / src, dest / out)
    run("run_stage.py", R, name)
    done.add(name)
    pending.remove(eligible)
run("prepare_display_tables.py", R)
run("prepare_render_inputs.py", R)
run("run_figures.py", "all")
run("document_provenance.py", R)
run("verify_inputs.py")
(B / "reports/LAST_RUN.json").write_text(
    json.dumps(
        {
            "analysis_stages": sum(s["name"] in done for s in stages),
            "run_directory": str(R),
            "status": "PASS",
            "scope": "Fresh analyses from declared processed inputs plus all figure drawing jobs",
        },
        indent=2,
    )
)
record_run_integrity(B, R)
print("All declared analysis and rendering stages passed.")
