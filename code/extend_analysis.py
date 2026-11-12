"""Recalculate post-v272 analyses inside the isolated processed-data workspace.

The original numerical implementations are preserved. Only filesystem roots and
obsolete live-manuscript guards are changed. No archived results substitute for
failed analysis stages.
"""

from pathlib import Path
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time

PACKAGE = Path(__file__).resolve().parents[1]
ENGINE = PACKAGE / "pipeline"
WORK = ENGINE / "rebuild"
CURRENT = WORK / "current"
OUTPUTS = CURRENT / "outputs"
ROOTS = json.loads((ENGINE / "root_map.json").read_text())


def copy_file(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def run(script):
    """Execute a producer offline and record success only after process exit."""
    name = script.parent.name + "__" + script.stem
    status_path = WORK / "status" / (name + ".json")
    script_hash = hashlib.sha256(script.read_bytes()).hexdigest()
    if status_path.exists():
        previous = json.loads(status_path.read_text())
        if (
            previous.get("status") == "PASS"
            and previous.get("script_sha256") == script_hash
        ):
            print("Already generated in this run:", name, flush=True)
            return
    log = WORK / "logs" / (name + ".log")
    start = time.time()
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with log.open("w") as stream:
        result = subprocess.run(
            [sys.executable, str(ENGINE / "stage_guard.py"), str(WORK), str(script)],
            cwd=CURRENT,
            env=env,
            stdout=stream,
            stderr=subprocess.STDOUT,
        )
    status = {
        "stage": name,
        "script_sha256": script_hash,
        "status": "PASS" if result.returncode == 0 else "FAIL",
        "seconds": round(time.time() - start, 2),
        "log": str(log.relative_to(PACKAGE)),
    }
    (WORK / "status" / (name + ".json")).write_text(json.dumps(status, indent=2))
    print(status, flush=True)
    if result.returncode:
        raise RuntimeError(f"Analysis failed; inspect {log}")


def main():
    """Build flat tables, refit models, then prepare current figure inputs."""
    sys.path.insert(0, str(ENGINE))
    from run_integrity import verify as verify_run_integrity

    verify_run_integrity(ENGINE, WORK)
    # Copy and relocate only code; never copy cached fits into the workspace.
    mapping = {old: str(WORK / key) for key, old in ROOTS.items()}
    pattern = "|".join(re.escape(key) for key in sorted(mapping, key=len, reverse=True))
    for source in (ENGINE / "extensions").rglob("*.py"):
        dest = OUTPUTS / source.relative_to(ENGINE / "extensions")
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(
            re.sub(pattern, lambda m: mapping[m.group(0)], source.read_text())
        )
    for name in [
        "ICM_naive_correspondence_check",
        "ICM_blastoid_correspondence_check",
        "ICM_iyer_correspondence_check",
        "ICM_culture_correspondence",
        "ICM_addback_controls",
    ]:
        if name == "ICM_addback_controls" and "--correspondence-only" in sys.argv:
            return
        run(OUTPUTS / name / "run.py")
    run(OUTPUTS / "ICM_blastoid_marker_check" / "check.py")

    # Prepare Figure 1 single-cell gene values from newly made reference counts.
    panel = OUTPUTS / "figure1_extension"
    for filename in ["pseudotime_cells.csv", "paga_edges.csv"]:
        copy_file(WORK / "legacy/results_v242" / filename, panel / "inputs" / filename)
    source = panel / "analyze.py"
    run(source)
    run(panel / "render.py")

    # Figure 3 uses new matched-cell losses, not the archived published estimates.
    panel = OUTPUTS / "figure3_current"
    for source, filename in [
        (
            OUTPUTS / "results_v264/source_data/Figure3_displayed_values.csv",
            "Figure3_displayed_values.csv",
        ),
        (
            OUTPUTS / "ICM_neighbor_population_test/matched/summary.csv",
            "neighbor_summary.csv",
        ),
        (OUTPUTS / "ICM_addback_controls/summary.csv", "ICM_addback_controls.csv"),
    ]:
        copy_file(source, panel / "inputs" / filename)
    run(panel / "render.py")

    # Figure 5: all fractions derive from freshly fitted embryo-reference models.
    panel = OUTPUTS / "figure5_current"
    mapping = {
        "naive_summary.csv": OUTPUTS
        / "ICM_naive_correspondence_check/culture_summary.csv",
        "validation_summary.csv": OUTPUTS
        / "ICM_naive_correspondence_check/validation_summary.csv",
        "blastoid_summary.csv": OUTPUTS
        / "ICM_blastoid_correspondence_check/culture_summary.csv",
        "iyer_summary.csv": OUTPUTS
        / "ICM_iyer_correspondence_check/culture_summary.csv",
    }
    for filename, source in mapping.items():
        copy_file(source, panel / "inputs" / filename)
    run(panel / "render.py")
    panel = OUTPUTS / "figure5_markers"
    for filename in ["gene_summary.csv", "threshold_sensitivity.csv"]:
        copy_file(
            OUTPUTS / "ICM_blastoid_marker_check" / filename,
            panel / "inputs" / filename,
        )
    run(panel / "render.py")

    # Figure 6 gene values are extracted from deposited processed/scaled h5ad.
    panel = OUTPUTS / "figure6_current"
    (panel / "inputs").mkdir(exist_ok=True)
    for filename in [
        "culture_summary.csv",
        "marker_profiles.csv",
        "validation_summary.csv",
    ]:
        copy_file(
            OUTPUTS / "ICM_culture_correspondence" / filename,
            panel / "inputs" / filename,
        )
    for filename in ["projection_summary.csv", "culture_metadata.csv"]:
        copy_file(
            OUTPUTS / "results_v253" / filename,
            panel / "inputs/results_v253" / filename,
        )
    run(panel / "extract_gene_timecourses.py")
    run(panel / "render_figure6.py")

    # v281 systematic analysis: freshly reconstructed primary analytic matrix.
    target = OUTPUTS / "ICM_systematic_panel_check"
    for filename in ["expression.npz", "metadata.csv"]:
        copy_file(
            WORK / "archive/theory_v230/network_analysis/inputs" / filename,
            target / "inputs" / filename,
        )
    copy_file(
        PACKAGE / "data/analysis/systematic_panel/gene_sets.json",
        target / "gene_sets.json",
    )
    run(OUTPUTS / "ICM_joint_daughter_remodeling/analyze.py")
    run(OUTPUTS / "ICM_joint_daughter_remodeling/report.py")
    run(OUTPUTS / "render_joint_panels.py")
    (PACKAGE / "results").mkdir(exist_ok=True)
    (PACKAGE / "results/extension_status.json").write_text(
        json.dumps({"status": "PASS"}, indent=2)
    )


if __name__ == "__main__":
    main()
    from run_integrity import record as record_run_integrity

    record_run_integrity(ENGINE, WORK)
