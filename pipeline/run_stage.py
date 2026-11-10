"""Run one analysis stage in the isolated workspace and record its status.

Usage: python run_stage.py RUN_DIRECTORY STAGE_NAME

Run materialize.py and prepare_primary_inputs.py first. A failing producer is
never replaced by an archived result: the failure is recorded in
status/STAGE_NAME.json and returned as this script's exit code.
"""

import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

PIPELINE = Path(__file__).resolve().parent
LOG_TAIL_LINES = 30

# The tested R; keep in step with rproj.toml and ENVIRONMENT.json. R 4.4.2
# bundles Matrix 1.7.1, the only package the R stage uses.
R_VERSION = "4.4.2"
MATRIX_VERSION = "1.7.1"
# Prints R.version.string, then the Matrix version. The integrity receipt
# records this output, so changing it invalidates existing receipts.
R_VERSION_PROBE = (
    'cat(R.version.string,"\\n"); cat(as.character(packageVersion("Matrix")))'
)


def check_dependencies(run_dir, stage):
    """Raise unless every stage this one depends on has passed."""
    for dependency in stage["depends"]:
        path = run_dir / "status" / (dependency + ".json")
        if not path.exists() or json.loads(path.read_text()).get("status") != "PASS":
            raise RuntimeError("Prerequisite has not passed: " + dependency)


def stage_environment(run_dir):
    """Return the stage environment: workspace roots, one thread, local caches."""
    env = os.environ.copy()
    env.update(
        ICM_ARCHIVE_ROOT=str(run_dir / "archive"),
        ICM_WORK_ROOT=str(run_dir / "archive/work"),
        ICM_OUTPUT_ROOT=str(run_dir / "legacy"),
        OPENBLAS_NUM_THREADS="1",
        OMP_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        MPLCONFIGDIR=str(run_dir / "cache/mpl"),
        NUMBA_CACHE_DIR=str(run_dir / "cache/numba"),
        XDG_CACHE_HOME=str(run_dir / "cache"),
    )
    return env


def r_environment(run_dir):
    """Return the stage environment limited to R's bundled package library.

    Every user and site library variable points at one folder in the run's
    cache that holds no packages (R creates the R_LIBS_USER folder if it is
    missing), so only the library R ships with can supply a package.
    """
    env = stage_environment(run_dir)
    empty = str(run_dir / "cache/R/empty-library")
    env.update(R_LIBS=empty, R_LIBS_USER=empty, R_LIBS_SITE=empty)
    return env


def r_command(*args):
    """Return an Rscript command that ignores all .Rprofile and .Renviron files."""
    return ["Rscript", "--vanilla", *args]


def r_versions(run_dir):
    """Return the R and Matrix versions exactly as the R stage would see them."""
    result = subprocess.run(
        r_command("-e", R_VERSION_PROBE),
        cwd=run_dir / "current",
        env=r_environment(run_dir),
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.strip()


def check_r(run_dir):
    """Raise unless the Rscript on PATH is the tested R with its bundled Matrix."""
    found = r_versions(run_dir)
    lines = found.splitlines()
    r_ok = bool(lines) and lines[0].startswith(f"R version {R_VERSION} ")
    if not r_ok or lines[-1] != MATRIX_VERSION:
        raise RuntimeError(
            f"Expected R {R_VERSION} with Matrix {MATRIX_VERSION}; the Rscript on "
            f"PATH reports {found!r}. With rig, run `rig default 4.4`."
        )


def stage_command(run_dir, script):
    """Return the command that runs an R script, or a Python script under the guard."""
    if script.suffix == ".R":
        return r_command(str(script))
    guard = PIPELINE / "stage_guard.py"
    return [sys.executable, str(guard), str(run_dir), str(script)]


def print_log_tail(run_dir, log):
    """Print the end of a failed stage's log to stderr, where callers see it."""
    lines = (run_dir / log).read_text(errors="replace").splitlines()
    tail = lines[-LOG_TAIL_LINES:]
    print(f"--- last {len(tail)} lines of {log} ---", file=sys.stderr)
    print("\n".join(tail), file=sys.stderr, flush=True)


def main():
    """Run the named stage, write its log and status, and return its exit code."""
    run_dir = Path(sys.argv[1]).resolve()
    name = sys.argv[2]
    stages = json.loads((PIPELINE / "stages.json").read_text())
    stage = {stage["name"]: stage for stage in stages}[name]
    check_dependencies(run_dir, stage)
    script = run_dir / stage["script"]
    if not script.is_file():
        raise FileNotFoundError(f"Stage script not found: {script}")
    if script.suffix == ".R":
        check_r(run_dir)
        env = r_environment(run_dir)
    else:
        env = stage_environment(run_dir)
    (run_dir / "logs").mkdir(exist_ok=True)
    (run_dir / "status").mkdir(exist_ok=True)

    start = time.time()
    print("START", name, flush=True)
    with (run_dir / "logs" / (name + ".log")).open("w") as log:
        result = subprocess.run(
            stage_command(run_dir, script),
            cwd=run_dir / "current",
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    record = {
        "stage": name,
        "status": "PASS" if result.returncode == 0 else "FAIL",
        "exit_code": result.returncode,
        "seconds": round(time.time() - start, 2),
        "script": stage["script"],
        "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
        "log": "logs/" + name + ".log",
        "scope": stage["scope"],
    }
    (run_dir / "status" / (name + ".json")).write_text(json.dumps(record, indent=2))
    print(json.dumps(record), flush=True)
    if result.returncode != 0:
        print_log_tail(run_dir, record["log"])
    return result.returncode


if __name__ == "__main__":
    sys.exit(main())
