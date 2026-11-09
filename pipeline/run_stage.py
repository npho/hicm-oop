"""Run an explicit analysis stage in the isolated workspace and retain its status.
Use materialize.py + prepare_primary_inputs.py first. This command never substitutes
an archived result if a producer fails. Failure is recorded and propagated.
"""

from pathlib import Path
import sys, json, os, subprocess, time, hashlib

B = Path(__file__).resolve().parent
R = Path(sys.argv[1]).resolve()
name = sys.argv[2]
stages = {s["name"]: s for s in json.loads((B / "stages.json").read_text())}
s = stages[name]
for dependency in s["depends"]:
    p = R / "status" / (dependency + ".json")
    if not p.exists() or json.loads(p.read_text()).get("status") != "PASS":
        raise RuntimeError("Prerequisite has not passed: " + dependency)
script = R / s["script"]
assert script.is_file(), script
(R / "logs").mkdir(exist_ok=True)
(R / "status").mkdir(exist_ok=True)
env = os.environ.copy()
env.update(
    ICM_ARCHIVE_ROOT=str(R / "archive"),
    ICM_WORK_ROOT=str(R / "archive/work"),
    ICM_OUTPUT_ROOT=str(R / "legacy"),
    OPENBLAS_NUM_THREADS="1",
    OMP_NUM_THREADS="1",
    MKL_NUM_THREADS="1",
    MPLCONFIGDIR=str(R / "cache/mpl"),
    NUMBA_CACHE_DIR=str(R / "cache/numba"),
    XDG_CACHE_HOME=str(R / "cache"),
)
command = (
    ["Rscript", str(script)]
    if script.suffix == ".R"
    else [sys.executable, str(B / "stage_guard.py"), str(R), str(script)]
)
start = time.time()
print("START", name, flush=True)
with (R / "logs" / (name + ".log")).open("w") as f:
    res = subprocess.run(
        command, cwd=R / "current", env=env, stdout=f, stderr=subprocess.STDOUT
    )
record = dict(
    stage=name,
    status="PASS" if res.returncode == 0 else "FAIL",
    exit_code=res.returncode,
    seconds=round(time.time() - start, 2),
    script=s["script"],
    script_sha256=hashlib.sha256(script.read_bytes()).hexdigest(),
    log="logs/" + name + ".log",
    scope=s["scope"],
)
(R / "status" / (name + ".json")).write_text(json.dumps(record, indent=2))
print(json.dumps(record), flush=True)
sys.exit(res.returncode)
