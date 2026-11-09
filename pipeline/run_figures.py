"""Render v272 figures using freshly generated analysis tables. Usage: python run_figures.py [job|all].
Input resolution is explicit and fails on ambiguity. --discover is an author-only
archive construction option; reviewers use the default, which never reads outside
the archive. This is display reproduction, not a fresh raw-data analysis.
"""

from pathlib import Path
import os, sys, json, hashlib, shutil, traceback, time

ROOT = Path(__file__).resolve().parent
# Isolate jobs: Matplotlib rcParams and random-generator state must not leak
# from one historical plotting script into the next.
_requested = [v for v in sys.argv[1:] if not v.startswith("--")]
_jobs = json.loads((ROOT / "render_jobs.json").read_text())
_selected = list(_jobs) if not _requested or _requested == ["all"] else _requested
if len(_selected) > 1:
    import subprocess

    _results = []
    for _job in _selected:
        _run = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), _job]
            + (["--discover"] if "--discover" in sys.argv else [])
        )
        _status = ROOT / "audit" / ("render_status_" + _job + ".json")
        _results.extend(
            json.loads(_status.read_text())
            if _status.exists()
            else [dict(job=_job, status="FAIL", error="No status file")]
        )
    (ROOT / "audit/render_status_all.json").write_text(json.dumps(_results, indent=2))
    sys.exit(1 if any(r["status"] != "PASS" for r in _results) else 0)
PROJECT = ROOT / "project"
os.chdir(PROJECT)
os.environ.update(
    MPLCONFIGDIR=str(ROOT / "cache/matplotlib"),
    XDG_CACHE_HOME=str(ROOT / "cache"),
    OPENBLAS_NUM_THREADS="1",
    OMP_NUM_THREADS="1",
)
(ROOT / "cache/matplotlib").mkdir(parents=True, exist_ok=True)
JOBS = json.loads((ROOT / "render_jobs.json").read_text())
DISCOVER = False
ALIAS = ROOT / "input_aliases.json"
_PACKAGE = "{PACKAGE}"


def _portable(key):
    """Replace this machine's package path with the {PACKAGE} placeholder."""
    if key == str(ROOT) or key.startswith(str(ROOT) + os.sep):
        return _PACKAGE + key[len(str(ROOT)) :]
    return key


aliases = (
    {
        k.replace(_PACKAGE, str(ROOT)): v
        for k, v in json.loads(ALIAS.read_text()).items()
    }
    if ALIAS.exists()
    else {}
)
used = []
BASE = Path("/Users/adelosan/Documents/Codex/2026-09-26/what-i")


def resolve(path):
    if not isinstance(path, (str, Path)):
        return path
    p = Path(path)
    key = str(p)
    if key in aliases:
        q = ROOT / aliases[key]
    elif p.is_file() and (not p.is_absolute() or ROOT in p.parents):
        q = p
    else:
        parts = list(p.parts)
        candidates = []
        for i, part in enumerate(parts):
            if part.startswith(("results_v", "ICM_", "v258_exploration")):
                rel = Path(*parts[i:])
                candidates += [
                    PROJECT / "outputs" / rel,
                    PROJECT / "source_data" / rel,
                    PROJECT / rel,
                ]
        found = [x for x in candidates if x.is_file()]
        if found:
            q = found[0]
        else:
            # Only use a basename match when every candidate has identical bytes.
            found = list((PROJECT / "source_data").rglob(p.name))
            found = [f for f in found if f.is_file()]
            if (
                found
                and len({hashlib.sha256(f.read_bytes()).hexdigest() for f in found})
                == 1
            ):
                q = found[0]
            elif DISCOVER:
                choices = [
                    p,
                    BASE / p,
                    Path("/Users/adelosan/Documents/New project/icm_v234") / p,
                ]
                if parts and parts[0] == "outputs":
                    choices += [
                        Path("/Users/adelosan/Documents/New project/icm_v234")
                        / Path(*parts[1:])
                    ]
                src = next((f for f in choices if f.is_file()), None)
                if src is None:
                    raise FileNotFoundError("Unresolved input: " + key)
                digest = hashlib.sha256(src.read_bytes()).hexdigest()
                q = ROOT / "inputs" / digest[:12] / src.name
                q.parent.mkdir(parents=True, exist_ok=True)
                _copy(src, q)
            else:
                raise FileNotFoundError("Input absent/ambiguous in archive: " + key)
    aliases[key] = str(q.resolve().relative_to(ROOT))
    used.append(
        dict(
            request=_portable(key),
            file=aliases[key],
            sha256=hashlib.sha256(q.read_bytes()).hexdigest(),
        )
    )
    return q


import pandas as pd, numpy as np

_read_csv = pd.read_csv
_load = np.load
_copy = shutil.copy2
_to_csv = pd.DataFrame.to_csv
pd.read_csv = lambda path, *a, **k: _read_csv(resolve(path), *a, **k)
np.load = lambda path, *a, **k: _load(resolve(path), *a, **k)


def copy(src, dst, *a, **k):
    src = resolve(src)
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.resolve() == dst.resolve():
        return dst
    return _copy(src, dst, *a, **k)


shutil.copy2 = copy


def to_csv(self, path=None, *a, **k):
    if isinstance(path, (str, Path)):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    return _to_csv(self, path, *a, **k)


pd.DataFrame.to_csv = to_csv
import matplotlib

matplotlib.use("Agg")
from matplotlib.figure import Figure

if not DISCOVER:

    def offline_guard(event, args):
        if event.startswith("socket.connect") or event in ("socket.getaddrinfo",):
            raise RuntimeError("Network is disabled in figure replay")
        if event == "open" and args and isinstance(args[0], (str, bytes)):
            p = Path(os.fsdecode(args[0])).resolve()
            if str(p).startswith("/Users/adelosan/Documents/") and not any(
                root == p or root in p.parents
                for root in (
                    ROOT,
                    ROOT.parent / "data",
                    Path(sys.prefix).resolve(),
                    Path(sys.base_prefix).resolve(),
                )
            ):
                raise RuntimeError("External historical data access blocked: " + str(p))

    sys.addaudithook(offline_guard)
_savefig = Figure.savefig
written = []


def savefig(self, path, *a, **k):
    if isinstance(path, (str, Path)):
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        if p.is_absolute() and ROOT not in p.parents:
            raise ValueError("Output outside archive: " + str(p))
        written.append(str(p.resolve()))
    return _savefig(self, path, *a, **k)


Figure.savefig = savefig
names = [s for s in sys.argv[1:] if not s.startswith("--")]
names = list(JOBS) if not names or names == ["all"] else names
results = []
for job in names:
    used.clear()
    written.clear()
    start = time.time()
    status = "PASS"
    error = ""
    dest = ROOT / "figures" / job
    dest.mkdir(parents=True, exist_ok=True)
    try:
        src = ROOT / "renderers" / (job + ".py")
        code = src.read_text()
        env = {
            "__file__": str(PROJECT / "scripts" / (job + ".py")),
            "__name__": "__main__",
        }
        exec(compile(code, str(src), "exec"), env)
        for name in JOBS[job]["outputs"]:
            fs = [Path(f) for f in written if Path(f).name == name and Path(f).exists()]
            if not fs:
                fs = [f for f in PROJECT.rglob(name) if f.stat().st_mtime >= start - 1]
            if not fs:
                raise FileNotFoundError("Expected plot not generated: " + name)
            _copy(fs[-1], dest / name)
    except Exception as ex:
        status = "FAIL"
        error = repr(ex)
        traceback.print_exc()
    finally:
        ALIAS.write_text(
            json.dumps({_portable(k): v for k, v in aliases.items()}, indent=2)
        )
        (ROOT / "audit" / ("inputs_" + job + ".json")).write_text(
            json.dumps(used, indent=2)
        )
        import matplotlib.pyplot as plt

        plt.close("all")
    row = dict(
        job=job,
        status=status,
        error=error,
        seconds=round(time.time() - start, 2),
        input_count=len(used),
    )
    results.append(row)
    print(json.dumps(row), flush=True)
(
    ROOT
    / "audit"
    / ("render_status_" + ("_".join(names) if len(names) < 3 else "all") + ".json")
).write_text(json.dumps(results, indent=2))
if any(r["status"] != "PASS" for r in results):
    sys.exit(1)
