"""Run one stage script under a read guard and record the files it opens.

Usage: python stage_guard.py RUN_DIRECTORY SCRIPT

The guard blocks network connections, the igraph user configuration, reads
under a user's home directory that belong neither to the package, the run nor
the Python installation, and reads of the author's historical workspace roots.
Every file the script opens inside the run is recorded in
RUN_DIRECTORY/access_logs.

Only Python open() calls are blocked and recorded. Reads by C libraries such as
h5py, by subprocesses such as R, and directory listings are neither, so the
logs cannot show that a declared input went unused.
"""

import atexit
import json
import os
import runpy
import sys
from pathlib import Path

PIPELINE = Path(__file__).resolve().parent
HISTORICAL_ROOTS = list(json.loads((PIPELINE / "root_map.json").read_text()).values())


def install_guard(run_dir):
    """Install the audit hook that enforces the guard.

    Args:
        run_dir: Run directory whose file accesses are recorded.

    Returns:
        The set of ``(path, mode)`` pairs the guard has recorded, filled as the
        script runs.
    """
    seen = set()
    allowed = [
        PIPELINE,
        PIPELINE.parent / "data",
        run_dir,
        Path(sys.prefix).resolve(),
        Path(sys.base_prefix).resolve(),
    ]

    def guard(event, args):
        if event in ["socket.connect", "socket.getaddrinfo"]:
            raise RuntimeError("Network is disabled during the clean analysis run")
        if event != "open" or not args or not isinstance(args[0], str | bytes):
            return
        path = Path(os.fsdecode(args[0])).absolute()
        text = str(path)
        # Permit package inputs, this run and the Python installation; reject other
        # user-home data even when an old script contains an undeclared fallback.
        resolved = path.resolve()
        if resolved.name == ".igraphrc":
            raise FileNotFoundError(
                "Use igraph defaults; user configuration is not an input"
            )
        if str(resolved).startswith(("/Users/", "/home/")) and not any(
            resolved == root or root in resolved.parents for root in allowed
        ):
            raise RuntimeError("Undeclared user-home access blocked: " + str(resolved))
        if (
            any(text.startswith(root + "/") for root in HISTORICAL_ROOTS)
            and PIPELINE not in path.parents
            and not any(text.startswith(p + "/") for p in [sys.prefix, sys.base_prefix])
        ):
            raise RuntimeError("Historical workspace access blocked: " + text)
        if run_dir in path.parents:
            seen.add((text, str(args[1]) if len(args) > 1 else ""))

    sys.addaudithook(guard)
    return seen


def save_access_log(run_dir, script, seen):
    """Write the recorded accesses to the run's access_logs folder."""
    out = run_dir / "access_logs"
    out.mkdir(exist_ok=True)
    name = script.parent.name + "__" + script.name + ".json"
    (out / name).write_text(json.dumps(sorted(seen), indent=2))


def main():
    """Run the script named on the command line under the guard."""
    run_dir = Path(sys.argv[1]).resolve()
    script = Path(sys.argv[2]).resolve()
    seen = install_guard(run_dir)
    atexit.register(save_access_log, run_dir, script, seen)
    sys.path.insert(0, str(script.parent))
    sys.argv = [str(script)]
    runpy.run_path(str(script), run_name="__main__")


if __name__ == "__main__":
    main()
