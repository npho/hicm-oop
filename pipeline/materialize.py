"""Create an isolated analysis workspace from code and declared inputs only.

No reference result table, fitted dictionary or archived figure is copied in.
Historical absolute paths are rewritten mechanically, with every replacement
logged.

Usage: python materialize.py RUN_DIRECTORY
"""

import hashlib
import json
import os
import sys
from collections import Counter
from pathlib import Path

from source_patch import require_roots, rewrite_roots, unmapped_paths

PIPELINE = Path(__file__).resolve().parent


def copy_code(run_dir):
    """Copy pipeline/code into the run with historical roots rewritten.

    Args:
        run_dir: Run directory that receives the rewritten copy.

    Returns:
        The files whose text changed, for path_replacements.json.
    """
    roots = json.loads((PIPELINE / "root_map.json").read_text())
    mapping = {old: str(run_dir / name) for name, old in roots.items()}
    matches = Counter()
    changes = []
    unmapped = {}
    for source in (PIPELINE / "code").rglob("*"):
        if not source.is_file() or source.name.startswith("."):
            continue
        destination = run_dir / source.relative_to(PIPELINE / "code")
        destination.parent.mkdir(parents=True, exist_ok=True)
        original = source.read_text()
        text = rewrite_roots(original, mapping, matches)
        destination.write_text(text)
        relative = str(source.relative_to(PIPELINE))
        if text != original:
            changes.append(
                dict(
                    file=relative,
                    operation=(
                        "Replace historical root paths with isolated workspace roots"
                    ),
                )
            )
        if left := unmapped_paths(text, PIPELINE.parent):
            unmapped[relative] = left
    require_roots(matches, roots.values(), "pipeline/code")
    for name, paths in unmapped.items():
        print("Unmapped historical path left in", name, paths)
    return changes


def link_inputs(run_dir):
    """Link every declared input into the run, checking its hash first.

    Args:
        run_dir: Run directory that receives the links.
    """
    for row in json.loads((PIPELINE / "input_manifest.json").read_text()):
        if not row["available"]:
            continue
        source = PIPELINE / row["packaged"]
        assert hashlib.sha256(source.read_bytes()).hexdigest() == row["sha256"], source
        destination = run_dir / row["root"] / row["path"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() or destination.is_symlink():
            assert (
                destination.resolve() == source.resolve()
                or hashlib.sha256(destination.read_bytes()).hexdigest() == row["sha256"]
            )
            continue
        # A relative symlink avoids doubling large files and stays portable with
        # the package directory.
        destination.symlink_to(os.path.relpath(source, destination.parent))


def main():
    """Build the workspace named on the command line."""
    run_dir = Path(sys.argv[1]).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    changes = copy_code(run_dir)
    link_inputs(run_dir)
    directories = json.loads((PIPELINE / "workspace_directories.json").read_text())
    for relative in directories:
        (run_dir / relative).mkdir(parents=True, exist_ok=True)
    (run_dir / "path_replacements.json").write_text(json.dumps(changes, indent=2))
    print(
        "Workspace ready:",
        run_dir.relative_to(PIPELINE.parent),
        "— only code and declared inputs supplied",
    )


if __name__ == "__main__":
    main()
