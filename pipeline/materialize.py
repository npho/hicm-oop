"""Create an isolated analysis workspace using code + declared input files only.
No reference result table, fitted dictionary or archived figure is copied in.
Historical absolute paths are rewritten mechanically, with every replacement logged.
"""

from pathlib import Path
import sys, json, shutil, os, hashlib
from collections import Counter

from source_patch import require_roots, rewrite_roots, unmapped_paths

B = Path(__file__).resolve().parent
R = Path(sys.argv[1]).resolve()
R.mkdir(parents=True, exist_ok=True)
roots = json.loads((B / "root_map.json").read_text())
changes = []
mapping = {v: str(R / k) for k, v in roots.items()}
matches = Counter()
unmapped = {}
for p in (B / "code").rglob("*"):
    if not p.is_file() or p.name.startswith("."):
        continue
    dest = R / p.relative_to(B / "code")
    dest.parent.mkdir(parents=True, exist_ok=True)
    text = p.read_text()
    old = text
    text = rewrite_roots(text, mapping, matches)
    dest.write_text(text)
    if text != old:
        changes.append(
            dict(
                file=str(p.relative_to(B)),
                operation="Replace historical root paths with isolated workspace roots",
            )
        )
    if left := unmapped_paths(text, B.parent):
        unmapped[str(p.relative_to(B))] = left
require_roots(matches, roots.values(), "pipeline/code")
for name, paths in unmapped.items():
    print("Unmapped historical path left in", name, paths)
for row in json.loads((B / "input_manifest.json").read_text()):
    if not row["available"]:
        continue
    p = B / row["packaged"]
    assert hashlib.sha256(p.read_bytes()).hexdigest() == row["sha256"], p
    dest = R / row["root"] / row["path"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        assert (
            dest.resolve() == p.resolve()
            or hashlib.sha256(dest.read_bytes()).hexdigest() == row["sha256"]
        )
        continue
    # Relative symlink avoids doubling large files; portable with the package directory.
    dest.symlink_to(os.path.relpath(p, dest.parent))
for rel in json.loads((B / "workspace_directories.json").read_text()):
    (R / rel).mkdir(parents=True, exist_ok=True)
(R / "path_replacements.json").write_text(json.dumps(changes, indent=2))
print(
    "Workspace ready:",
    R.relative_to(B.parent),
    "— only code and declared inputs supplied",
)
