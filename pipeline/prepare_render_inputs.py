"""Resolve renderer inputs exclusively to newly generated analysis tables."""

from pathlib import Path
import sys, json, shutil, os

B = Path(__file__).resolve().parent
R = Path(sys.argv[1]).resolve()
m = json.loads((B / "display_input_map.json").read_text())
req = json.loads((B / "render_requests.json").read_text())
aliases = {}
missing = []
for ref, src in m.items():
    p = R / src
    if not p.is_file():
        missing.append(str(p))
        continue
    if ref.startswith("project/"):
        dest = B / ref
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, dest)
for key, ref in req.items():
    # Keep the {PACKAGE} placeholder on disk; run_figures.py expands it on load.
    aliases[key] = os.path.relpath(R / m[ref], B)
(B / "input_aliases.json").write_text(json.dumps(aliases, indent=2))
(B / "reports/render_input_preparation.json").write_text(
    json.dumps({"missing": missing, "mapping": m}, indent=2)
)
if missing:
    raise FileNotFoundError("\n".join(missing))
