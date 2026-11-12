"""Block accidental reads of the author's historical workspace and network calls.
Record project file accesses. Libraries/runtime resources remain readable.
"""

from pathlib import Path
import sys, os, json, runpy, atexit

B = Path(__file__).resolve().parent
R = Path(sys.argv[1]).resolve()
script = Path(sys.argv[2]).resolve()
roots = list(json.loads((B / "root_map.json").read_text()).values())
seen = set()


def guard(event, args):
    if event in ["socket.connect", "socket.getaddrinfo"]:
        raise RuntimeError("Network is disabled during the clean analysis run")
    if event == "open" and args and isinstance(args[0], (str, bytes)):
        p = Path(os.fsdecode(args[0])).absolute()
        s = str(p)
        # Permit package inputs, this run and the Python installation; reject other
        # user-home data even when an old script contains an undeclared fallback.
        resolved = p.resolve()
        if resolved.name == ".igraphrc":
            raise FileNotFoundError(
                "Use igraph defaults; user configuration is not an input"
            )
        allowed = [
            B,
            B.parent / "data",
            R,
            Path(sys.prefix).resolve(),
            Path(sys.base_prefix).resolve(),
        ]
        if str(resolved).startswith(("/Users/", "/home/")) and not any(
            resolved == a or a in resolved.parents for a in allowed
        ):
            raise RuntimeError("Undeclared user-home access blocked: " + str(resolved))
        if (
            any(s.startswith(v + "/") for v in roots)
            and B not in p.parents
            and not any(s.startswith(v + "/") for v in [sys.prefix, sys.base_prefix])
        ):
            raise RuntimeError("Historical workspace access blocked: " + s)
        if R in p.parents:
            seen.add((s, str(args[1]) if len(args) > 1 else ""))


sys.addaudithook(guard)


def save():
    out = R / "access_logs"
    out.mkdir(exist_ok=True)
    (out / (script.parent.name + "__" + script.name + ".json")).write_text(
        json.dumps(sorted(seen), indent=2)
    )


atexit.register(save)
sys.path.insert(0, str(script.parent))
sys.argv = [str(script)]
runpy.run_path(str(script), run_name="__main__")
