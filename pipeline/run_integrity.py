"""Fail closed before reusing a completed analysis workspace.

A receipt records actual code, declared input and generated numerical-file
hashes, plus installed package versions. It is written after successful work.
It is not a proof of correctness and does not validate upstream provenance.
"""

from pathlib import Path
import hashlib
import importlib.metadata
import json
import platform
import shutil

from run_stage import r_versions

NUMERICAL_SUFFIXES = {
    ".csv",
    ".tsv",
    ".npy",
    ".npz",
    ".h5ad",
    ".h5",
    ".mtx",
    ".rds",
    ".gz",
    ".json",
    ".txt",
    ".gmt",
    ".py",
    ".r",
    ".sh",
}


def digest(path, memo):
    """Hash bytes; reuse only within this call for duplicate resolved inputs."""
    resolved = path.resolve()
    if resolved not in memo:
        h = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                h.update(block)
        memo[resolved] = h.hexdigest()
    return memo[resolved]


def snapshot(engine, workspace):
    engine, workspace = Path(engine).resolve(), Path(workspace).resolve()
    memo = {}
    code = {}
    for base in [
        engine / "code",
        engine / "extensions",
        engine / "renderers",
        engine.parent / "code",
    ]:
        for path in sorted(base.rglob("*")):
            if (
                path.is_file()
                and path.suffix in {".py", ".R", ".sh", ".json"}
                and "__pycache__" not in path.parts
            ):
                code[str(path.relative_to(engine.parent))] = digest(path, memo)
    for path in sorted(engine.iterdir()):
        if (
            path.is_file()
            and path.suffix in {".py", ".json"}
            and path.name not in {"input_aliases.json"}
        ):
            code[str(path.relative_to(engine.parent))] = digest(path, memo)
    inputs = {}
    for row in json.loads((engine / "input_manifest.json").read_text()):
        if row["available"]:
            path = engine / row["packaged"]
            key = str(path.resolve().relative_to(engine.parent))
            inputs[key] = digest(path, memo)
    for path in sorted((engine.parent / "data/source/atlas").glob("*")):
        if path.is_file() and not path.name.startswith("."):
            inputs[str(path.relative_to(engine.parent))] = digest(path, memo)
    gene_sets = engine.parent / "data/analysis/systematic_panel/gene_sets.json"
    if gene_sets.exists():
        inputs[str(gene_sets.relative_to(engine.parent))] = digest(gene_sets, memo)
    outputs = {}
    for path in sorted(workspace.rglob("*")):
        rel = path.relative_to(workspace)
        if (
            not path.is_file()
            or path.is_symlink()
            or path.suffix.lower() not in NUMERICAL_SUFFIXES
        ):
            continue
        if any(part in {"__pycache__", "logs", "access_logs"} for part in rel.parts):
            continue
        if path.name == "RUN_INTEGRITY.json":
            continue
        outputs[str(rel)] = digest(path, memo)
    runtime = {"python": platform.python_version()}
    for name in [
        "numpy",
        "scipy",
        "pandas",
        "scikit-learn",
        "scanpy",
        "anndata",
        "h5py",
        "python-igraph",
        "umap-learn",
        "matplotlib",
    ]:
        try:
            runtime[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            runtime[name] = "not installed"
    if shutil.which("Rscript"):
        # Same command, directory and libraries as the R stage (run_stage.py).
        runtime["R_and_Matrix"] = r_versions(workspace)
    else:
        runtime["R_and_Matrix"] = "not installed"
    return {
        "code": code,
        "inputs": inputs,
        "generated_numerical_files": outputs,
        "runtime": runtime,
    }


def record(engine, workspace):
    receipt = snapshot(engine, workspace)
    (Path(workspace) / "RUN_INTEGRITY.json").write_text(json.dumps(receipt, indent=2))
    return receipt


def verify(engine, workspace):
    path = Path(workspace) / "RUN_INTEGRITY.json"
    if not path.exists():
        raise RuntimeError(
            "This workspace has no completed-run integrity receipt. Use a fresh empty run directory; do not reuse unverified cached results."
        )
    before = json.loads(path.read_text())
    after = snapshot(engine, workspace)
    changed = []
    for category in after:
        old, new = before.get(category, {}), after[category]
        changed.extend(
            f"{category}: {name}"
            for name in sorted(set(old) | set(new))
            if old.get(name) != new.get(name)
        )
    if changed:
        raise RuntimeError(
            "Cached run differs from its recorded inputs, code, runtime or outputs. Start a fresh run. Differences: "
            + "; ".join(changed[:12])
        )
    return True
