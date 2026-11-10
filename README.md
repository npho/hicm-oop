# ICM Supporting Code

Reproduction of the manuscript figures, supplementary figures, and molecular table.

## Python

The Python environment is managed with [uv](https://docs.astral.sh/uv/). Three files define it:

| File              | Purpose                                                                 |
| ----------------- | ----------------------------------------------------------------------- |
| `.python-version` | Pins the interpreter to Python 3.10.14.                                 |
| `pyproject.toml`  | Declares the exact-pinned dependencies and holds transitive pins.       |
| `uv.lock`         | Records the fully resolved package set, with hashes, for every platform. |

The interpreter is pinned to the patch version on purpose: run fingerprints
(`pipeline/run_integrity.py`) record the exact Python version, so a different
interpreter prevents reuse of previously verified runs. The locked package set matches
the tested environment recorded in `ENVIRONMENT.json`.

### 1. Install uv

```sh
curl -LsSf https://astral.sh/uv/install.sh | sh   # or: brew install uv
```

### 2. Create the environment

From the repository root:

```sh
uv sync --locked
```

This downloads Python 3.10.14 if it is not already available, creates `.venv/`, and
installs exactly the packages in `uv.lock`. `--locked` makes the command fail rather
than re-resolve if `pyproject.toml` and `uv.lock` disagree. Add `--no-dev` to skip the
formatting and linting tool (`ruff`), which the analysis does not use.

### 3. Check the environment

```sh
uv run python -c "import platform; print(platform.python_version())"   # 3.10.14
uv pip check
```

### Running commands

Prefix commands with `uv run` to use the project environment without activating it:

```sh
uv run python execute_notebook.py
```

To open the notebook in an editor such as VS Code or JupyterLab, select the
interpreter at `.venv/bin/python`; its `python3` kernel is installed with the
environment.

### Changing dependencies

Edit the pins in `pyproject.toml`, then run `uv lock` and `uv sync --locked`. Do not
change versions casually: the pins reproduce the environment the numerical results
were validated against.

## R

One pipeline stage, `culture_export`
(`pipeline/code/legacy/scripts_v233/export_cultures.R`), uses R to extract counts from
the supplied Seurat object. It needs only the `Matrix` package. The R version is
managed with [rig](https://rig.r-lib.org/), and its project files
([`rig proj`](https://rig.r-lib.org/reference/proj.html), currently experimental)
define the environment:

| File         | Purpose                                                                       |
| ------------ | ----------------------------------------------------------------------------- |
| `rproj.toml` | Pins R to exactly 4.4.2.                                                      |
| `rproj.lock` | Records the resolved environment; no packages beyond those bundled with R.    |
| `.Renviron`  | Points R sessions started in the repository root at the project library.      |

`Matrix` is not declared because R 4.4.2 ships `Matrix` 1.7.1 as a recommended package,
which is the tested version recorded in `ENVIRONMENT.json`. Pinning R pins `Matrix`.

### 1. Install rig

```sh
brew tap r-lib/rig
brew install --cask rig
```

### 2. Install R 4.4.2

```sh
rig add 4.4.2
```

rig installs the official CRAN build alongside any other R versions (one patch release
per minor version on macOS). Installing asks for an administrator password.

### 3. Create the environment

From the repository root:

```sh
rig proj sync
```

This creates the machine-specific `.rvenv/` and fills in `.rvenvlib/`, rig's small
startup package that `.Renviron` loads. Until it has run, R sessions started in the
repository warn that package `rvenv` was not found.

### 4. Check the environment

```sh
Rscript -e 'cat(R.version.string, "\n"); packageVersion("Matrix")'   # R 4.4.2, 1.7.1
```

### How the pipeline finds R

The pipeline runs the `Rscript` that is first on `PATH` (with rig, its default R), from
inside its run directory. It runs it with `--vanilla` and with only the package library
bundled with R, so no `.Rprofile`, `.Renviron`, project library or user-installed
package can change the result. Before the R stage starts, the pipeline checks for R
4.4.2 with `Matrix` 1.7.1 and stops otherwise. Make sure rig's default is 4.4:

```sh
rig list          # the default is marked with *
rig default 4.4   # if another version is the default
```

### Changing the R version

Edit `rproj.toml`, then run `rig proj lock` and `rig proj sync`. Use the `= x.y.z` form
for an exact pin: a bare `x.y.z` means "x.y.z or any later compatible version", so
`R = "4.4.2"` would also accept R 4.6.
