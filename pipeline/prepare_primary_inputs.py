"""Create the analysis subset from the published full count matrix.

Selects the ICM, preimplantation epiblast and hypoblast cells from the deposited
atlas and writes the counts and metadata that the network analysis reads. No
archived fit or cached model input is used.

Usage: python prepare_primary_inputs.py RUN_DIRECTORY
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PIPELINE = Path(__file__).resolve().parent
ATLAS = PIPELINE.parent / "data/source/atlas"
CELL_TYPES = ["ICM", "preIm-Epi", "Hypo"]


def load_primary_cells():
    """Read the primary cells' metadata and counts and check the counts.

    Returns:
        A tuple ``(metadata, counts)`` with one row per selected cell.

    Raises:
        ValueError: If the counts are not usable as they are stored.
    """
    metadata = pd.read_csv(ATLAS / "Human_Sample_Info.csv", index_col=0)
    metadata = metadata[metadata.Manual_Annotations.isin(CELL_TYPES)].copy()
    counts = pd.read_csv(ATLAS / "Human_Embryo_Counts.csv", index_col=0)
    counts = counts.loc[metadata.index]
    if not (counts.index.is_unique and counts.columns.is_unique):
        raise ValueError("Count matrix has duplicate cell or gene names")
    if not (counts.to_numpy() >= 0).all():
        raise ValueError("Count matrix has negative or missing counts")
    if not (counts.sum(axis=1) > 0).all():
        raise ValueError("A selected cell has no counts")
    if not (counts.to_numpy().astype(np.float32) == counts.to_numpy()).all():
        raise ValueError("Counts are not exactly representable as float32")
    return metadata, counts


def main():
    """Write the primary expression matrix and metadata into the run directory."""
    run_dir = Path(sys.argv[1]).resolve()
    output = run_dir / "archive/theory_v230/network_analysis/inputs"
    output.mkdir(parents=True, exist_ok=True)
    metadata, counts = load_primary_cells()
    np.savez_compressed(
        output / "expression.npz",
        counts=counts.values.astype(np.float32),
        genes=counts.columns.to_numpy(dtype=str),
        cells=counts.index.to_numpy(dtype=str),
    )
    metadata.to_csv(output / "metadata.csv")
    print(
        json.dumps(
            dict(
                cells=len(metadata),
                genes=len(counts.columns),
                ICM=int(metadata.Manual_Annotations.eq("ICM").sum()),
                source="Published complete count matrix; no cached model inputs",
            )
        )
    )


if __name__ == "__main__":
    main()
