"""Create the analysis subset from the published full count matrix, not archived fits."""

from pathlib import Path
import sys, numpy as np, pandas as pd, json, hashlib

B = Path(__file__).resolve().parent
R = Path(sys.argv[1]).resolve()
O = R / "archive/theory_v230/network_analysis/inputs"
O.mkdir(parents=True, exist_ok=True)
m = pd.read_csv(B.parent / "data/source/atlas/Human_Sample_Info.csv", index_col=0)
m = m[m.Manual_Annotations.isin(["ICM", "preIm-Epi", "Hypo"])].copy()
c = pd.read_csv(
    B.parent / "data/source/atlas/Human_Embryo_Counts.csv", index_col=0
).loc[m.index]
if not (c.index.is_unique and c.columns.is_unique):
    raise ValueError("Count matrix has duplicate cell or gene names")
if not (c.to_numpy() >= 0).all():
    raise ValueError("Count matrix has negative or missing counts")
if not (c.sum(axis=1) > 0).all():
    raise ValueError("A selected cell has no counts")
if not (c.to_numpy().astype(np.float32) == c.to_numpy()).all():
    raise ValueError("Counts are not exactly representable as float32")
np.savez_compressed(
    O / "expression.npz",
    counts=c.values.astype(np.float32),
    genes=c.columns.to_numpy(dtype=str),
    cells=c.index.to_numpy(dtype=str),
)
m.to_csv(O / "metadata.csv")
print(
    json.dumps(
        dict(
            cells=len(m),
            genes=len(c.columns),
            ICM=int(m.Manual_Annotations.eq("ICM").sum()),
            source="Published complete count matrix; no cached model inputs",
        )
    )
)
