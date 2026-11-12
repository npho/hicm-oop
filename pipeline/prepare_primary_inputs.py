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
assert c.index.is_unique and c.columns.is_unique and (c.to_numpy() >= 0).all()
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
