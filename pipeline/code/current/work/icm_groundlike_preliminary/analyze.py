from pathlib import Path
import json, hashlib
import numpy as np, pandas as pd
from scipy.stats import rankdata

O = Path("outputs/ICM_groundlike_preliminary")
P = Path(
    "/Users/adelosan/Documents/New project/icm_v234/results_v237/developmental_projection"
)
states = ["ICM/TE Branch", "ICM", "Epi/Hypo Branch", "preIm-Epi", "Hypo"]
lowgenes = ["MYC", "DNMT3A", "DNMT3B", "UHRF1"]
pluri = ["NANOG", "SOX2"]
hypo = ["GATA6", "PDGFRA"]
cycle = [
    "MCM5",
    "PCNA",
    "TYMS",
    "FEN1",
    "MCM2",
    "MCM4",
    "RRM1",
    "UNG",
    "GINS2",
    "MCM6",
    "CDK1",
    "NUSAP1",
    "UBE2C",
    "BIRC5",
    "TPX2",
    "TOP2A",
    "NDC80",
    "CKS2",
    "MKI67",
    "CCNB2",
]
protocol = {
    "status": "Exploratory; specified before inspecting joint cell results. No manuscript edits.",
    "primary_coexpression": "At least one of NANOG/SOX2 AND at least one of GATA6/PDGFRA exceed 1 CPM in the SAME cell.",
    "primary_low": "Mean weighted expression percentile of MYC/DNMT3A/DNMT3B/UHRF1 <= 0.5. Gene percentiles within study; reference weights equalize five compartments, embryos within compartment, then cells within embryo-compartment. This is relative expression, not biochemical activity.",
    "primary_joint": "Primary coexpression AND primary low.",
    "sensitivities": [
        "Co-detection at >0 and >10 CPM",
        "NANOG AND GATA6 >1 CPM",
        "At least two of NANOG/SOX2/KLF17/TFCP2L1 AND GATA6 or PDGFRA >1 CPM",
        "All four low genes at/below their weighted median",
        "Exclude MYC from low score",
        "Require both GATA6 and PDGFRA",
        "Inspect GATA4/SOX17 and MYCN/MYCL separately",
    ],
    "summaries": "Equal-embryo average of within-compartment cell fractions, separate studies.",
    "uncertainty": "1000 whole-embryo bootstrap replicates, joint sampling of all compartments in an embryo; thresholds held fixed; descriptive conditional percentile intervals, no significance selection.",
    "technical_control": "Weighted linear regressions of low score and joint indicator on compartment, log1p(total counts), detected-gene fraction, fixed-panel cell-cycle score. Equal-compartment/equal-embryo weights. Contrasts ICM minus each neighbor, no causal interpretation.",
    "cycle_genes": cycle,
    "seed": 254,
}
(O / "protocol.json").write_text(json.dumps(protocol, indent=2))
z = np.load(P / "atlas_queries.npz")
m = pd.read_csv(P / "atlas_metadata.csv", index_col=0)
assert np.array_equal(m.index, z["cells"])
g = pd.Index(z["genes"])
x = z["counts"]
assert (x >= 0).all()
keep = m.Manual_Annotations.isin(states)
m = m.loc[keep].copy()
x = x[keep]
tot = x.sum(1)
cpm = x / tot[:, None] * 1e6
log = np.log2(1 + cpm)
m["total_counts"] = tot
m["detected_fraction"] = (x > 0).mean(1)
m["study"] = m.Dataset
m["state"] = m.Manual_Annotations
m["embryo"] = m["group"]
for gene in sorted(
    set(lowgenes + pluri + hypo + ["KLF17", "TFCP2L1", "GATA4", "SOX17", "MYCN", "MYCL"])
):
    m[gene] = cpm[:, g.get_loc(gene)]
# Library-relative cycle score: average within-cell gene percentile, no use of target genes.
cc = [g.get_loc(a) for a in cycle if a in g]
m["cycle_score"] = rankdata(log, axis=1, method="average")[:, cc].mean(1) / log.shape[1]
weights = np.zeros(len(m))
scores = {a: np.zeros(len(m)) for a in lowgenes}
threshold = []
for study in m.study.unique():
    sel = np.where(m.study.to_numpy() == study)[0]
    d = m.iloc[sel]
    n = d.groupby(["state", "embryo"]).state.transform("size")
    ne = d.groupby("state").embryo.transform("nunique")
    w = 1 / (len(states) * n.to_numpy() * ne.to_numpy())
    weights[sel] = w
    for gene in lowgenes:
        vals = d[gene].to_numpy()
        order = np.argsort(vals)
        v = vals[order]
        ww = w[order]
        u, starts = np.unique(v, return_index=True)
        mass = np.add.reduceat(ww, starts)
        cdf = np.cumsum(mass) / mass.sum()
        mid = cdf - mass / mass.sum() / 2
        scores[gene][sel] = mid[np.searchsorted(u, vals)]
        med = u[np.searchsorted(cdf, 0.5)]
        threshold.append(dict(study=study, gene=gene, median_CPM=med))
for a in lowgenes:
    m[a + "_percentile"] = scores[a]
m["weight"] = weights
m["low_score"] = np.mean([scores[a] for a in lowgenes], axis=0)
m["low_noMYC"] = np.mean([scores[a] for a in lowgenes[1:]], axis=0) <= 0.5
m["low"] = m.low_score <= 0.5
m["low_all4"] = True
for row in threshold:
    sel = m.study == row["study"]
    m.loc[sel, "low_all4"] = m.loc[sel, "low_all4"] & (
        m.loc[sel, row["gene"]] <= row["median_CPM"]
    )
for t, label in [(0, "detected"), (1, "primary"), (10, "stronger")]:
    m["co_" + label] = (m[pluri].max(axis=1) > t) & (m[hypo].max(axis=1) > t)
    m["joint_" + label] = m["co_" + label] & m.low
m["co_NANOG_GATA6"] = (m.NANOG > 1) & (m.GATA6 > 1)
m["joint_NANOG_GATA6"] = m.co_NANOG_GATA6 & m.low
m["co_extended"] = ((m[["NANOG", "SOX2", "KLF17", "TFCP2L1"]] > 1).sum(axis=1) >= 2) & (
    m[hypo].max(axis=1) > 1
)
m["joint_extended"] = m.co_extended & m.low
m["joint_all4"] = m.co_primary & m.low_all4
m["joint_noMYC"] = m.co_primary & m.low_noMYC
m["joint_bothHYPO"] = (m[pluri].max(axis=1) > 1) & (m[hypo].min(axis=1) > 1) & m.low
metrics = [
    "co_primary",
    "low",
    "joint_primary",
    "joint_detected",
    "joint_stronger",
    "joint_NANOG_GATA6",
    "joint_extended",
    "joint_all4",
    "joint_noMYC",
    "joint_bothHYPO",
    "low_score",
]
e = m.groupby(["study", "embryo", "state"])[metrics].mean().reset_index()
rows = []
for (study, state), d in e.groupby(["study", "state"]):
    for metric in metrics:
        rows.append(
            dict(
                study=study,
                state=state,
                metric=metric,
                value=d[metric].mean(),
                n_embryos=len(d),
                n_cells=int(((m.study == study) & (m.state == state)).sum()),
            )
        )
summary = pd.DataFrame(rows)
rng = np.random.default_rng(254)
boot = []
for study in m.study.unique():
    d = e[e.study == study]
    emb = d.embryo.unique()
    pivot = {
        met: d.pivot(index="embryo", columns="state", values=met)
        .reindex(index=emb, columns=states)
        .to_numpy()
        for met in metrics
    }
    ix = rng.integers(0, len(emb), (1000, len(emb)))
    for metric, a in pivot.items():
        with np.errstate(invalid="ignore"):
            b = np.nanmean(a[ix], axis=1)
        for j, state in enumerate(states):
            sel = (
                (summary.study == study)
                & (summary.state == state)
                & (summary.metric == metric)
            )
            summary.loc[sel, "lo"], summary.loc[sel, "hi"] = np.nanquantile(
                b[:, j], [0.025, 0.975]
            )
        if metric in ["joint_primary", "low_score", "co_primary"]:
            for j, state in enumerate(states):
                if state == "ICM":
                    continue
                diff = b[:, 1] - b[:, j]
                boot.append(
                    dict(
                        study=study,
                        metric=metric,
                        comparison="ICM minus " + state,
                        effect=np.nanmean(a[:, 1]) - np.nanmean(a[:, j]),
                        lo=np.nanquantile(diff, 0.025),
                        hi=np.nanquantile(diff, 0.975),
                        valid=int(np.isfinite(diff).sum()),
                    )
                )
# Technical/cycle adjustment; bootstrap full regression with embryo blocks.
adj = []
for study in m.study.unique():
    d = m[m.study == study].copy()
    w = d.weight.to_numpy()
    cov = np.column_stack(
        [np.log1p(d.total_counts), d.detected_fraction, d.cycle_score]
    )
    cov = (cov - np.average(cov, axis=0, weights=w)) / np.sqrt(
        np.average((cov - np.average(cov, axis=0, weights=w)) ** 2, axis=0, weights=w)
    )
    A = np.column_stack(
        [(d.state == state).astype(float) for state in states]
        + [cov[:, i] for i in range(3)]
    )
    emb = d.embryo.unique()
    ei = pd.Index(emb).get_indexer(d.embryo)
    for metric in ["joint_primary", "low_score"]:
        y = d[metric].to_numpy(float)

        def fit(wt):
            if np.linalg.matrix_rank(A * np.sqrt(wt[:, None])) < A.shape[1]:
                return np.full(A.shape[1], np.nan)
            return np.linalg.lstsq(
                A * np.sqrt(wt[:, None]), y * np.sqrt(wt), rcond=None
            )[0]

        coef = fit(w)
        bs = []
        for _ in range(1000):
            freq = np.bincount(rng.integers(0, len(emb), len(emb)), minlength=len(emb))
            bs.append(fit(w * freq[ei]))
        bs = np.array(bs)
        for j, state in enumerate(states):
            if state == "ICM":
                continue
            delta = bs[:, 1] - bs[:, j]
            adj.append(
                dict(
                    study=study,
                    metric=metric,
                    comparison="ICM minus " + state,
                    effect=coef[1] - coef[j],
                    lo=np.nanquantile(delta, 0.025),
                    hi=np.nanquantile(delta, 0.975),
                    valid=int(np.isfinite(delta).sum()),
                )
            )
m.to_csv(O / "cell_scores.csv")
e.to_csv(O / "embryo_fractions.csv", index=False)
summary.to_csv(O / "summary.csv", index=False)
pd.DataFrame(boot).to_csv(O / "unadjusted_contrasts.csv", index=False)
pd.DataFrame(adj).to_csv(O / "adjusted_contrasts.csv", index=False)
pd.DataFrame(threshold).to_csv(O / "thresholds.csv", index=False)
(O / "inputs.json").write_text(
    json.dumps(
        {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [P / "atlas_queries.npz", P / "atlas_metadata.csv"]
        },
        indent=2,
    )
)
print(
    summary[
        summary.metric.isin(
            ["co_primary", "low", "joint_primary", "joint_all4", "joint_noMYC"]
        )
    ]
    .pivot(index=["study", "state"], columns="metric", values="value")
    .round(3)
    .to_string()
)
print("\nADJUSTED")
print(pd.DataFrame(adj).round(3).to_string(index=False))
