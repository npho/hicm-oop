# ============================================================================
# Figure S20  --  unbiased AUC marker ranking
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   Specificity of the ICM signature, and why TFCP2L1 cannot localise it.
#   NOTE: this figure number currently COLLIDES with an existing 'Figure
#   S20' (Stirparo pooling artifact) - see UNRESOLVED.md.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   unbiased AUC marker ranking
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S20. The Stirparo bimodal ICM is a pooling artifact. φ3
#   (multigene epiblast-hypoblast order parameter) distributions for the
#   independent Stirparo annotation. Stirparo Early ICM (n=58; almost all
#   cESFW ICM/TE-branch; mean φ3 −0.17) and ICM (n=63; a later, epiblast-
#   shifted mix; mean +0.25) are each unimodal by Gaussian-mixture BIC (ΔBIC
#   −12.3 and −2.8); the apparent second mode (ΔBIC +11.7) appears only when
#   the two developmentally-offset bins are pooled (right; black curve, the
#   preferred two-component fit). The human ICM does not resolve into two
#   populations. Same φ3 and BIC as the main analysis.
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: authored-here)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /private/tmp/claude-501/-Users-adelosan-Downloads-threshold-ptd/b69a1676-2239-4389-9e30-373e2af378bb/scratchpad/icm_staining_markers.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
What would stain human ICM most specifically?
Unbiased search over all genes: ICM-family vs the compartments an IF experiment must
actually distinguish. Two framings:
  (A) WITHIN A BLASTOCYST (the real experiment): confusables = TE, preIm-Epi/Epi-Hypo, Hypo.
      Morula/8-cell are excluded by morphology, not by antibody.
  (B) ALL STAGES (strict): include Morula/8-cell as confusables too.
Reports AUC + a detection-based precision, and flags nuclear TFs vs surface proteins.
"""
import numpy as np, pandas as pd, sys
import warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, "/Users/adelosan/cESFW_analysis")
import _load_data_yanagida as ldy
from sklearn.metrics import roc_auc_score
from scipy.stats import mannwhitneyu

Sample_Info, Counts, cESFW_Genes, Embedding, _, _ = ldy.load_embryo_data()
df = Sample_Info[["Dataset", "Embryo", "Manual_Annotations"]].copy()
lib = Counts.sum(axis=1); lib[lib == 0] = 1
lc = np.log2(Counts.div(lib, axis=0) * 1e6 + 1).reindex(df.index)

ICM_FAM = ["ICM", "ICM/TE Branch"]
BLASTO_OTHER = ["Early TE", "Mid TE", "Mural TE", "Polar TE", "preIm-Epi",
                "Epi/Hypo Branch", "Hypo"]
ALL_OTHER = BLASTO_OTHER + ["Morula", "8-Cell"]

icm = df.Manual_Annotations.isin(ICM_FAM).values

def search(other_stages, tag, topn=30):
    oth = df.Manual_Annotations.isin(other_stages).values
    sub = lc[icm | oth]
    y = icm[icm | oth].astype(int)
    # keep genes detected in a reasonable fraction so results are stainable, not noise
    det = (sub > 0).mean(axis=0)
    keep = det[(det > 0.05)].index
    X = sub[keep]
    aucs = {}
    for g in keep:
        v = X[g].values
        if v.std() == 0:
            continue
        aucs[g] = roc_auc_score(y, v)
    a = pd.Series(aucs).sort_values(ascending=False)
    res = []
    for g in a.head(topn).index:
        icm_v = lc.loc[icm, g]; oth_v = lc.loc[oth, g]
        res.append({
            "gene": g, "AUC": round(a[g], 3),
            "ICM_mean": round(icm_v.mean(), 2),
            "other_mean": round(oth_v.mean(), 2),
            "ICM_pct_det": round(100 * (icm_v > 0).mean(), 1),
            "other_pct_det": round(100 * (oth_v > 0).mean(), 1),
            "log2FC": round(icm_v.mean() - oth_v.mean(), 2),
        })
    print(f"\n{'='*100}\n### {tag} — top {topn} genes by AUC (ICM-family vs {len(other_stages)} other compartments)\n{'='*100}")
    print(pd.DataFrame(res).to_string(index=False))
    return a

a_blasto = search(BLASTO_OTHER, "A) WITHIN A BLASTOCYST (TE / epiblast / hypoblast as confusables)")
a_all    = search(ALL_OTHER,    "B) ALL STAGES (Morula + 8-cell also counted as confusables)")

# how do the paper's genes rank?
print("\n\n=== Where the paper's genes rank ===")
panel = ["ESRRB", "TFCP2L1", "GATA6", "NANOG", "SOX2", "KLF17", "POU5F1", "OCT4",
         "SUSD2", "DNMT3L", "PDGFRA", "GATA3", "CDX2", "KLF4", "ZFP42", "DPPA3",
         "PRDM14", "TET2", "UHRF1"]
rows = []
for g in panel:
    if g in a_blasto.index:
        rows.append({"gene": g,
                     "AUC_in_blastocyst": round(a_blasto[g], 3),
                     "rank_in_blastocyst": int((a_blasto > a_blasto[g]).sum() + 1),
                     "AUC_all_stages": round(a_all[g], 3) if g in a_all.index else np.nan,
                     "ICM_mean": round(lc.loc[icm, g].mean(), 2)})
print(pd.DataFrame(rows).sort_values("AUC_in_blastocyst", ascending=False).to_string(index=False))

# per-compartment mean for the top candidates, so we can see WHAT each one excludes
print("\n\n=== Per-compartment mean log2CPM for top candidates (what does each exclude?) ===")
COMP = {"ICM-family": ICM_FAM, "Morula": ["Morula"], "8-Cell": ["8-Cell"],
        "preIm-Epi": ["preIm-Epi"], "Epi/Hypo Br": ["Epi/Hypo Branch"], "Hypo": ["Hypo"],
        "Early TE": ["Early TE"], "Mid TE": ["Mid TE"], "Mural TE": ["Mural TE"],
        "Polar TE": ["Polar TE"]}
grp = pd.Series(index=df.index, dtype=object)
for c, stgs in COMP.items():
    grp[df.Manual_Annotations.isin(stgs).values] = c
cands = list(dict.fromkeys(list(a_blasto.head(20).index) + ["ESRRB", "TFCP2L1", "GATA6", "NANOG", "KLF17", "SOX2"]))
cands = [g for g in cands if g in lc.columns]
tbl = lc[cands].groupby(grp).mean().reindex(list(COMP.keys())).round(2)
print(tbl.T.to_string())

# embryo-level reproducibility of the best few (is it driven by one embryo/dataset?)
print("\n\n=== Reproducibility of top candidates across datasets (mean log2CPM in ICM-family) ===")
for g in list(a_blasto.head(10).index) + ["ESRRB"]:
    if g not in lc.columns: continue
    s = lc.loc[icm, g].groupby(df.loc[icm, "Dataset"]).agg(["mean", "count"]).round(2)
    print(f"\n{g}:")
    print(s.to_string())
