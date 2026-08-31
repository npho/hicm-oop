# ============================================================================
# Figure S20  --  underlying specificity statistics
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   Specificity of the ICM signature, and why TFCP2L1 cannot localise it.
#   NOTE: this figure number currently COLLIDES with an existing 'Figure
#   S20' (Stirparo pooling artifact) - see UNRESOLVED.md.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   underlying specificity statistics
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
# ORIGINAL PATH: /private/tmp/claude-501/-Users-adelosan-Downloads-threshold-ptd/b69a1676-2239-4389-9e30-373e2af378bb/scratchpad/icm_signature_specificity.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
Specificity test for candidate ICM signatures in the user's cESFW Manual_Annotations
framework. Question: is ESRRB+GATA6+ ICM-specific, or do pre-ICM (Morula) / TE cells
also score positive? Does adding TFCP2L1 (or swapping ESRRB->TFCP2L1) localise it?
"""
import numpy as np, pandas as pd, sys
import warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, "/Users/adelosan/cESFW_analysis")
import _load_data_yanagida as ldy

Sample_Info, Counts, cESFW_Genes, Embedding, _, _ = ldy.load_embryo_data()

print("=== Manual_Annotations stages available ===")
print(Sample_Info.Manual_Annotations.value_counts().to_string())
print()
print("=== Datasets ===")
print(Sample_Info.Dataset.value_counts().to_string())
print()

GENES = ["ESRRB", "TFCP2L1", "GATA6", "NANOG", "KLF17"]
missing = [g for g in GENES if g not in Counts.columns]
print("missing genes:", missing)
GENES = [g for g in GENES if g in Counts.columns]

# raw molecule counts (normalisation-independent presence/absence), as used for the
# blastoid comparison in the paper
raw = Counts[GENES]

df = Sample_Info[["Dataset", "Embryo", "Manual_Annotations"]].copy()
for g in GENES:
    df[g] = raw[g].reindex(df.index).values

# ordered stage list (developmental order)
ORDER = ["8-Cell", "8-cell", "Morula", "ICM/TE Branch", "ICM", "Epi/Hypo Branch",
         "preIm-Epi", "Hypo", "Early TE", "Medium TE", "Late TE", "TE"]
stages = [s for s in ORDER if s in set(df.Manual_Annotations)]
extra = [s for s in df.Manual_Annotations.unique() if s not in stages and pd.notna(s)]
stages = stages + sorted(extra)

def pos(col, thr=1):
    return df[col] >= thr

for THR in (1, 3):
    print("=" * 100)
    print(f"### DETECTION THRESHOLD >= {THR} molecule(s)  (all datasets pooled) ###")
    print("=" * 100)
    E, T, G = pos("ESRRB", THR), pos("TFCP2L1", THR), pos("GATA6", THR)
    sig = pd.DataFrame({
        "ESRRB+":            E,
        "TFCP2L1+":          T,
        "GATA6+":            G,
        "ESRRB+GATA6+":      E & G,
        "TFCP2L1+GATA6+":    T & G,
        "ESRRB+TFCP2L1+":    E & T,
        "E+T+G+ (triple)":   E & T & G,
    }, index=df.index)

    rows = []
    for s in stages:
        m = (df.Manual_Annotations == s).values
        if m.sum() == 0:
            continue
        r = {"stage": s, "n": int(m.sum())}
        for c in sig.columns:
            r[c] = round(100 * sig.loc[m, c].mean(), 1)
        rows.append(r)
    tab = pd.DataFrame(rows).set_index("stage")
    print(tab.to_string())
    print()

    # ICM-family (user's framework) vs each non-ICM stage: how specific is each signature?
    ICM_FAM = [s for s in ["ICM", "ICM/TE Branch"] if s in stages]
    icm_m = df.Manual_Annotations.isin(ICM_FAM).values
    print(f"ICM-family = {ICM_FAM}, n={icm_m.sum()}")
    print("\n--- Sensitivity / specificity / precision vs ALL other annotated cells ---")
    other_m = (~icm_m) & df.Manual_Annotations.notna().values
    out = []
    for c in sig.columns:
        tp = sig.loc[icm_m, c].sum(); fn = icm_m.sum() - tp
        fp = sig.loc[other_m, c].sum(); tn = other_m.sum() - fp
        out.append({
            "signature": c,
            "sens% (ICM+)": round(100 * tp / max(icm_m.sum(), 1), 1),
            "spec% (non-ICM-)": round(100 * tn / max(other_m.sum(), 1), 1),
            "precision%": round(100 * tp / max(tp + fp, 1), 1),
            "n_FP": int(fp),
        })
    print(pd.DataFrame(out).to_string(index=False))
    print()

    # worst offenders: which non-ICM stage contributes most false positives
    print("--- False-positive breakdown by stage (cells scoring positive that are NOT ICM-family) ---")
    fp_rows = []
    for c in ["ESRRB+GATA6+", "TFCP2L1+GATA6+", "E+T+G+ (triple)"]:
        d = df.loc[sig[c].values & other_m, "Manual_Annotations"].value_counts()
        fp_rows.append(pd.Series(d, name=c))
    print(pd.DataFrame(fp_rows).fillna(0).astype(int).to_string())
    print()

# per-dataset sanity check at thr>=1 for the key contrast Morula vs ICM-family
print("=" * 100)
print("### Per-dataset: Morula vs ICM-family, thr>=1 ###")
print("=" * 100)
E, T, G = pos("ESRRB", 1), pos("TFCP2L1", 1), pos("GATA6", 1)
sig2 = pd.DataFrame({"ESRRB+GATA6+": E & G, "TFCP2L1+GATA6+": T & G,
                     "E+T+G+": E & T & G, "ESRRB+": E, "TFCP2L1+": T}, index=df.index)
for ds in df.Dataset.dropna().unique():
    dm = (df.Dataset == ds).values
    sub = []
    for s in ["Morula", "ICM/TE Branch", "ICM", "preIm-Epi"]:
        m = dm & (df.Manual_Annotations == s).values
        if m.sum() < 3:
            continue
        r = {"stage": s, "n": int(m.sum())}
        for c in sig2.columns:
            r[c] = round(100 * sig2.loc[m, c].mean(), 1)
        sub.append(r)
    if sub:
        print(f"\n--- {ds} ---")
        print(pd.DataFrame(sub).set_index("stage").to_string())

# expression levels (log2 CPM) by stage, for context on "positive at what level"
print("\n" + "=" * 100)
print("### Mean log2CPM by stage (level, not just detection) ###")
print("=" * 100)
lib = Counts.sum(axis=1); lib[lib == 0] = 1
lc = np.log2(Counts[GENES].div(lib, axis=0) * 1e6 + 1)
lvl = lc.reindex(df.index).groupby(df.Manual_Annotations).mean().reindex(
    [s for s in stages if s in set(df.Manual_Annotations)])
print(lvl.round(2).to_string())
