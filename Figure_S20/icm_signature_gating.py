# ============================================================================
# Figure S20  --  candidate gate search
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   Specificity of the ICM signature, and why TFCP2L1 cannot localise it.
#   NOTE: this figure number currently COLLIDES with an existing 'Figure
#   S20' (Stirparo pooling artifact) - see UNRESOLVED.md.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   candidate gate search
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
# ORIGINAL PATH: /private/tmp/claude-501/-Users-adelosan-Downloads-threshold-ptd/b69a1676-2239-4389-9e30-373e2af378bb/scratchpad/icm_signature_gating.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
Follow-up: can any small gene combination localise the ICM signature to ICM-family,
i.e. exclude Morula (upstream) AND trophectoderm (off-lineage)?
Tests detection-based and level-based gates.
"""
import numpy as np, pandas as pd, sys
import warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, "/Users/adelosan/cESFW_analysis")
import _load_data_yanagida as ldy

Sample_Info, Counts, cESFW_Genes, Embedding, _, _ = ldy.load_embryo_data()

GENES = ["ESRRB", "TFCP2L1", "GATA6", "NANOG", "KLF17", "GATA3", "CDX2", "SOX2",
         "PDGFRA", "KLF4", "DPPA3", "ZFP42", "SUSD2", "DNMT3L"]
GENES = [g for g in GENES if g in Counts.columns]
print("genes found:", GENES)

df = Sample_Info[["Dataset", "Embryo", "Manual_Annotations"]].copy()
lib = Counts.sum(axis=1); lib[lib == 0] = 1
lc = np.log2(Counts[GENES].div(lib, axis=0) * 1e6 + 1).reindex(df.index)
raw = Counts[GENES].reindex(df.index)

ICM_FAM = ["ICM", "ICM/TE Branch"]
icm = df.Manual_Annotations.isin(ICM_FAM).values
other = (~icm) & df.Manual_Annotations.notna().values
print(f"\nICM-family n={icm.sum()}, other annotated n={other.sum()}")

# key confusable groups
CONF = {"Morula": ["Morula"], "8-Cell": ["8-Cell"],
        "TE(all)": ["Early TE", "Mid TE", "Mural TE", "Polar TE"],
        "Epi(preIm-Epi+Epi/Hypo)": ["preIm-Epi", "Epi/Hypo Branch"],
        "Hypo": ["Hypo"]}

def P(g, thr=1):   return (raw[g] >= thr).values
def L(g, cut):     return (lc[g] >= cut).values
def N(g, thr=1):   return (raw[g] < thr).values

CANDS = {
    "ESRRB+GATA6+ (current Fig4D)":        P("ESRRB") & P("GATA6"),
    "ESRRB+TFCP2L1+GATA6+ (triple)":       P("ESRRB") & P("TFCP2L1") & P("GATA6"),
    "TFCP2L1+GATA6+":                      P("TFCP2L1") & P("GATA6"),
    "ESRRB+GATA6+NANOG+":                  P("ESRRB") & P("GATA6") & P("NANOG"),
    "ESRRB+GATA6+NANOG+SOX2+":             P("ESRRB") & P("GATA6") & P("NANOG") & P("SOX2"),
    "ESRRB+GATA6+NANOG+ TFCP2L1+":         P("ESRRB") & P("GATA6") & P("NANOG") & P("TFCP2L1"),
}
if "GATA3" in GENES:
    CANDS["ESRRB+GATA6+NANOG+ GATA3-"] = P("ESRRB") & P("GATA6") & P("NANOG") & N("GATA3")
if "CDX2" in GENES:
    CANDS["ESRRB+GATA6+NANOG+ CDX2-"] = P("ESRRB") & P("GATA6") & P("NANOG") & N("CDX2")
if "KLF17" in GENES:
    CANDS["ESRRB+GATA6+NANOG+ KLF17-lo(<4)"] = P("ESRRB") & P("GATA6") & P("NANOG") & (lc["KLF17"] < 4).values

# level-gated versions (log2CPM cutoffs chosen from stage means)
CANDS["ESRRB>=3 & GATA6>=3 (level)"] = L("ESRRB", 3) & L("GATA6", 3)
CANDS["ESRRB>=3 & GATA6>=3 & NANOG>=3"] = L("ESRRB", 3) & L("GATA6", 3) & L("NANOG", 3)
CANDS["ESRRB>=3 & GATA6>=3 & NANOG>=3 & TFCP2L1>=6"] = (
    L("ESRRB", 3) & L("GATA6", 3) & L("NANOG", 3) & L("TFCP2L1", 6))

rows = []
for name, mask in CANDS.items():
    tp = (mask & icm).sum(); fp = (mask & other).sum()
    tn = other.sum() - fp
    r = {"signature": name,
         "sens%": round(100 * tp / icm.sum(), 1),
         "spec%": round(100 * tn / other.sum(), 1),
         "prec%": round(100 * tp / max(tp + fp, 1), 1),
         "F1": round(200 * tp / max(2 * tp + fp + (icm.sum() - tp), 1), 1)}
    for cname, stgs in CONF.items():
        m = df.Manual_Annotations.isin(stgs).values
        r[cname] = round(100 * mask[m].mean(), 1) if m.sum() else np.nan
    rows.append(r)

res = pd.DataFrame(rows).set_index("signature")
print("\n=== Candidate ICM signatures: ICM-family vs all other annotated cells ===")
print("(last columns = % of that confusable group scoring POSITIVE; lower is better)\n")
print(res.to_string())

print("\n\n=== Why TFCP2L1 cannot localise: mean log2CPM by compartment ===")
grp = pd.Series(index=df.index, dtype=object)
for cname, stgs in {**{"ICM-family": ICM_FAM}, **CONF}.items():
    grp[df.Manual_Annotations.isin(stgs).values] = cname
print(lc.groupby(grp).mean().round(2)[["ESRRB", "TFCP2L1", "GATA6", "NANOG", "SOX2", "KLF17"]].to_string())
