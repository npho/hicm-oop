# ============================================================================
# Figure S20  --  A-E: per-stage specificity, sens/spec/precision, AUC ranking, panel heatmap
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   Specificity of the ICM signature, and why TFCP2L1 cannot localise it.
#   NOTE: this figure number currently COLLIDES with an existing 'Figure
#   S20' (Stirparo pooling artifact) - see UNRESOLVED.md.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   A-E: per-stage specificity, sens/spec/precision, AUC ranking, panel
#   heatmap
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
# ORIGINAL PATH: /private/tmp/claude-501/-Users-adelosan-Downloads-threshold-ptd/b69a1676-2239-4389-9e30-373e2af378bb/scratchpad/make_figS20.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""Figure S20: specificity of the ICM signature; why TFCP2L1 cannot localise it;
and the proposed NANOG/GATA6/ESRRB staining panel."""
import numpy as np, pandas as pd, sys
import warnings; warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0, "/Users/adelosan/cESFW_analysis")
import _load_data_yanagida as ldy
from sklearn.metrics import roc_auc_score

OUT = "/Users/adelosan/cESFW_analysis/plots/Fig_S20_ICM_signature_specificity.png"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica"],
                     "font.size": 7, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6})

Sample_Info, Counts, cESFW_Genes, Embedding, _, _ = ldy.load_embryo_data()
df = Sample_Info[["Dataset", "Embryo", "Manual_Annotations"]].copy()
lib = Counts.sum(axis=1); lib[lib == 0] = 1

ICM_FAM = ["ICM", "ICM/TE Branch"]
icm = df.Manual_Annotations.isin(ICM_FAM).values
BLASTO_OTHER = ["Early TE", "Mid TE", "Mural TE", "Polar TE", "preIm-Epi",
                "Epi/Hypo Branch", "Hypo"]

STAGES = ["8-Cell", "Morula", "ICM/TE Branch", "ICM", "Epi/Hypo Branch", "preIm-Epi",
          "Hypo", "Early TE", "Mid TE", "Polar TE", "Mural TE"]
SHORT = {"8-Cell": "8-cell", "Morula": "Morula", "ICM/TE Branch": "ICM/TE\nBranch",
         "ICM": "ICM", "Epi/Hypo Branch": "Epi/Hypo\nBranch", "preIm-Epi": "preIm-\nEpi",
         "Hypo": "Hypo", "Early TE": "Early TE", "Mid TE": "Mid TE",
         "Polar TE": "Polar TE", "Mural TE": "Mural TE"}

raw = Counts[["ESRRB", "TFCP2L1", "GATA6", "NANOG", "KLF17", "SOX2"]].reindex(df.index)
E = (raw.ESRRB >= 1).values; T = (raw.TFCP2L1 >= 1).values
G = (raw.GATA6 >= 1).values; N = (raw.NANOG >= 1).values

SIGS = {"ESRRB$^+$GATA6$^+$\n(current)": E & G,
        "ESRRB$^+$TFCP2L1$^+$GATA6$^+$": E & T & G,
        "TFCP2L1$^+$GATA6$^+$": T & G,
        "NANOG$^+$GATA6$^+$ESRRB$^+$\n(proposed panel)": N & G & E}
COLS = ["#4C72B0", "#55A868", "#C44E52", "#8172B2"]

fig = plt.figure(figsize=(11.0, 8.2))
gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 1.0, 1.0], hspace=0.62, wspace=0.26,
                      left=0.075, right=0.985, top=0.945, bottom=0.075)

# ---- A: per-stage % positive ----
axA = fig.add_subplot(gs[0, :])
x = np.arange(len(STAGES)); w = 0.2
for i, (name, mask) in enumerate(SIGS.items()):
    vals = [100 * mask[(df.Manual_Annotations == s).values].mean()
            if (df.Manual_Annotations == s).sum() else np.nan for s in STAGES]
    axA.bar(x + (i - 1.5) * w, vals, w, label=name, color=COLS[i], edgecolor="none")
for xi, s in enumerate(STAGES):
    if s in ICM_FAM:
        axA.axvspan(xi - 0.5, xi + 0.5, color="#FFD54F", alpha=0.22, zorder=0)
axA.set_xticks(x); axA.set_xticklabels([SHORT[s] for s in STAGES], fontsize=6.5)
axA.set_ylabel("% of cells positive"); axA.set_ylim(0, 108)
axA.legend(fontsize=6.2, ncol=4, frameon=False, loc="upper center",
           bbox_to_anchor=(0.5, -0.22), columnspacing=1.6, handlelength=1.3)
axA.set_title("A   The ESRRB$^+$GATA6$^+$ signature is not ICM-specific: morula and early\n"
              "      trophectoderm score as high as ICM (ICM-family shaded)",
              fontsize=7.8, loc="left", pad=8, fontweight="bold")
for sp in ("top", "right"): axA.spines[sp].set_visible(False)

# ---- B: sens/spec/precision ----
axB = fig.add_subplot(gs[1, 0])
other = (~icm) & df.Manual_Annotations.notna().values
met = []
for name, mask in SIGS.items():
    tp = (mask & icm).sum(); fp = (mask & other).sum(); tn = other.sum() - fp
    met.append([100 * tp / icm.sum(), 100 * tn / other.sum(), 100 * tp / max(tp + fp, 1)])
met = np.array(met)
xb = np.arange(3); wb = 0.2
for i, name in enumerate(SIGS):
    axB.bar(xb + (i - 1.5) * wb, met[i], wb, color=COLS[i], edgecolor="none")
    for j in range(3):
        axB.text(xb[j] + (i - 1.5) * wb, met[i, j] + 1.5, f"{met[i,j]:.0f}",
                 ha="center", fontsize=5.6)
axB.set_xticks(xb); axB.set_xticklabels(["Sensitivity\n(ICM detected)",
                                         "Specificity\n(non-ICM excluded)", "Precision"])
axB.set_ylabel("%"); axB.set_ylim(0, 112)
axB.set_title("B   Adding TFCP2L1 does not improve specificity;\n      NANOG-gating does",
              fontsize=7.6, loc="left", fontweight="bold")
for sp in ("top", "right"): axB.spines[sp].set_visible(False)

# ---- C: why TFCP2L1 fails — expression by compartment ----
axC = fig.add_subplot(gs[1, 1])
lc_sel = np.log2(raw.div(lib.reindex(df.index), axis=0) * 1e6 + 1)
COMP = {"ICM-family": ICM_FAM, "Morula": ["Morula"], "preIm-Epi": ["preIm-Epi"],
        "Hypo": ["Hypo"], "TE (all)": ["Early TE", "Mid TE", "Mural TE", "Polar TE"]}
gene_ord = ["ESRRB", "TFCP2L1", "GATA6", "NANOG"]
xc = np.arange(len(COMP))
mk = ["o", "s", "^", "D"]; cc = ["#C44E52", "#4C72B0", "#55A868", "#8172B2"]
for gi, g in enumerate(gene_ord):
    vals = [lc_sel.loc[df.Manual_Annotations.isin(st).values, g].mean() for st in COMP.values()]
    axC.plot(xc, vals, marker=mk[gi], color=cc[gi], lw=1.4, ms=4.5, label=g)
axC.set_xticks(xc); axC.set_xticklabels(list(COMP.keys()), fontsize=6.5)
axC.set_ylabel("mean log$_2$CPM")
axC.legend(fontsize=6.0, frameon=False, ncol=2, loc="lower left", handlelength=1.4)
axC.set_title("C   TFCP2L1 is as high in trophectoderm as in ICM\n      (flat line = no "
              "discriminating power)", fontsize=7.6, loc="left", fontweight="bold")
for sp in ("top", "right"): axC.spines[sp].set_visible(False)

# ---- D: AUC ranking ----
axD = fig.add_subplot(gs[2, 0])
oth_b = df.Manual_Annotations.isin(BLASTO_OTHER).values
sel = icm | oth_b
y = icm[sel].astype(int)
PANEL = ["KLF4", "GATA6", "ESRRB", "SOX2", "SUSD2", "NANOG", "PRDM14", "DPPA3",
         "KLF17", "DNMT3L", "ZFP42", "POU5F1", "TFCP2L1"]
aucs = {}
for g in PANEL:
    if g in Counts.columns:
        v = np.log2(Counts.loc[df.index, g] / lib.reindex(df.index) * 1e6 + 1).values[sel]
        if v.std() > 0:
            aucs[g] = roc_auc_score(y, v)
a = pd.Series(aucs).sort_values()
colors = ["#C44E52" if g == "TFCP2L1" else ("#8172B2" if g in ("NANOG", "GATA6", "ESRRB")
          else "#B0B0B0") for g in a.index]
axD.barh(np.arange(len(a)), a.values, color=colors, edgecolor="none", height=0.72)
axD.axvline(0.5, color="k", ls=":", lw=0.8)
axD.set_yticks(np.arange(len(a))); axD.set_yticklabels(a.index, fontsize=6.3)
axD.set_xlabel("AUC, ICM-family vs rest of blastocyst"); axD.set_xlim(0.4, 1.0)
for i, v in enumerate(a.values):
    axD.text(v + 0.006, i, f"{v:.2f}", va="center", fontsize=5.6)
axD.set_title("D   TFCP2L1 ranks last (AUC 0.57, ~2,800$^{th}$ genome-wide);\n"
              "      purple = proposed staining panel", fontsize=7.6, loc="left",
              fontweight="bold")
for sp in ("top", "right"): axD.spines[sp].set_visible(False)

# ---- E: proposed panel heatmap ----
axE = fig.add_subplot(gs[2, 1])
COMP2 = {"ICM-family": ICM_FAM, "Morula": ["Morula"], "preIm-Epi": ["preIm-Epi"],
         "Epi/Hypo Br": ["Epi/Hypo Branch"], "Hypo": ["Hypo"],
         "TE (all)": ["Early TE", "Mid TE", "Mural TE", "Polar TE"]}
genes5 = ["NANOG", "GATA6", "ESRRB", "KLF17"]
M = np.array([[lc_sel.loc[df.Manual_Annotations.isin(st).values, g].mean()
               for st in COMP2.values()] for g in genes5])
cmap = LinearSegmentedColormap.from_list("bl", ["#F7F7F7", "#9EC9E2", "#2166AC"])
im = axE.imshow(M, cmap=cmap, aspect="auto", vmin=0, vmax=9)
axE.set_xticks(np.arange(len(COMP2))); axE.set_xticklabels(list(COMP2.keys()),
                                                           rotation=35, ha="right", fontsize=6.3)
axE.set_yticks(np.arange(len(genes5))); axE.set_yticklabels(genes5, fontsize=6.8)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        axE.text(j, i, f"{M[i,j]:.1f}", ha="center", va="center", fontsize=5.8,
                 color="white" if M[i, j] > 5.6 else "#333333")
axE.add_patch(plt.Rectangle((-0.5, -0.5), 1, 3, fill=False, ec="#D62728", lw=1.6))
cb = fig.colorbar(im, ax=axE, fraction=0.035, pad=0.02); cb.ax.tick_params(labelsize=5.6)
cb.set_label("mean log$_2$CPM", fontsize=6.0)
axE.set_title("E   Proposed in-situ panel: NANOG$^+$GATA6$^+$ESRRB$^+$ (red box);\n"
              "      KLF17-low separates ICM from morula", fontsize=7.6, loc="left",
              fontweight="bold")

fig.savefig(OUT, dpi=330, bbox_inches="tight", facecolor="white")
print("wrote", OUT)

# receipts CSV
rows = []
for s in STAGES:
    m = (df.Manual_Annotations == s).values
    if not m.sum(): continue
    r = {"stage": s, "n": int(m.sum())}
    for name, mask in SIGS.items():
        r[name.replace("$^+$", "+").replace("\n", " ")] = round(100 * mask[m].mean(), 1)
    rows.append(r)
pd.DataFrame(rows).to_csv("/Users/adelosan/cESFW_analysis/plots/ICM_signature_specificity_byStage.csv",
                          index=False)
print(pd.DataFrame(rows).to_string(index=False))
