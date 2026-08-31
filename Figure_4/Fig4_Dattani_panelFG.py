# ============================================================================
# Figure 4  --  F,G: Dattani naive-to-hypoblast timecourse
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   No in-vitro system stabilises the ICM programme, but rare genuine ESRRB+
#   cells arise. Two independent measurements (reference mapping; direct
#   signature scoring).
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   F,G: Dattani naive-to-hypoblast timecourse
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure 4. No in-vitro system builds the ICM-like state, but rare genuine
#   ESRRB+ cells arise. (A, B) cESFW reference-mapping of naive/reset PSCs
#   and blastoids onto the human embryo atlas (cosine k-NN on cESFW genes;
#   leave-one-out ICM recall 55%); Naive PSCs map predominantly to the
#   preimplantation epiblast (or, for e4CL/8CLC, the 8-cell/morula stage),
#   whereas blastoids map to epiblast and trophectoderm. The embryo profiled
#   alongside serves as a positive control, with 22.8% ICM. (C) Reference-
#   mapped (Measurement 1): all naive-PSC and blastoid systems contain 2.5%
#   or less ICM versus 22.8% for the embryo. (D) Direct scoring of the dual-
#   competent ICM signature (Measurement 2): fraction of cells co-expressing
#   ESRRB+GATA6+ (solid) and ESRRB+TFCP2L1+GATA6+ (hatched). 93% of cESFW
#   embryo ICM-family cells versus under 5% in every steady-state in-vitro
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/Fig4_Dattani_panelFG.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""Figure 4 panels F,G,H: a directed naive->hypoblast protocol (Dattani et al.
2024) transiently and PARTIALLY reaches the ICM-like state. Three panels
mirroring C/D/E: dual-competent signature (F), cESFW reference-mapping (G),
ICM Identity Score (H). Rendered at Figure-4 width for compositing below A-E.
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE, GREEN, VIOLET = "#2a78d6", "#eb6834", "#008300", "#4a3aa7"
INK, INK_2, INK_MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8880", "#e4e3de", "#fcfcfb"
D = pd.read_csv("plots/Fig_Dattani_ICMlike_data.csv")
S = pd.read_csv("plots/dattani_icmscore_bycond.csv")
cell = pd.read_csv("plots/dattani_icm_similarity.csv")
lab = ["PXGL 0h\n(naïve)", "PDA83 24h\n(intermediate)", "FAX 36h", "FAX 48h", "FAX 72h\n(hypoblast)"]
x = np.arange(5)
AUTH, THRESH = 0.703, 0.287
order = ["PXGL / 0h", "PDA83 / 24h", "PDA83 FAX / 36h", "PDA83 FAX / 48h", "PDA83 FAX / 72h"]

fig, (axF, axG, axH) = plt.subplots(1, 3, figsize=(15.6, 3.7), dpi=300,
                                    gridspec_kw=dict(width_ratios=[1, 1.15, 1], wspace=0.28))
fig.patch.set_facecolor(SURFACE)

# ── F: transient ESRRB+GATA6+ ignition ───────────────────────────────────────
axF.set_facecolor(SURFACE)
axF.bar(x, D["EG"], width=0.62, color=[ORANGE if i == 1 else INK_MUTED for i in range(5)], zorder=3)
for i, v in enumerate(D["EG"]):
    axF.text(i, v + 0.4, f"{v:.1f}%", ha="center", fontsize=8,
             color=INK, fontweight="semibold" if i == 1 else "normal")
axF.set_ylim(0, 21); axF.set_ylabel("cells ESRRB⁺GATA6⁺ (%)", fontsize=9, color=INK_2)
axF.set_title("F   Dual-competent signature ignites\n     transiently (Dattani et al. 2024)",
              loc="left", fontsize=10.4, fontweight="semibold", color=INK, pad=7)

# ── G: cESFW reference-mapping composition ───────────────────────────────────
LIN = ["ICM", "ICM/TE Branch", "preIm-Epi", "Epi/Hypo Branch", "Hypo", "Mid TE", "Mural TE", "Polar TE", "8-Cell"]
COL = {"ICM": ORANGE, "ICM/TE Branch": "#f2a07a", "preIm-Epi": BLUE, "Epi/Hypo Branch": "#9085e9",
       "Hypo": GREEN, "Mid TE": "#cfc9bd", "Mural TE": "#b8b2a6", "Polar TE": "#a49e92", "8-Cell": "#dcd7cc"}
comp = (cell[cell["cond"].isin(order)].groupby(["cond", "cESFW_pre"]).size()
        .unstack(fill_value=0).reindex(order))
comp = comp.div(comp.sum(1), axis=0) * 100
axG.set_facecolor(SURFACE)
bottom = np.zeros(5)
for lin in LIN:
    if lin not in comp.columns: continue
    v = comp[lin].values
    axG.bar(x, v, bottom=bottom, width=0.66, color=COL.get(lin, "#ccc"), label=lin,
            edgecolor=SURFACE, linewidth=0.4, zorder=3)
    bottom += v
axG.set_ylim(0, 100); axG.set_ylabel("cESFW identity (% of cells)", fontsize=9, color=INK_2)
axG.set_title("G   cESFW mapping: preimplantation EPIBLAST →\n     HYPOBLAST, essentially no ICM (<0.5%)",
              loc="left", fontsize=10.4, fontweight="semibold", color=INK, pad=7)
# legend only for identities that actually appear, placed below panel G
present_lin = [l for l in LIN if l in comp.columns and comp[l].max() > 0.3]
handles = [plt.Rectangle((0, 0), 1, 1, color=COL[l]) for l in present_lin]
axG.legend(handles, present_lin, frameon=False, fontsize=6.8, labelcolor=INK_2,
           ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.16),
           handlelength=1.1, columnspacing=1.0, handletextpad=0.4)

# ── H: ICM Identity Score vs authentic ICM ───────────────────────────────────
axH.set_facecolor(SURFACE)
axH.axhline(AUTH, color=ORANGE, lw=1.5, ls=(0, (4, 3)), zorder=2)
axH.text(2.0, AUTH + 0.02, "authentic embryo ICM (+0.70)", ha="center", fontsize=7.4, color=ORANGE, va="bottom")
axH.axhline(THRESH, color=INK_MUTED, lw=1, ls=(0, (2, 2)), zorder=2)
axH.text(2.0, THRESH + 0.012, "authentic-ICM threshold (+0.29)", ha="center", fontsize=6.9, color=INK_MUTED, va="bottom")
axH.plot(x, S["ICM_score"], "-o", color=INK, lw=1.9, ms=6, zorder=5)
axH.scatter([1], [S["ICM_score"][1]], s=130, color=ORANGE, zorder=6, edgecolor=SURFACE, lw=1.2)
for i, v in enumerate(S["ICM_score"]):
    axH.text(i, v + 0.035, f"{v:+.2f}", ha="center", fontsize=8,
             color=INK, fontweight="semibold" if i == 1 else "normal")
axH.set_ylim(0, AUTH + 0.16); axH.set_ylabel("ICM Identity Score (40-gene)", fontsize=8.8, color=INK_2)
axH.set_title("H   ...but ~4× below authentic ICM (only 5%\n     of cells reach its range)",
              loc="left", fontsize=10.4, fontweight="semibold", color=INK, pad=7)

for ax in (axF, axG, axH):
    ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=7.6, color=INK_2)
    ax.grid(True, axis="y", color=GRID, lw=0.6, zorder=0); ax.set_axisbelow(True)
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)

fig.subplots_adjust(left=0.05, right=0.985, top=0.83, bottom=0.28)
fig.savefig("plots/Fig4_Dattani_panelFG.png", facecolor=SURFACE, dpi=300)
print("wrote plots/Fig4_Dattani_panelFG.png (F,G,H)")
