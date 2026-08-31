# ============================================================================
# Figure S20  --  the OTHER Figure S20 (Stirparo pooling artifact)
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   Specificity of the ICM signature, and why TFCP2L1 cannot localise it.
#   NOTE: this figure number currently COLLIDES with an existing 'Figure
#   S20' (Stirparo pooling artifact) - see UNRESOLVED.md.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   the OTHER Figure S20 (Stirparo pooling artifact)
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
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/Fig_TableS2_Stirparo_pooling.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""Table S2 companion: the Stirparo 'bimodal' ICM is a pooling artifact.

phi3 (multigene epiblast-hypoblast order parameter) distributions for the
Stirparo independent annotation: 'Early ICM' and 'ICM' are each unimodal, and
only their pooling looks bimodal, because the two annotation bins sit at
different points on the epiblast-hypoblast axis. Same phi3 and GMM/BIC as
pilot_multigene_phi.py.

Output: plots/revision_results/Fig_TableS2_Stirparo_pooling.{png,pdf}
"""
import os, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from sklearn.mixture import GaussianMixture
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from _load_data import load_embryo_data

BLUE, ORANGE, INK, INK_2, INK_MUTED = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#8a8880"
SURFACE, GRID = "#fcfcfb", "#e4e3de"
z = lambda v: (v - v.mean()) / (v.std() + 1e-9)
EPI = ["NANOG","SOX2","POU5F1","KLF17","TDGF1","DPPA5","PRDM14"]
HYP = ["GATA6","GATA4","SOX17","PDGFRA","FOXA2","HNF1B"]

def dbic(x):
    X = np.asarray(x).reshape(-1, 1)
    return (GaussianMixture(1, n_init=3, random_state=0).fit(X).bic(X)
            - GaussianMixture(2, n_init=3, random_state=0).fit(X).bic(X))

LIN = ["ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Hypo"]
SI, C, cG, E, *_ = load_embryo_data()
ann = SI["Manual_Annotations"].astype(str).values
m = np.isin(ann, LIN)
ad = sc.AnnData(X=C.values[m].astype(float)); ad.var_names = C.columns; ad.obs_names = C.index[m]
ad.obs["stir"] = SI["Stirparo_Labels"].astype(str).values[m]
sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
sc.tl.score_genes(ad, [g for g in EPI if g in ad.var_names], score_name="ep")
sc.tl.score_genes(ad, [g for g in HYP if g in ad.var_names], score_name="hy")
phi3 = z(ad.obs["ep"].values) - z(ad.obs["hy"].values)
stir = ad.obs["stir"].values
early = phi3[stir == "Early ICM"]
icm = phi3[stir == "ICM"]
pooled = phi3[np.isin(stir, ["ICM", "Early ICM"])]

fig, ax = plt.subplots(1, 3, figsize=(11.4, 3.7), dpi=300, sharex=True, sharey=True)
fig.patch.set_facecolor(SURFACE)
bins = np.linspace(-3.5, 3.5, 22)
xs = np.linspace(-3.5, 3.5, 400)

bw = bins[1] - bins[0]

def gmm_curve(x, k):
    g = GaussianMixture(k, n_init=3, random_state=0).fit(np.asarray(x).reshape(-1, 1))
    return np.exp(g.score_samples(xs.reshape(-1, 1))) * len(x) * bw  # scaled to counts

def panel(a, data, title, color, k_fit):
    a.set_facecolor(SURFACE)
    a.hist(data, bins=bins, color=color, alpha=0.5,
           edgecolor=SURFACE, linewidth=0.6, zorder=3)
    a.plot(xs, gmm_curve(data, k_fit), color=INK, lw=1.6, zorder=5)
    db = dbic(data)
    verdict = "bimodal" if db > 0 else "unimodal"
    a.set_title(f"{title}\nΔBIC {db:+.1f} → {verdict}", fontsize=9.5,
                color=INK, fontweight="semibold", linespacing=1.5)
    a.axvline(0, color=INK_MUTED, lw=0.8, ls=(0, (2, 2)), zorder=1)
    a.grid(True, color=GRID, lw=0.6, zorder=0); a.set_axisbelow(True)
    a.tick_params(length=0, labelsize=8, colors=INK_2)
    for s in a.spines.values(): s.set_visible(False)

panel(ax[0], early, f"Stirparo “Early ICM”  (n={len(early)})", BLUE, 1)
panel(ax[1], icm, f"Stirparo “ICM”  (n={len(icm)})", ORANGE, 1)

# panel C: pooled, two bins stacked so the "two modes" are visibly the two bins;
# black curve = the 2-component GMM that BIC prefers
ax[2].set_facecolor(SURFACE)
ax[2].hist([early, icm], bins=bins, stacked=True,
           color=[BLUE, ORANGE], alpha=0.55, edgecolor=SURFACE, linewidth=0.6,
           label=["Early ICM", "ICM"], zorder=3)
ax[2].plot(xs, gmm_curve(pooled, 2), color=INK, lw=1.6, zorder=5)
dbp = dbic(pooled)
ax[2].set_title(f"Pooled ICM + Early ICM  (n={len(pooled)})\nΔBIC {dbp:+.1f} → bimodal",
                fontsize=9.5, color=INK, fontweight="semibold", linespacing=1.5)
ax[2].axvline(0, color=INK_MUTED, lw=0.8, ls=(0, (2, 2)), zorder=1)
ax[2].grid(True, color=GRID, lw=0.6, zorder=0); ax[2].set_axisbelow(True)
ax[2].tick_params(length=0, labelsize=8, colors=INK_2)
for s in ax[2].spines.values(): s.set_visible(False)
ax[2].legend(frameon=False, fontsize=7.5, labelcolor=INK_2, loc="upper right",
             handlelength=1.1, borderpad=0.2)

ax[0].set_ylabel("cells", fontsize=8.5, color=INK_2)
for a in ax:
    a.set_xlabel("φ3  (epiblast − hypoblast score)", fontsize=8.5, color=INK_2)

fig.suptitle("The Stirparo “bimodal” ICM is a pooling artifact — each annotation bin is itself unimodal",
             x=0.012, y=1.03, ha="left", fontsize=11.5, fontweight="semibold", color=INK)
fig.text(0.012, -0.06,
         "φ3 (multigene epiblast−hypoblast order parameter) for the independent Stirparo annotation. “Early ICM” (≈ cESFW ICM/TE-branch, mean φ3 −0.17) and “ICM” "
         "(a later, epiblast-shifted mix, mean +0.25) are each\nunimodal by Gaussian-mixture BIC; the apparent second mode appears only when the two developmentally-offset bins are pooled. "
         "The ICM is not two resolved populations. Same φ3/BIC as the main analysis (n as shown).",
         ha="left", va="top", fontsize=6.6, color=INK_MUTED, linespacing=1.6)

fig.subplots_adjust(left=0.06, right=0.99, top=0.80, bottom=0.30, wspace=0.12)
for ext in ("png", "pdf"):
    p = f"/Users/adelosan/cESFW_analysis/plots/revision_results/Fig_TableS2_Stirparo_pooling.{ext}"
    fig.savefig(p, facecolor=SURFACE, bbox_inches="tight")
    print("wrote", p)
print(f"Early ICM ΔBIC={dbic(early):+.1f}  ICM ΔBIC={dbic(icm):+.1f}  pooled ΔBIC={dbp:+.1f}")
