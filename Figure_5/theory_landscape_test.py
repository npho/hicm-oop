# ============================================================================
# Figure 5  --  A-D: empirical quasi-potential, normal-form fit, criticality test
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   Quasi-potential reconstruction of the ICM->epiblast/hypoblast decision:
#   the ICM sits in a single flattened basin that only later splits in two.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   A-D: empirical quasi-potential, normal-form fit, criticality test
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure 5. A quasi-potential reconstruction of the human
#   ICM→epiblast/hypoblast decision. Order parameter
#   ϕ=z(NANOG)-z(GATA6)along the branching lineage (n=441). (A) Quasi-
#   potential U(ϕ;μ)=-logP(ϕ∣μ)reconstructed across pseudotime (early at
#   bottom, late at top; curves offset for clarity): a single narrow basin
#   broadens and flattens through the ICM. (B) Var(ϕ)peaks at the ICM under
#   this two-gene order parameter (bootstrap peak-at-ICM 63%;
#   a(μ)=1/Var(ϕ)minimal there); this localisation is specific to ϕ1and does
#   not survive the multigene composite (Figure S10). (C) NANOG versus GATA6
#   resolves from early co-expression to mutually exclusive
#   epiblast/hypoblast states. (D) Two-component mixture support (ΔBIC)
#   along pseudotime: unimodal through the ICM in this binning, bimodal only
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: ms+hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/theory_landscape_test.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""EMPIRICAL TEST of the quasi-potential / bifurcation Theory model on real data.
Order parameter phi = z(NANOG) - z(GATA6) along the branching lineage
ICM/TE Branch -> ICM -> Epi/Hypo Branch -> {preIm-Epi (epiblast), Hypo}.
Tests: (1) landscape U(phi;mu) = -log P flattens then splits;
(2) Var(phi) peaks at the transition; (3) the intermediate is UNIMODAL-broad
(true criticality) vs already BIMODAL (trivial epiblast+hypoblast mixture)."""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import gaussian_kde, spearmanr
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _load_data import load_embryo_data

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots", "revision_results")
os.makedirs(OUT, exist_ok=True)
LIN = ["ICM/TE Branch", "ICM", "Epi/Hypo Branch", "preIm-Epi", "Hypo"]; ROOT = "ICM/TE Branch"

SI, C, cG, E, *_ = load_embryo_data()
ann = SI["Manual_Annotations"].astype(str).values
m = np.isin(ann, LIN)
ad = sc.AnnData(X=C.values[m].astype(float)); ad.var_names = C.columns; ad.obs_names = C.index[m]
ad.obs["ann"] = ann[m]
sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
al = pd.DataFrame(ad.X, index=ad.obs_names, columns=ad.var_names)
sc.pp.neighbors(ad, n_neighbors=15, use_rep="X", metric="correlation"); sc.tl.diffmap(ad)
ad.uns["iroot"] = int(np.where(ad.obs["ann"].values == ROOT)[0][0]); sc.tl.dpt(ad)
pt = ad.obs["dpt_pseudotime"].values
n = ad.n_obs
print("lineage cells: %d" % n, "| per stage:", {s: int((ad.obs['ann'].values==s).sum()) for s in LIN})

nan = al["NANOG"].values; gat = al["GATA6"].values
z = lambda v: (v - v.mean()) / (v.std() + 1e-9)
phi = z(nan) - z(gat)          # order parameter: >0 epiblast, <0 hypoblast

# ---- bin pseudotime; per-bin variance, modes, bimodality ----
nb = 8
edges = np.quantile(pt, np.linspace(0, 1, nb + 1)); edges[-1] += 1e-9
grid = np.linspace(phi.min() - 0.5, phi.max() + 0.5, 200)
def n_modes(x):
    if len(x) < 15: return np.nan, None
    kde = gaussian_kde(x); d = kde(grid)
    peaks = [i for i in range(1, len(d)-1) if d[i] > d[i-1] and d[i] > d[i+1] and d[i] > 0.05*d.max()]
    return len(peaks), d
def bimod_gain(x):  # BIC improvement of 2-comp over 1-comp GMM (positive => bimodal preferred)
    if len(x) < 20: return np.nan
    X = x.reshape(-1, 1)
    b1 = GaussianMixture(1, n_init=2, random_state=0).fit(X).bic(X)
    b2 = GaussianMixture(2, n_init=2, random_state=0).fit(X).bic(X)
    return b1 - b2

rows = []; dens = []
for i in range(nb):
    sel = (pt >= edges[i]) & (pt < edges[i+1])
    x = phi[sel]; nm, d = n_modes(x)
    stage = pd.Series(ad.obs["ann"].values[sel]).mode()
    rows.append(dict(bin=i, mid=pt[sel].mean(), nn=int(sel.sum()), var=x.var(),
                     var_nanog=nan[sel].var(), var_gata6=gat[sel].var(),
                     modes=nm, bimod=bimod_gain(x), stage=stage.iloc[0] if len(stage) else "?"))
    dens.append(d)
df = pd.DataFrame(rows)
print("\n=== per pseudotime bin ===")
print(df[["bin","mid","nn","stage","var","modes","bimod"]].round(3).to_string(index=False))

# ---- variance peak: bootstrap CI on which bin has max Var(phi) ----
rng = np.random.RandomState(0)
peakbin = []
for _ in range(500):
    idx = rng.randint(0, n, n)
    vb = [phi[idx][(pt[idx]>=edges[i])&(pt[idx]<edges[i+1])].var() for i in range(nb)]
    peakbin.append(int(np.nanargmax(vb)))
from collections import Counter
pc = Counter(peakbin)
print("\nVar(phi) peak-bin bootstrap distribution:", dict(sorted(pc.items())))
print("modal peak bin = %d (stage %s), in %.0f%% of bootstraps" %
      (max(pc, key=pc.get), df.iloc[max(pc, key=pc.get)]["stage"], 100*pc[max(pc,key=pc.get)]/500))

# ---- criticality vs mixture: is the peak-variance bin UNIMODAL (dip) not bimodal? ----
pk = max(pc, key=pc.get)
print("\n=== is the transition genuinely unimodal-broad (criticality) or bimodal (mixture)? ===")
for i in range(nb):
    tag = "  <-- variance peak" if i == pk else ""
    verdict = "BIMODAL(mixture-like)" if (df.iloc[i]["modes"]>=2 or df.iloc[i]["bimod"]>6) else "unimodal-broad"
    print("  bin %d (%s): modes=%s bimodBIC=%.1f -> %s%s" %
          (i, df.iloc[i]["stage"], df.iloc[i]["modes"], df.iloc[i]["bimod"], verdict, tag))

# ---- FIGURE ----
fig = plt.figure(figsize=(13, 9))
# A: quasi-potential U(phi) waterfall
axA = fig.add_subplot(2, 2, 1)
cmap = plt.cm.viridis(np.linspace(0, 0.92, nb))
for i in range(nb):
    if dens[i] is None: continue
    U = -np.log(dens[i] + 1e-6); U -= U.min()
    axA.plot(grid, U + i*1.2, color=cmap[i], lw=2)
axA.set_title("A  Quasi-potential U($\\phi$;$\\mu$) = -log P along pseudotime\n(bottom=ICM/TE early → top=late; offset for clarity)", fontsize=10)
axA.set_xlabel("$\\phi$ = z(NANOG) - z(GATA6)"); axA.set_ylabel("U + offset (early→late upward)")
axA.spines[["top","right"]].set_visible(False)
# B: variance of phi vs pseudotime
axB = fig.add_subplot(2, 2, 2)
axB.plot(df["mid"], df["var"], "o-", color="#C0392B", lw=2, label="Var($\\phi$)")
axB.axvline(df.iloc[pk]["mid"], color="#888", ls=":", label="variance peak (%s)"%df.iloc[pk]["stage"])
axB.set_title("B  Order-parameter variance peaks at the transition", fontsize=10)
axB.set_xlabel("diffusion pseudotime"); axB.set_ylabel("Var($\\phi$)")
axB.legend(frameon=False, fontsize=8); axB.spines[["top","right"]].set_visible(False)
# C: 2D NANOG-GATA6 at early / mid / late
axC = fig.add_subplot(2, 2, 3)
epochs = [(pt <= edges[2], "#2E86C1", "early (ICM/TE→ICM)"),
          ((pt > edges[3]) & (pt <= edges[5]), "#E8873A", "mid (Epi/Hypo Br.)"),
          (pt > edges[6], "#28B463", "late (Epi/Hypo)")]
for sel, c, lab in epochs:
    axC.scatter(nan[sel], gat[sel], s=12, color=c, alpha=0.5, edgecolors="none", label=lab)
axC.set_title("C  NANOG vs GATA6 resolves from co-expression to exclusive", fontsize=10)
axC.set_xlabel("NANOG (log-norm)"); axC.set_ylabel("GATA6 (log-norm)")
axC.legend(frameon=False, fontsize=8); axC.spines[["top","right"]].set_visible(False)
# D: bimodality (2-comp BIC gain) vs pseudotime
axD = fig.add_subplot(2, 2, 4)
axD.plot(df["mid"], df["bimod"], "s-", color="#6C3483", lw=2)
axD.axhline(6, color="#888", ls=":", label="BIC gain=6 (bimodal threshold)")
axD.set_title("D  Bifurcation: bimodal split emerges only late", fontsize=10)
axD.set_xlabel("diffusion pseudotime"); axD.set_ylabel("2-comp BIC gain (>0 = bimodal)")
axD.legend(frameon=False, fontsize=8); axD.spines[["top","right"]].set_visible(False)
fig.suptitle("Quasi-potential landscape of the human ICM→epiblast/hypoblast bifurcation (n=%d cells)"%n, fontsize=12, y=1.0)
fig.tight_layout()
out = os.path.join(OUT, "Fig_theory_landscape.png"); fig.savefig(out, dpi=200, bbox_inches="tight")
print("\nSaved:", out)
