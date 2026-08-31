# ============================================================================
# Figure 5  --  bifurcation diagram (fixed points vs pseudotime)
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   Quasi-potential reconstruction of the ICM->epiblast/hypoblast decision:
#   the ICM sits in a single flattened basin that only later splits in two.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   bifurcation diagram (fixed points vs pseudotime)
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
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: ms)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/theory_bifurcation.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""Formal bifurcation analysis (#1) + cross-species generality (#3).
HUMAN: fit the canonical normal form U(phi)=c - h*phi + (a/2)phi^2 + (b/4)phi^4 per
pseudotime window; extract stability a(mu), draw the bifurcation diagram (stable
fixed points vs pseudotime), confirm variance peak + unimodality at the critical point.
CYNO: repeat the variance/stability test on the cynomolgus ICM->epi/hypo decision.
phi = z(NANOG) - z(GATA6)."""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import gaussian_kde
from sklearn.mixture import GaussianMixture
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _load_data import load_embryo_data
np.seterr(all="ignore")

def order_param_pt(al, pt):
    z = lambda v: (v - v.mean()) / (v.std() + 1e-9)
    return z(al["NANOG"].values) - z(al["GATA6"].values), pt

def fit_normalform(phi):
    """linear LSQ of U=-log density on basis [1,phi,phi^2,phi^4]; return a,b,h."""
    if len(phi) < 20: return np.nan, np.nan, np.nan
    g = np.linspace(np.percentile(phi,2)-0.3, np.percentile(phi,98)+0.3, 120)
    d = gaussian_kde(phi)(g); U = -np.log(d + 1e-6)
    w = d / d.sum()                                   # fit where cells are
    Xb = np.vstack([np.ones_like(g), g, g**2, g**4]).T
    W = np.diag(w)
    beta = np.linalg.lstsq(W@Xb, W@U, rcond=None)[0]
    a = 2*beta[2]; b = 4*beta[3]; h = -beta[1]
    return a, b, h

def stable_fixed_points(a, b, h):
    if not np.isfinite([a,b,h]).all() or b <= 1e-6:
        return [ -a/1.0 ] if False else []          # degenerate
    roots = np.roots([b, 0, a, -h])
    real = roots[np.abs(roots.imag) < 1e-6].real
    return sorted([r for r in real if (a + 3*b*r**2) > 0])   # minima (stable)

def modes(phi, grid):
    if len(phi) < 12: return []
    d = gaussian_kde(phi)(grid)
    return [grid[i] for i in range(1,len(d)-1) if d[i]>d[i-1] and d[i]>d[i+1] and d[i]>0.05*d.max()]

def bimod_bic(phi):
    if len(phi) < 20: return np.nan
    X = phi.reshape(-1,1)
    return (GaussianMixture(1,n_init=2,random_state=0).fit(X).bic(X)
            - GaussianMixture(2,n_init=2,random_state=0).fit(X).bic(X))

def analyse(al, pt, label, nb=8):
    phi, pt = order_param_pt(al, pt)
    n = len(phi); edges = np.quantile(pt, np.linspace(0,1,nb+1)); edges[-1]+=1e-9
    grid = np.linspace(phi.min()-0.5, phi.max()+0.5, 200)
    rows = []
    for i in range(nb):
        sel = (pt>=edges[i]) & (pt<edges[i+1]); x = phi[sel]
        a,b,h = fit_normalform(x)
        rows.append(dict(mid=pt[sel].mean(), nn=int(sel.sum()), var=x.var(),
                         a=a, b=b, h=h, nmodes=len(modes(x,grid)), bic=bimod_bic(x)))
    df = pd.DataFrame(rows)
    pk = int(np.nanargmax(df["var"].values))
    print(f"\n===== {label} (n={n}) =====")
    print(df[["mid","nn","var","a","nmodes","bic"]].round(3).to_string(index=False))
    _amsg = (f"min a(mu)={np.nanmin(df.a):.3f} at bin {int(np.nanargmin(df.a))}"
             if np.isfinite(df.a.values).any() else "a(mu) fit unavailable (too few cells/bin)")
    print(f"  variance peak: bin {pk} (pt={df.mid[pk]:.3f}); a(mu) at peak={df.a[pk]:.3f}; {_amsg}")
    print(f"  unimodal at variance peak? nmodes={df.nmodes[pk]}, bic={df.bic[pk]:.2f} "
          f"({'UNIMODAL (mixture rejected)' if df.bic[pk]<0 else 'bimodal'})")
    print(f"  bimodal (split) first at bin: "
          f"{next((i for i in range(nb) if df.nmodes[i]>=2 and df.bic[i]>0), 'none')}")
    return df, phi, pt

# ---------------- HUMAN ----------------
LIN = ["ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Hypo"]; ROOT="ICM/TE Branch"
SI,C,cG,E,*_ = load_embryo_data()
ann = SI["Manual_Annotations"].astype(str).values; m = np.isin(ann, LIN)
ad = sc.AnnData(X=C.values[m].astype(float)); ad.var_names=C.columns; ad.obs_names=C.index[m]; ad.obs["ann"]=ann[m]
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
alH = pd.DataFrame(ad.X, index=ad.obs_names, columns=ad.var_names)
sc.pp.neighbors(ad,n_neighbors=15,use_rep="X",metric="correlation"); sc.tl.diffmap(ad)
ad.uns["iroot"]=int(np.where(ad.obs["ann"].values==ROOT)[0][0]); sc.tl.dpt(ad)
dfH, phiH, ptH = analyse(alH, ad.obs["dpt_pseudotime"].values, "HUMAN")

# sliding-window bifurcation diagram (finer)
def bifurcation_curve(phi, pt, win=70, step=8):
    o = np.argsort(pt); phi, pt = phi[o], pt[o]; grid=np.linspace(phi.min()-.5,phi.max()+.5,200)
    xs, brs = [], []
    for s in range(0, len(phi)-win, step):
        seg = phi[s:s+win]; xs.append(pt[s:s+win].mean()); brs.append(modes(seg, grid))
    return xs, brs
xsH, brsH = bifurcation_curve(phiH, ptH)
nsplit = sum(1 for b in brsH if len(b)>=2)
print(f"  [human bifurcation diagram] {len(xsH)} windows; {nsplit} show 2 branches (split).")

# ---------------- CYNOMOLGUS ----------------
print("\n---------------- CYNO ----------------")
Dc = "data/Nakamura2016_Cyno"
raw = pd.read_csv(os.path.join(Dc,"Cy_Processed.txt.gz"),sep="\t").set_index("macFas5_gene_symbol")
cells=[c for c in raw.columns if c!="macFas5_entrez_id"]
X = raw[cells].apply(pd.to_numeric,errors="coerce").fillna(0); X=X[~X.index.duplicated()]
lab = pd.read_csv("primate_cesfw/cyno_embedding.csv",index_col=0)["lineage"]
keep = [c for c in cells if c in lab.index and lab[c] in ("ICM/early","Epiblast","Hypoblast")]
adc = sc.AnnData(X=X[keep].T.values.astype(float)); adc.var_names=X.index; adc.obs_names=keep
adc.obs["lin"]=lab[keep].values
sc.pp.normalize_total(adc,target_sum=1e4); sc.pp.log1p(adc)
alC = pd.DataFrame(adc.X, index=adc.obs_names, columns=adc.var_names)
sc.pp.neighbors(adc,n_neighbors=min(12,adc.n_obs-1),use_rep="X",metric="correlation"); sc.tl.diffmap(adc)
adc.uns["iroot"]=int(np.where(adc.obs["lin"].values=="ICM/early")[0][0]); sc.tl.dpt(adc)
print("cyno lineage counts:", {k:int((adc.obs['lin'].values==k).sum()) for k in ['ICM/early','Epiblast','Hypoblast']})
for g in ["NANOG","GATA6"]:
    print(f"  {g} present in cyno: {g in alC.columns}")
dfC, phiC, ptC = analyse(alC, adc.obs["dpt_pseudotime"].values, "CYNO", nb=5)

# save arrays for figure
np.savez("/private/tmp/claude-501/-Users-adelosan/73724e25-f2ba-48fa-a729-bbdde91c906b/scratchpad/v24/theory_bif.npz",
         phiH=phiH, ptH=ptH, dfH=dfH.values, dfH_cols=dfH.columns.values,
         xsH=np.array(xsH,dtype=object), brsH=np.array(brsH,dtype=object),
         phiC=phiC, ptC=ptC, dfC=dfC.values, dfC_cols=dfC.columns.values, allow_pickle=True)
print("\nsaved arrays for figure.")
