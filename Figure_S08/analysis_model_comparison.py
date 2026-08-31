# ============================================================================
# Figure S8  --  three-model BIC comparison
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   Transparent model selection: narrow state vs mixture vs broad basin, by
#   BIC over bootstrap resamples.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   three-model BIC comparison
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S8. Transparent model selection for the ICM. For each species,
#   the frequency with which each of three models best explains the ICM
#   order-parameter distribution over 500 bootstrap resamples: a narrow
#   stable state (variance fixed to the mean within-population dispersion of
#   the committed epiblast and hypoblast), a two-component mixture (already-
#   committed cells), and a broad single (near-critical) basin. Results are
#   shown for the canonical two-gene order parameter ϕ1. Under ϕ1the narrow
#   stable-state model is disfavoured in every species (the ICM is 3–9× more
#   dispersed than a resolved lineage and it never wins in human or mouse);
#   this dispersion signature is, however, ϕ1-specific and does not hold
#   under the dropout-insensitive multigene composite, under which the human
#   ICM is no broader than the committed lineages and a narrow-state model
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: ms+hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/analysis_model_comparison.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""#1 Three-model comparison: does the ICM behave as a STABLE ATTRACTOR, a
MIXTURE of committed cells, or a NEAR-CRITICAL (broad unimodal) basin?
Objective model selection by BIC on phi = z(NANOG)-z(GATA6) of the ICM:
  attractor  = 1 Gaussian with variance fixed to the committed-cell dispersion (tight basin)
  mixture    = 2-component Gaussian mixture (already-committed cells)
  near-crit  = 1 Gaussian with free (inflated) variance (broad flat basin)
#2 Per-dataset replication: does the human ICM signature hold in each of the 5
integrated datasets separately (integration robustness)?"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _load_data import load_embryo_data
z=lambda v:(v-v.mean())/(v.std()+1e-9)
def bic_gmm(x,k): X=np.asarray(x).reshape(-1,1); return GaussianMixture(k,n_init=3,random_state=0).fit(X).bic(X)
def bic_attractor(x,var_ref):        # 1 Gaussian, mean MLE, variance fixed to committed dispersion
    x=np.asarray(x); ll=norm.logpdf(x,x.mean(),np.sqrt(var_ref)).sum(); return -2*ll+1*np.log(len(x))
def three_model(icm, committed_vars):
    var_ref=float(np.mean(committed_vars))
    b_att=bic_attractor(icm,var_ref); b_mix=bic_gmm(icm,2); b_cri=bic_gmm(icm,1)
    winner=min([("attractor",b_att),("mixture",b_mix),("near-critical",b_cri)],key=lambda t:t[1])[0]
    return dict(n=len(icm),var_icm=np.var(icm),var_ref=var_ref,
                BIC_attractor=b_att,BIC_mixture=b_mix,BIC_nearcrit=b_cri,winner=winner)

# ================= load all three species =================
LIN=["ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Hypo"]
SI,C,cG,E,*_=load_embryo_data(); ann=SI["Manual_Annotations"].astype(str).values; m=np.isin(ann,LIN)
adH=sc.AnnData(X=C.values[m].astype(float)); adH.var_names=C.columns; adH.obs_names=C.index[m]
adH.obs["a"]=ann[m]; adH.obs["ds"]=SI["Dataset"].astype(str).values[m]
adH.obs["stir"]=SI["Stirparo_Labels"].astype(str).values[m]
sc.pp.normalize_total(adH,target_sum=1e4); sc.pp.log1p(adH); alH=pd.DataFrame(adH.X,index=adH.obs_names,columns=adH.var_names)
phiH=z(alH["NANOG"].values)-z(alH["GATA6"].values); aH=adH.obs["a"].values
icmH=phiH[aH=="ICM"]; commH=[np.var(phiH[aH=="preIm-Epi"]),np.var(phiH[aH=="Hypo"])]

# marmoset
D="data/nhp/marmoset/RawCounts.csv"; meta=pd.read_csv(D,index_col=0,nrows=2); lin=meta.iloc[1]
cnt=pd.read_csv(D,skiprows=[1,2],index_col=0); cnt=cnt[~cnt.index.duplicated()]
cells=[c for c in cnt.columns if lin.get(c) in ("ICM_CS3","Epi_CS3","Hyp_CS3")]
adM=sc.AnnData(X=cnt[cells].T.values.astype(float)); adM.var_names=cnt.index; adM.obs_names=cells; adM.obs["l"]=[lin[c] for c in cells]
sc.pp.normalize_total(adM,target_sum=1e4); sc.pp.log1p(adM); alM=pd.DataFrame(adM.X,index=adM.obs_names,columns=adM.var_names)
phiM=z(alM["NANOG"].values)-z(alM["GATA6"].values); lM=np.array(adM.obs["l"])
icmM=phiM[lM=="ICM_CS3"]; commM=[np.var(phiM[lM=="Epi_CS3"]),np.var(phiM[lM=="Hyp_CS3"])]

# mouse
DATA="data/mouse_dev"; F={"E3.5":["GSM3494334_Lib1-3_E3.5_counts.csv.gz","GSM3494335_Lib1-4_E3.5_counts.csv.gz"],
                          "E4.5":["GSM3494336_Lib1-1_E4.5_counts.csv.gz","GSM3494337_Lib1-2_E4.5_counts.csv.gz"]}
fr,st=[],[]
for s,fs in F.items():
    for f in fs:
        d=pd.read_csv(os.path.join(DATA,f),index_col=0); d.index=[f"{f[:10]}_{i}" for i in range(len(d))]; fr.append(d); st+=[s]*len(d)
mt=pd.concat(fr).fillna(0); adK=sc.AnnData(X=mt.values.astype(float)); adK.var_names=[g.upper() for g in mt.columns]
adK.var_names_make_unique(); adK.obs_names=mt.index; adK.obs["stage"]=st
sc.pp.filter_cells(adK,min_genes=1500); sc.pp.normalize_total(adK,target_sum=1e4); sc.pp.log1p(adK)
alK=pd.DataFrame(adK.X,index=adK.obs_names,columns=adK.var_names)
sc.tl.score_genes(adK,[g for g in ["CDX2","GATA3","KRT8","KRT18","TEAD4","EOMES"] if g in adK.var_names],score_name="TE")
phiK=z(alK["NANOG"].values)-z(alK["GATA6"].values); innK=(adK.obs["TE"].values<0.15)
icmK=phiK[innK&(adK.obs["stage"].values=="E3.5")]
e45=phiK[innK&(adK.obs["stage"].values=="E4.5")]; commK=[np.var(e45[e45>0]),np.var(e45[e45<0])]  # resolved poles

print("========== #1 THREE-MODEL COMPARISON (BIC; lowest wins) ==========")
print("%-10s %5s %8s %8s | %10s %10s %11s | %s"%("species","n","var_ICM","var_ref","BIC_attr","BIC_mix","BIC_nearcr","WINNER"))
for nm,icm,cv in [("Human",icmH,commH),("Marmoset",icmM,commM),("Mouse",icmK,commK)]:
    r=three_model(icm,cv)
    print("%-10s %5d %8.2f %8.2f | %10.1f %10.1f %11.1f | %s"%(nm,r["n"],r["var_icm"],r["var_ref"],
          r["BIC_attractor"],r["BIC_mixture"],r["BIC_nearcrit"],r["winner"].upper()))
print("(attractor rejected when var_ICM >> var_ref; mixture rejected when unimodal; near-critical = broad unimodal)")

print("\n========== #2 PER-DATASET REPLICATION (human ICM, by integrated study) ==========")
def dbic(x): X=np.asarray(x).reshape(-1,1); return bic_gmm(x,1)-bic_gmm(x,2)
print("Manual_Annotations=='ICM' split by Dataset:")
for ds in sorted(set(adH.obs["ds"].values)):
    x=phiH[(aH=="ICM")&(adH.obs["ds"].values==ds)]
    if len(x)>=12:
        print("  %-18s n=%3d  var=%.2f  DeltaBIC=%6.2f -> %s"%(ds,len(x),np.var(x),dbic(x),
              "unimodal(near-critical)" if dbic(x)<0 else "bimodal"))
    else:
        print("  %-18s n=%3d  (too few ICM cells to test)"%(ds,len(x)))
print("\nCross-check with INDEPENDENT annotation (Stirparo_Labels 'ICM'+'Early ICM'):")
xs=phiH[np.isin(adH.obs["stir"].values,["ICM","Early ICM"])]
print("  Stirparo ICM: n=%d  var=%.2f  DeltaBIC=%.2f -> %s"%(len(xs),np.var(xs),dbic(xs),
      "unimodal(near-critical)" if dbic(xs)<0 else "bimodal"))

# ---- figure: bootstrap model-selection win-frequency (transparent evidence strength) ----
import matplotlib.pyplot as plt
from collections import Counter as _Ctr
def _boot(icm,var_ref,B=500):
    wins=_Ctr(); rng=np.random.RandomState(0)
    for _ in range(B):
        x=np.asarray(icm)[rng.randint(0,len(icm),len(icm))]
        b=[("attractor",bic_attractor(x,var_ref)),("mixture",bic_gmm(x,2)),("near-critical",bic_gmm(x,1))]
        wins[min(b,key=lambda t:t[1])[0]]+=1
    return {k:wins.get(k,0)/B for k in ["attractor","mixture","near-critical"]}
spf=[("Human",icmH,commH),("Marmoset",icmM,commM),("Mouse",icmK,commK)]
freqs={nm:_boot(icm,float(np.mean(cv))) for nm,icm,cv in spf}
models=[("narrow stable state","attractor","#8a8a8a"),("mixture (resolved into two)","mixture","#2E86C1"),
        ("broad single basin (near-critical)","near-critical","#C0392B")]
x=np.arange(3); w=0.26
fig,ax=plt.subplots(figsize=(8.8,4.7))
for j,(lab,key,c) in enumerate(models):
    vals=[freqs[nm][key]*100 for nm,_,_ in spf]
    ax.bar(x+(j-1)*w,vals,w,color=c,label=lab,edgecolor="white")
    for k,(nm,_,_) in enumerate(spf):
        v=freqs[nm][key]*100
        if v>3: ax.text(x[k]+(j-1)*w,v+1.5,f"{v:.0f}",ha="center",fontsize=8,color=c)
ax.set_xticks(x); ax.set_xticklabels([nm for nm,_,_ in spf]); ax.set_ylim(0,108)
ax.set_ylabel("bootstrap win-frequency (%, 500 resamples)")
ax.set_title("Which model best describes the ICM? Transparent model selection\n"
             "narrow stable state disfavoured everywhere; mouse = resolved mixture; human ambiguous ($\\approx$50/50)",fontsize=10)
ax.legend(frameon=False,fontsize=8.5,loc="upper center"); ax.spines[["top","right"]].set_visible(False)
out="plots/revision_results/Fig_model_comparison.png"; fig.savefig(out,dpi=200,bbox_inches="tight"); print("Saved:",out)

# ---- bootstrap model-selection win frequencies (evidence strength) ----
from collections import Counter
def boot_winfreq(icm, var_ref, B=500):
    wins=Counter(); rng=np.random.RandomState(0)
    for _ in range(B):
        x=np.asarray(icm)[rng.randint(0,len(icm),len(icm))]
        b=[("attractor",bic_attractor(x,var_ref)),("mixture",bic_gmm(x,2)),("near-critical",bic_gmm(x,1))]
        wins[min(b,key=lambda t:t[1])[0]]+=1
    return {k:round(wins.get(k,0)/B,3) for k in ["attractor","mixture","near-critical"]}
print("\n===== BOOTSTRAP win-frequency (500 resamples) =====")
for nm,icm,cv in [("Human",icmH,commH),("Marmoset",icmM,commM),("Mouse",icmK,commK)]:
    print(f"  {nm:9s}: {boot_winfreq(icm, float(np.mean(cv)))}")
