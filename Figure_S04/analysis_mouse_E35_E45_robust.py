# ============================================================================
# Figure S4  --  covariate-adjusted robustness
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The mTOR-sterol programme also rises across normal mouse E3.5 -> E4.5
#   (cross-species support that the axis is not a human artefact).
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   covariate-adjusted robustness
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S4. The mTOR–sterol programme rises across normal mouse
#   E3.5 → E4.5 (secondary, cross-species support). Sterol/SREBF (top) and
#   mTORC1 (bottom) delta-scores from the E3.5 ICM (common origin) into the
#   E4.5 epiblast and primitive-endoderm (PrE) fates, in Mohammed 2017
#   (Smart-seq2) and Nowotschin 2019 (10x). The conserved sterol increase is
#   carried by the PrE (hypoblast) arm; the epiblast-specific sterol rise is
#   weak and not robust across both datasets, so this supports the overall
#   ICM maturation rather than the epiblast lineage specifically.
#
# ATTRIBUTION CONFIDENCE: MEDIUM   (evidence: name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/analysis_mouse_E35_E45_robust.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
ROBUST cross-species test: does sterol/SREBF (and mTORC1) rise across normal
mouse E3.5 (ICM) -> E4.5 (epiblast/PrE), after the same controls used for the
human claims — cell-cycle (S/G2M) and sequencing depth adjustment?

Two independent single-cell datasets (Mohammed 2017 Smart-seq2; Nowotschin 2019
10x). E3.5 ICM = common origin (denominator); TE excluded. Lineage by markers
(validated). Test = OLS  score ~ isE4.5 + S_score + G2M_score + log10_counts;
report the adjusted E4.5 coefficient. "Holds" = sterol E4.5 coefficient > 0 with
p<0.05 in BOTH datasets after adjustment.

Run: conda run -n cesfw_env python analysis_mouse_E35_E45_robust.py
"""
import os, re, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import scanpy as sc
import statsmodels.formula.api as smf

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "mouse_dev")
OUT = os.path.join(HERE, "plots", "revision_results", "mouse_E35_E45")
os.makedirs(OUT, exist_ok=True)

U = str.upper
STEROL = [U(x) for x in ["Hmgcs1","Hmgcr","Mvk","Pmvk","Mvd","Idi1","Fdps","Fdft1","Sqle","Lss",
          "Cyp51","Tm7sf2","Msmo1","Nsdhl","Hsd17b7","Sc5d","Dhcr7","Dhcr24","Ldlr","Insig1",
          "Srebf1","Srebf2","Scap","Fasn","Acat2"]]
MTORC1 = [U(x) for x in ["Slc7a5","Slc2a1","Hk2","Ldha","Eif4ebp1","Rps6kb1","Ddit4","Hspa5",
          "Aldoa","Pgk1","Pkm","Slc1a5","Cad","Psat1","Shmt2"]]
EPI = [U(x) for x in ["Nanog","Sox2","Pou5f1","Fgf4","Tdgf1","Klf4","Esrrb","Zfp42"]]
PRE = [U(x) for x in ["Gata6","Gata4","Sox17","Pdgfra","Fn1","Col4a1","Lamb1","Foxa2","Dab2"]]
TE  = [U(x) for x in ["Cdx2","Gata3","Krt8","Krt18","Tead4","Eomes","Elf5","Krt7"]]
S_G  = [U(x) for x in ["Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6",
        "Cdca7","Dtl","Prim1","Uhrf1","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76",
        "Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1",
        "Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8"]]
G2M = [U(x) for x in ["Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2",
        "Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb",
        "Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Cdc20","Ttk","Cdc25c",
        "Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka",
        "Psrc1","Anln","Lbr","Ckap5","Cenpe","Nek2","G2e3","Cbx5","Cenpa"]]


def prep(ad):
    ad.var_names = [U(v) for v in ad.var_names]; ad.var_names_make_unique()
    ad.obs["log10_counts"] = np.log10(np.asarray(ad.X.sum(1)).ravel() + 1)
    sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
    sc.tl.score_genes_cell_cycle(ad, s_genes=[g for g in S_G if g in ad.var_names],
                                 g2m_genes=[g for g in G2M if g in ad.var_names])
    for prog, nm in [(STEROL,"sterol"),(MTORC1,"mtorc1"),(EPI,"epi"),(PRE,"pre"),(TE,"te")]:
        sc.tl.score_genes(ad, [g for g in prog if g in ad.var_names], score_name=nm, ctrl_size=50)
    ad.obs["lineage"] = np.array(["Epiblast","PrE","TE"])[ad.obs[["epi","pre","te"]].values.argmax(1)]
    return ad


def load_mohammed():
    raw = pd.read_csv(os.path.join(DATA,"GSE100597_Mohammed.txt.gz"), sep="\t", index_col=0)
    raw.index = [i.rsplit("_",2)[0] for i in raw.index]; raw = raw[~raw.index.duplicated()]
    st = pd.Series([re.search(r"E[0-9]\.[0-9]+",c).group(0) for c in raw.columns], index=raw.columns)
    m = st.isin(["E3.5","E4.5"]); raw=raw.loc[:,m]; st=st[m]
    ad = sc.AnnData(X=raw.T.values.astype(float)); ad.var_names=raw.index; ad.obs_names=raw.columns
    ad.obs["stage"]=st.values
    return prep(ad)


def load_nowotschin():
    F={"E3.5":["GSM3494334_Lib1-3_E3.5_counts.csv.gz","GSM3494335_Lib1-4_E3.5_counts.csv.gz"],
       "E4.5":["GSM3494336_Lib1-1_E4.5_counts.csv.gz","GSM3494337_Lib1-2_E4.5_counts.csv.gz"]}
    fr,stg=[],[]
    for s,fs in F.items():
        for f in fs:
            d=pd.read_csv(os.path.join(DATA,f),index_col=0); d.index=[f"{f[:10]}_{i}" for i in range(len(d))]
            fr.append(d); stg+=[s]*len(d)
    mat=pd.concat(fr).fillna(0)
    ad=sc.AnnData(X=mat.values.astype(float)); ad.var_names=mat.columns; ad.obs_names=mat.index; ad.obs["stage"]=stg
    sc.pp.filter_cells(ad,min_genes=1500); sc.pp.filter_genes(ad,min_cells=3)
    return prep(ad)


def adjusted(ad, prog, tag):
    d = ad.obs[ad.obs.lineage!="TE"].copy()          # ICM-derived only
    res=[]
    for glabel, sel in [("E4.5 all", d.stage=="E4.5"),
                        ("E4.5 Epi", (d.stage=="E4.5")&(d.lineage=="Epiblast")),
                        ("E4.5 PrE", (d.stage=="E4.5")&(d.lineage=="PrE"))]:
        sub = d[(d.stage=="E3.5") | sel].copy()
        sub["isE45"]=(sub.stage=="E4.5").astype(int)
        if sub.isE45.sum()<4:
            res.append((glabel,np.nan,np.nan,np.nan)); continue
        # raw diff
        raw_d = sub[sub.isE45==1][prog].median()-sub[sub.isE45==0][prog].median()
        model = smf.ols(f"{prog} ~ isE45 + S_score + G2M_score + log10_counts", data=sub).fit()
        beta, p = model.params["isE45"], model.pvalues["isE45"]
        res.append((glabel, raw_d, beta, p))
        print(f"  [{tag}] {prog:6s} E3.5 ICM -> {glabel:9s}: raw Δmed={raw_d:+.3f} | "
              f"CC+depth-adj β={beta:+.3f} | p={p:.1e}")
    return res


print("Loading + scoring (cell-cycle, depth, programs)…")
DS = {"Mohammed 2017": load_mohammed(), "Nowotschin 2019": load_nowotschin()}

# validate lineage calls
print("\n=== Lineage-marker validation (mean expr; Epi should be Nanog-high, PrE Gata6-high) ===")
for name, ad in DS.items():
    d=ad.obs
    print(f"{name}: median epi-score  Epi={d[d.lineage=='Epiblast'].epi.median():+.2f} "
          f"PrE={d[d.lineage=='PrE'].epi.median():+.2f} | pre-score Epi={d[d.lineage=='Epiblast'].pre.median():+.2f} "
          f"PrE={d[d.lineage=='PrE'].pre.median():+.2f}")
    print(f"   lineage x stage: {pd.crosstab(d.lineage,d.stage).to_dict()}")

print("\n=== CELL-CYCLE + DEPTH ADJUSTED tests (E3.5 ICM = reference) ===")
rows=[]
for name, ad in DS.items():
    for prog in ["sterol","mtorc1"]:
        for glabel, raw_d, beta, p in adjusted(ad, prog, name):
            rows.append(dict(dataset=name, program=prog, contrast=glabel,
                             raw_dmedian=raw_d, adj_beta=beta, adj_p=p))
summary=pd.DataFrame(rows)
summary.to_csv(os.path.join(OUT,"mouse_E35_E45_ADJUSTED_summary.csv"), index=False)

# verdict
def holds(prog, contrast):
    s=summary[(summary.program==prog)&(summary.contrast==contrast)]
    return (s.adj_beta>0).all() and (s.adj_p<0.05).all()
print("\n=== VERDICT (survives CC+depth adjustment in BOTH datasets) ===")
for prog in ["sterol","mtorc1"]:
    for con in ["E4.5 all","E4.5 Epi","E4.5 PrE"]:
        print(f"  {prog:6s} E3.5 ICM -> {con:9s}: {'HOLDS' if holds(prog,con) else 'does NOT hold'}")
print("\nSaved ->", os.path.join(OUT,"mouse_E35_E45_ADJUSTED_summary.csv"))
