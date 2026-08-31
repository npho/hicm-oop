# ============================================================================
# Figure S5  --  paused ESC (GSE81285) + diapause (GSE241462) analysis
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The sterol programme tracks mTOR activity in mouse paused ESCs and
#   diapause - links the axis to a known mTOR-low state.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   paused ESC (GSE81285) + diapause (GSE241462) analysis
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S5. The sterol programme tracks mTOR activity in mouse paused
#   ESCs and diapause (secondary support). Left: sterol delta-score is lower
#   in mTOR-inhibition-paused ESCs than in active (2i/serum) ESCs across
#   three lines. Right: sterol delta-score is lower in diapaused than
#   reactivated mouse embryos. Both are low-mTOR states, consistent with
#   mTOR acting upstream of the sterol programme.
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cesfw_runs/mouse_sterol_fast.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
DIA="/Users/adelosan/cESFW_analysis/data/diapause"
sterol=["Hmgcr","Hmgcs1","Sqle","Ldlr","Fdft1","Dhcr7","Dhcr24","Mvd","Mvk","Idi1","Cyp51","Lss","Sc5d","Insig1","Fasn","Acaca","Scd1","Scd2","Elovl6","Acly"]
def delta_vec(ld):  # ld: genes x samples log2; returns sample-wise delta (pathway mean - bg mean)
    p=[g for g in sterol if g in ld.index]
    pm=ld.loc[p].mean(axis=0); bgm=ld.mean(axis=0)  # bg≈all genes (pathway is tiny, negligible bias)
    return pm-bgm

print("1. PAUSED MOUSE ESCs (GSE81285): sterol delta-score by line/condition")
df=pd.read_csv(f"{DIA}/GSE81285_RNA-seq_expression_mm10.tsv.gz",sep="\t",index_col=0)
norm=[c for c in df.columns if "Normalized" in c]; df=df[norm]
df.columns=[c.replace(" RNA-seq Normalized Expression","") for c in df.columns]
ld=np.log2(df+1); dv=delta_vec(ld)
low=0
for L in ["E14","Oct4 GiP","v6.5"]:
    v={c:dv[f"{L} {c}"] for c in ["2i","Serum","Paused"]}
    d=v["Paused"]-np.mean([v["2i"],v["Serum"]]); low+=d<0
    print(f"  {L:10s} 2i={v['2i']:.2f} Serum={v['Serum']:.2f} Paused={v['Paused']:.2f}  -> Paused {d:+.2f} {'LOWER' if d<0 else 'higher'}")
print(f"  => sterol LOWER in paused in {low}/3 lines")
for tf in ["Srebf1","Srebf2"]:
    if tf in ld.index:
        av={c:np.mean([ld.loc[tf,f'{L} {c}'] for L in ['E14','Oct4 GiP','v6.5']]) for c in ['2i','Serum','Paused']}
        print(f"  {tf}: 2i={av['2i']:.2f} Serum={av['Serum']:.2f} Paused={av['Paused']:.2f}")

print("\n2. MOUSE DIAPAUSE (GSE241462): diapause vs reactivated")
def load_pb(fn):
    d=pd.read_csv(f"{DIA}/{fn}",sep="\t",usecols=[0,2],names=["gene","count"],header=0)
    pb=d.groupby("gene")["count"].sum(); return np.log2(pb/pb.sum()*1e6+1)
m=pd.DataFrame({"EDG7.5":load_pb("GSM7729308_EDG75_gene_counts.tsv.gz"),
                "EDG9.5":load_pb("GSM7729309_EDG95_gene_counts.tsv.gz"),
                "Rea":load_pb("GSM7729310_Rea_gene_counts.tsv.gz")}).fillna(0)
dv2=delta_vec(m)
print(f"  sterol delta: EDG7.5={dv2['EDG7.5']:+.2f} EDG9.5={dv2['EDG9.5']:+.2f} Rea={dv2['Rea']:+.2f}")
dia=np.mean([dv2['EDG7.5'],dv2['EDG9.5']])
print(f"  => diapause={dia:+.2f} vs reactivated={dv2['Rea']:+.2f}  ({'LOWER in diapause' if dia<dv2['Rea'] else 'higher'})")
for tf in ["Srebf1","Srebf2"]:
    if tf in m.index: print(f"  {tf}: EDG7.5={m.loc[tf,'EDG7.5']:.2f} EDG9.5={m.loc[tf,'EDG9.5']:.2f} Rea={m.loc[tf,'Rea']:.2f}")
print("DONE")
