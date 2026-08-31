# ============================================================================
# Figure 7  --  C,D: intestinal absorptive/secretory decision
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The transition-versus-state framework generalises to four further human
#   fate decisions, each with a synthetic 50/50 mixture as positive control.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   C,D: intestinal absorptive/secretory decision
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure 7. The transition-versus-state framework generalises across four
#   human fate transitions, from development to regeneration. Order-
#   parameter modality test applied to (A,B) the astrocyte/oligodendrocyte
#   choice of glial progenitors in first-trimester human brain (Braun et
#   al., 2023; thalamus), (C,D) the absorptive/secretory choice of
#   intestinal stem cells in adult small intestine (Burclaff et al., 2022),
#   (E,F) the syncytiotrophoblast/extravillous choice of villous
#   cytotrophoblast in the placenta (Arutyunyan et al., 2023), and (G,H) the
#   AT2-to-AT1 conversion of the KRT8-positive alveolar intermediate in
#   fibrotic adult lung (Habermann et al., 2020). (A,C,E,G) φ = z(fate-A) −
#   z(fate-B): in each system the progenitor or intermediate is a single,
#   unimodal population positioned between its two committed fates.
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/figure_reconstruction/gut_ISC_transition_test.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
Generalisation test #3 (adult human GUT): is the intestinal stem cell (ISC) a
single unresolved TRANSITION between its absorptive and secretory fates, or a
resolved MIXTURE of already-committed cells?  Same order-parameter modality test
used for the ICM (human embryo) and glial progenitors (Braun 2023 brain).

Dataset: UCSC cell browser 'intestine-epithelium' = Burclaff et al. 2022
  (Cell Mol Gastroenterol Hepatol; healthy adult human small intestine + colon
   epithelium, 3 donors, proximal->distal). 12,590 cells, 10x mtx.
Download (once, into scratchpad):
  base=https://cells.ucsc.edu/intestine-epithelium
  curl -O $base/matrix.mtx.gz $base/barcodes.tsv.gz $base/features.tsv.gz $base/meta.tsv
  features are Ensembl IDs -> markers mapped via mygene.info (see gut_markers.json).

METHOD NOTE (important): on this axis the committed secretory (goblet) pole is
extreme, so phi is skewed; raw dBIC(1 vs 2 GMM) is SKEW-sensitive and spuriously
calls the (single-mode) ISC 'bimodal' -- the 2-component fit places both means at
the SAME location (SI: -0.97/-0.92, weights .88/.12). The skew-robust HARTIGAN
DIP TEST is therefore the primary modality test here, with a matched-n synthetic
50/50 (absorptive+goblet) mixture as a positive control.

RESULT (n=150 matched, median over 80 bootstraps):
  Small intestine:  ISC dip p=0.94 (UNIMODAL); synthetic mix dip p=0.000 (multimodal)
  Colon:            ISC dip p=0.94 (UNIMODAL); synthetic mix dip p=0.014 (multimodal)
  means (SI): absorptive +0.80, ISC -0.96 (intermediate, single peak), goblet -4.84
  => The ISC is a single unresolved population between its two fates -- the same
     transition (not mixture) signature as the human ICM and the Braun glial
     progenitor. Reproduced across two intestinal regions.
"""
import gzip, json, subprocess, numpy as np, pandas as pd, diptest
from sklearn.mixture import GaussianMixture
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

SC = "/private/tmp/claude-501/-Users-adelosan/73724e25-f2ba-48fa-a729-bbdde91c906b/scratchpad"
ABS_SYM = ['FABP1','FABP2','APOA1','APOA4','APOB','APOC3','RBP2','ALPI','ANPEP','SLC5A1','MTTP','SI','TREH','SLC2A2']
SEC_SYM = ['MUC2','TFF3','FCGBP','CLCA1','SPINK4','AGR2','ZG16','REG4','SPDEF','ATOH1','TFF1','MUC5AC']

def load():
    mk = json.load(open(f'{SC}/gut_markers.json')); ABS=list(mk['abs'].values()); SEC=list(mk['sec'].values())
    feat=[l.strip() for l in gzip.open(f'{SC}/gut_features.tsv.gz','rt')]
    rowpos={e:i for i,e in enumerate(feat)}
    order={rowpos[e]+1:k for k,e in enumerate(ABS+SEC)}
    ncells=12590; M=np.zeros((len(ABS)+len(SEC),ncells)); colsum=np.zeros(ncells)
    p=subprocess.Popen(['gunzip','-c',f'{SC}/gut_matrix.mtx.gz'],stdout=subprocess.PIPE,text=True,bufsize=1<<20)
    it=p.stdout; [next(it) for _ in range(3)]
    for line in it:
        r,c,v=line.split(); ci=int(c)-1; v=float(v); colsum[ci]+=v
        k=order.get(int(r));
        if k is not None: M[k,ci]=v
    p.wait()
    L=np.log1p(M/np.where(colsum<=0,1,colsum)*1e4)
    m=pd.read_parquet(f'{SC}/gut_meta.parquet')
    return pd.DataFrame(L.T,columns=ABS+SEC), m, ABS, SEC

def main():
    Ld,m,ABS,SEC=load(); Ld.index=m.index
    z=lambda v:(v-np.nanmean(v))/(np.nanstd(v)+1e-9); rng=np.random.RandomState(0)
    def dipboot(x,n=150,r=80):
        x=np.asarray(x); return np.median([diptest.diptest(x if len(x)<=n else x[rng.choice(len(x),n,0)])[1] for _ in range(r)])
    for ctx,name in [('SI_','SMALL INTESTINE'),('C_','COLON')]:
        sub=m[m['lineage'].astype(str).str.startswith(ctx)]; Zc=Ld.loc[sub.index].apply(z)
        phi=pd.Series(z(Zc[ABS].mean(1).values)-z(Zc[SEC].mean(1).values),index=sub.index)
        gg=lambda t: phi[sub['type'].values==t].values
        isc,ab,gob=gg('ISC'),gg('absorptive'),gg('goblet')
        synth=lambda n:np.concatenate([rng.choice(ab,n//2),rng.choice(gob,n-n//2)])
        psyn=np.median([diptest.diptest(synth(150))[1] for _ in range(80)])
        print(f'{name}: ISC dip p={dipboot(isc):.3f} (n={len(isc)}); synthetic-mix dip p={psyn:.3f}; '
              f'means abs={ab.mean():+.2f} ISC={isc.mean():+.2f} goblet={gob.mean():+.2f}')

if __name__=='__main__':
    main()
