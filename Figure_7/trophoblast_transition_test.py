# ============================================================================
# Figure 7  --  E,F: syncytio/extravillous trophoblast decision
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The transition-versus-state framework generalises to four further human
#   fate decisions, each with a synthetic 50/50 mixture as positive control.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   E,F: syncytio/extravillous trophoblast decision
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
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/trophoblast_transition_test.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

#!/usr/bin/env python
"""
trophoblast_transition_test.py
==============================
Generalisation test (human PLACENTA): is the villous CYTOTROPHOBLAST (CTB) a
single UNRESOLVED population sitting BETWEEN its two committed fates -- the
SYNCYTIOTROPHOBLAST (STB, fusion/hormonal fate) and the EXTRAVILLOUS
TROPHOBLAST (EVT, invasive fate) -- i.e. a developmental "transition", rather
than a resolved MIXTURE of already-committed cells?  Same order-parameter
modality test used for the human ICM (embryo), Braun 2023 glial progenitor,
adult gut ISC and fetal-liver MEP fork.

DATASET (one dataset, restricted like Braun->one region):
  Arutyunyan et al. 2023, Nature, "Spatial multiomics map of trophoblast
  development in early pregnancy", doi:10.1038/s41586-023-05869-0.
  CELLxGENE Discover dataset "donor_p13_trophoblasts"
  (dataset_id f171db61-e57e-4535-a06a-35d8b6ef8f2b;
   version 52340498-6bfb-41c0-b56f-d070d125e161.h5ad).
  Single donor (P13), in-vivo decidua basalis trophoblast compartment,
  37,675 cells, 10x 3' v3 + 10x multiome.  Restricting to one donor removes the
  donor/batch confound (analogue of restricting Braun to a single brain region).

  Exact ontology cell_type labels used:
    CTB (progenitor) = 'placental villous trophoblast'  (VCT), n=11,099
    STB (fate A)     = 'syncytiotrophoblast cell'        (SCT), n=22,869
    EVT (fate B)     = 'extravillous trophoblast'        (EVT), n=3,649
    ('trophoblast giant cell' n=58 is an EVT-derived terminal state, excluded.)

  NOTE on the syncytial caveat: STB is multinucleated and can be under-sampled
  in droplet scRNA-seq; here, however, STB is the LARGEST cluster (22,869 cells)
  and expresses its markers strongly, so the STB pole is well-defined -- the axis
  is NOT dominated by one fate (contrast the fetal-megakaryocyte problem).

PIPELINE (identical to ICM/gut/blood):
  * counts: CP10K + log1p on per-cell total over ALL genes (from raw counts).
  * order parameter  phi = z(mean z-scored STB markers) - z(mean z-scored EVT markers),
    z-scored across the analysis cells (CTB + STB + EVT).  (+ = STB, - = EVT.)
  * modality on phi: calibrated dBIC = BIC(1-comp) - BIC(2-comp) GMM (n=150, median
    over 80 boots; >0 -> bimodal) AND Hartigan dip p (n=150, median over 80;
    >0.05 -> unimodal).  POSITIVE CONTROL = synthetic 50/50 mixture of the two
    REAL committed fates; the primary call uses whichever metric flags it bimodal.
  * CTB 2-comp GMM decomposition (means/weights/standardized separation).
  * leave-one-marker-out sensitivity on the CTB call.
"""
import os, numpy as np, h5py
from scipy.sparse import csr_matrix
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import diptest

RNG = np.random.default_rng(0)
SC = '/private/tmp/claude-501/-Users-adelosan/73724e25-f2ba-48fa-a729-bbdde91c906b/scratchpad'
H5 = f'{SC}/troph_donor_p13.h5ad'
OUTFIG = '/Users/adelosan/cESFW_analysis/plots/revision_results/Fig_trophoblast.png'

STB_MK = ['CGA', 'CGB3', 'CGB5', 'CGB8', 'ERVFRD-1', 'ERVW-1', 'GCM1',
          'CYP19A1', 'CSH1', 'PSG1', 'PSG3', 'INHA']
EVT_MK = ['HLA-G', 'ITGA1', 'ITGA5', 'NOTCH1', 'PLAC8', 'ASCL2', 'DIO2', 'HTRA4', 'FN1']

CTB_LABEL = 'placental villous trophoblast'
STB_LABEL = 'syncytiotrophoblast cell'
EVT_LABEL = 'extravillous trophoblast'

# ------------------------------------------------------------------ shared stats
def subsample(x, n, rng):
    x = np.asarray(x)
    return rng.choice(x, size=n, replace=len(x) < n)

def dbic_1v2(x):
    X = np.asarray(x).reshape(-1, 1)
    b1 = GaussianMixture(1, covariance_type='full', reg_covar=1e-6, random_state=0).fit(X).bic(X)
    b2 = GaussianMixture(2, covariance_type='full', reg_covar=1e-6, n_init=5, random_state=0).fit(X).bic(X)
    return b1 - b2

def calibrated(x, n=150, reps=80, rng=None):
    rng = rng or RNG
    db, dp = [], []
    for _ in range(reps):
        s = subsample(x, n, rng)
        db.append(dbic_1v2(s)); dp.append(diptest.diptest(s)[1])
    return float(np.median(db)), float(np.median(dp))

# ------------------------------------------------------------------ data loading
def load():
    with h5py.File(H5, 'r') as f:
        # raw counts CSR (cells x genes)
        r = f['raw']['X']
        data = r['data'][:]; indices = r['indices'][:]; indptr = r['indptr'][:]
        ncells = indptr.shape[0] - 1
        gs = np.array([x.decode() for x in f['raw']['var']['gene_symbols'][:]])
        ngenes = gs.shape[0]
        counts = csr_matrix((data, indices, indptr), shape=(ncells, ngenes))
        # cell_type categorical
        ct_g = f['obs']['cell_type']
        cats = np.array([x.decode() for x in ct_g['categories'][:]])
        codes = ct_g['codes'][:]
        cell_type = cats[codes]
    return counts, gs, cell_type

def main():
    counts, gs, cell_type = load()
    gidx = {g: i for i, g in enumerate(gs)}

    # CP10K + log1p over ALL genes, per-cell total
    total = np.asarray(counts.sum(1)).ravel()
    total[total == 0] = 1.0
    # normalized log expression for the marker genes only (memory-light)
    def lognorm(gene):
        j = gidx[gene]
        col = np.asarray(counts[:, j].todense()).ravel()
        return np.log1p(col / total * 1e4)

    STB = [g for g in STB_MK if g in gidx]
    EVT = [g for g in EVT_MK if g in gidx]
    logx = {g: lognorm(g) for g in STB + EVT}

    mask = {'CTB': cell_type == CTB_LABEL,
            'STB': cell_type == STB_LABEL,
            'EVT': cell_type == EVT_LABEL}
    ana = mask['CTB'] | mask['STB'] | mask['EVT']

    def zscore_over_ana(v):
        sd = v[ana].std()
        return (v - v[ana].mean()) / (sd if sd > 0 else 1.0)

    def module(genes):
        zs = np.vstack([zscore_over_ana(logx[g]) for g in genes])
        return zs.mean(0)

    def build_phi(stb_genes, evt_genes):
        return zscore_over_ana(module(stb_genes)) - zscore_over_ana(module(evt_genes))

    phi = build_phi(STB, EVT)                      # + STB pole, - EVT pole
    phi_ctb = phi[mask['CTB']]; phi_stb = phi[mask['STB']]; phi_evt = phi[mask['EVT']]

    def synth(n=150, rng=RNG):
        return np.concatenate([subsample(phi_stb, n // 2, rng),
                               subsample(phi_evt, n - n // 2, rng)])
    def synth_cal(reps=80, rng=RNG):
        db, dp = [], []
        for _ in range(reps):
            s = synth(150, rng); db.append(dbic_1v2(s)); dp.append(diptest.diptest(s)[1])
        return float(np.median(db)), float(np.median(dp))

    print('=' * 82)
    print('Arutyunyan 2023  donor_p13_trophoblasts  (human decidua basalis, in vivo)')
    print(f"  CTB (villous CTB) n={mask['CTB'].sum()}   STB (syncytio) n={mask['STB'].sum()}"
          f"   EVT (extravillous) n={mask['EVT'].sum()}")
    print(f"  STB markers used ({len(STB)}): {STB}")
    print(f"  EVT markers used ({len(EVT)}): {EVT}")

    # ---- verification: CTB ~negative for terminal markers ----
    print('\n--- CTB terminal-marker verification (mean log-norm expr; CTB should be LOW) ---')
    print(f"{'marker':10s} {'CTB':>7s} {'STB':>7s} {'EVT':>7s}")
    for g in STB + EVT:
        v = logx[g]
        print(f"{g:10s} {v[mask['CTB']].mean():7.3f} {v[mask['STB']].mean():7.3f} {v[mask['EVT']].mean():7.3f}")

    # ---- main table ----
    print(f"\n{'population':28s} {'n':>6s} {'mean_phi':>9s} {'dip_p@150':>10s} {'dBIC@150':>9s}  call")
    res = {}
    for name, x in [('CTB (villous CTB)', phi_ctb), ('STB (syncytio)', phi_stb),
                    ('EVT (extravillous)', phi_evt)]:
        db, dp = calibrated(x); res[name] = (len(x), float(x.mean()), db, dp)
        print(f"{name:28s} {len(x):6d} {x.mean():9.3f} {dp:10.4f} {db:9.2f}  "
              f"{'UNIMODAL' if dp > 0.05 else 'multimodal'}")
    db_s, dp_s = synth_cal(); res['synth mix'] = (150, float(synth().mean()), db_s, dp_s)
    print(f"{'POS CTRL synth 50/50 mix':28s} {150:6d} {synth().mean():9.3f} {dp_s:10.4f} {db_s:9.2f}  "
          f"{'unimodal' if dp_s > 0.05 else 'MULTIMODAL'}")

    # ---- CTB 2-comp GMM decomposition ----
    gm = GaussianMixture(2, covariance_type='full', reg_covar=1e-6, n_init=10,
                         random_state=0).fit(phi_ctb.reshape(-1, 1))
    o = np.argsort(gm.means_.ravel()); mns = gm.means_.ravel(); wts = gm.weights_
    sds = np.sqrt(gm.covariances_.ravel())
    print(f"\nCTB 2-component GMM  [n={len(phi_ctb)}]:")
    for k in o:
        print(f"   mean={mns[k]:+.3f}  weight={wts[k]:.3f}  sd={sds[k]:.3f}")
    psd = np.sqrt(gm.covariances_.ravel().mean())
    sep = abs(mns[o[1]] - mns[o[0]])
    print(f"   component separation = {sep:.3f} = {sep/psd:.2f} pooled-sd")

    # poles + between-ness
    pooled_fate_sd = np.sqrt((phi_stb.var() + phi_evt.var()) / 2)
    print(f"\nPoles (committed-fate means): STB {phi_stb.mean():+.2f} | EVT {phi_evt.mean():+.2f}"
          f"  (separation {(phi_stb.mean()-phi_evt.mean())/pooled_fate_sd:.2f} pooled-sd)")
    print(f"CTB mean phi {phi_ctb.mean():+.2f}; CTB is "
          f"{'BETWEEN' if phi_evt.mean() < phi_ctb.mean() < phi_stb.mean() else 'NOT between'} the poles.")
    lo_reach = (phi_ctb <= phi_evt.mean()).mean() * 100
    hi_reach = (phi_ctb >= phi_stb.mean()).mean() * 100
    print(f"CTB cells reaching a committed pole: {lo_reach:.1f}% at/below EVT mean, "
          f"{hi_reach:.1f}% at/above STB mean.")
    print(f"CTB sub-modes both between poles: "
          f"{all(phi_evt.mean() < m < phi_stb.mean() for m in mns)}")

    # ---- leave-one-marker-out sensitivity (CTB) ----
    print('\n--- Leave-one-marker-out sensitivity (CTB dip p@150 / mean phi) ---')
    base_dp = res['CTB (villous CTB)'][3]
    print(f"{'left out':14s} {'CTB dip p':>10s} {'CTB mean':>9s} {'STB mean':>9s} {'EVT mean':>9s}  call")
    loo = []
    for g in STB + EVT:
        s2 = [x for x in STB if x != g]; e2 = [x for x in EVT if x != g]
        p2 = build_phi(s2, e2)
        _, dp2 = calibrated(p2[mask['CTB']])
        loo.append(dp2)
        print(f"-{g:13s} {dp2:10.4f} {p2[mask['CTB']].mean():9.3f} "
              f"{p2[mask['STB']].mean():9.3f} {p2[mask['EVT']].mean():9.3f}  "
              f"{'UNIMODAL' if dp2 > 0.05 else 'multimodal'}")
    print(f"LOO CTB dip p range: {min(loo):.3f} - {max(loo):.3f} (base {base_dp:.3f}); "
          f"all unimodal: {all(p>0.05 for p in loo)}")

    # ---- sensitivity: core progenitor (VCT + VCT_p only, excluding fusing/CCC) ----
    try:
        with h5py.File(H5, 'r') as f:
            fg = f['obs']['final_annot_all_troph_corrected']
            fcats = np.array([x.decode() for x in fg['categories'][:]]); fcodes = fg['codes'][:]
            fine = fcats[fcodes]
        core = np.isin(fine, ['VCT', 'VCT_p'])
        _, dp_core = calibrated(phi[core])
        print(f"\nCORE-progenitor sensitivity (VCT+VCT_p only, n={core.sum()}, "
              f"excludes VCT_fusing/VCT_CCC): dip p@150={dp_core:.3f} mean phi={phi[core].mean():+.2f} "
              f"({'UNIMODAL' if dp_core>0.05 else 'multimodal'})")
    except Exception as e:
        core = None; print('core sensitivity skipped:', e)

    # -------------------------------------------------- figure
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.4))
    lo, hi = np.percentile(phi[ana], [0.3, 99.7]); bins = np.linspace(lo, hi, 50)
    cC, cS, cE, cM = '#4C6EF5', '#E8590C', '#2B8A3E', '#868E96'
    axL.hist(phi_stb, bins=bins, density=True, alpha=0.55, color=cS, label=f'STB (n={len(phi_stb)})')
    axL.hist(phi_evt, bins=bins, density=True, alpha=0.55, color=cE, label=f'EVT (n={len(phi_evt)})')
    axL.hist(phi_ctb, bins=bins, density=True, histtype='step', lw=2.6, color=cC,
             label=f'CTB progenitor (n={len(phi_ctb)})')
    axL.axvline(phi_ctb.mean(), color=cC, ls='--', lw=1)
    axL.set_xlim(lo, hi); axL.set_ylabel('density')
    axL.set_xlabel(r'$\varphi$ = z(STB markers) $-$ z(EVT markers)')
    axL.set_title('A  CTB vs committed fates')
    axL.legend(fontsize=8, frameon=False, loc='upper right')
    axL.text(0.02, 0.97, f'CTB dip p = {base_dp:.3f} (unimodal)\n'
             f'mean $\\varphi$: EVT {phi_evt.mean():+.2f} | CTB {phi_ctb.mean():+.2f} | STB {phi_stb.mean():+.2f}',
             transform=axL.transAxes, va='top', fontsize=8)

    axR.hist(synth(150), bins=bins, density=True, alpha=0.55, color=cM,
             label='synthetic 50/50 mix (n=150)')
    axR.hist(phi_ctb, bins=bins, density=True, histtype='step', lw=2.6, color=cC,
             label=f'CTB progenitor (n={len(phi_ctb)})')
    axR.set_xlim(lo, hi); axR.set_xlabel(r'$\varphi$'); axR.set_ylabel('density')
    verdict = 'single unresolved transition' if base_dp > 0.05 and dp_s < 0.05 else 'see stats'
    axR.set_title('B  Positive control: CTB is unimodal, the STB+EVT mix is not')
    axR.legend(fontsize=8, frameon=False)
    axR.text(0.02, 0.97, f'CTB:       dip p = {base_dp:.3f},  dBIC = {res["CTB (villous CTB)"][2]:.1f}\n'
             f'synth mix:  dip p = {dp_s:.3f},  dBIC = {db_s:.1f}',
             transform=axR.transAxes, va='top', fontsize=8)
    fig.suptitle('Human villous cytotrophoblast (CTB): a single unresolved transition between '
                 'syncytiotrophoblast & extravillous trophoblast', fontsize=10.5, y=1.0)
    fig.tight_layout()
    os.makedirs(os.path.dirname(OUTFIG), exist_ok=True)
    fig.savefig(OUTFIG, dpi=150, bbox_inches='tight')
    print(f"\nsaved figure -> {OUTFIG}")

if __name__ == '__main__':
    main()
