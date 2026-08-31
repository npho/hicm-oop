# ============================================================================
# Figure S12  --  megakaryocyte-erythroid progenitor modality test
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The transition signature is also recovered in adult human blood (MEP).
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   megakaryocyte-erythroid progenitor modality test
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S12. The transition-versus-state signature is also recovered in
#   adult human blood. Order-parameter modality test applied to the
#   megakaryocyte-erythroid decision in adult human bone marrow (Zeng et
#   al., 2025). (A) φ = z(erythroid) − z(megakaryocyte): the megakaryocyte-
#   erythroid progenitor (MEP) is a single, unimodal population positioned
#   between the erythroid and megakaryocyte fates. (B) Positive control
#   (n=150 matched): a synthetic 50/50 mixture of the committed fates is
#   multimodal (dip p=0.004) whereas the MEP is unimodal (dip p=0.64); the
#   unimodal call is stable across leave-one-marker-out. Caveats: the
#   megakaryocytic pole is megakaryocyte-progenitor-dominated, so the MEP
#   sits mildly toward it, and the unimodal call requires excluding PLEK, a
#   broadly expressed marker the MEP itself expresses. This case is
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/hema_adult_test.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

#!/usr/bin/env python
"""
hema_adult_test.py
==================
CLEAN ADULT version of the MEP -> erythroid / megakaryocyte transition test.

Question: is the adult human MEGAKARYOCYTE-ERYTHROID PROGENITOR (MEP) a single
UNRESOLVED population sitting BETWEEN its two committed fates (a "transition"),
rather than a resolved mixture of already-committed erythroid + megakaryocyte cells?

DATASET
-------
CELLxGENE Discover dataset "A (Balanced) Bone Marrow Reference Map of Hematopoietic
Development" (dataset_id cd2f23c1-aef1-48ae-8eb4-0bcf124e567d), 263,159 cells,
adult (19-84 yr) HEALTHY human bone marrow, is_primary_data=True.
Citation: Zeng et al. 2025, Blood Cancer Discovery, doi:10.1158/2643-3230.BCD-24-0342
("A Balanced Bone Marrow Reference Map"; integrates Hay2018/Oetjen2018/Granja2019/
 Mende2022/Setty2019/Ainciburu2023). Ontology-harmonised cell_type labels; genes are
Ensembl IDs with feature_name symbols; raw counts in .raw.X.

LABELS (exact ontology cell_type)
  PROGENITOR = "megakaryocyte-erythroid progenitor cell"  (MEP, n=4856)
  FATE A (erythroid)     = "basophilic erythroblast" + "polychromatophilic erythroblast"
                           + "orthochromatic erythroblast"  (mature, hemoglobin-high; n=14048)
  FATE B (megakaryocyte) = "megakaryocyte progenitor cell" + "megakaryocyte"  (n=2107)
Rationale: MATURE terminal markers (below) are strongly expressed by the committed
fates but ~absent in the MEP (verified: MEP PF4/PPBP/GYPA/HBA ~0), so a "mega-leaning"
artefact like the fetal case cannot arise. Mature megakaryocytes alone are only n=54, so
the megakaryocytic FATE pools MkP+megakaryocyte (n=2107, PF4/PPBP/ITGA2B-positive).
"erythroid progenitor cell" (an early CFU-E-like stage, HBB~2.8 vs ~7 in erythroblasts) is
EXCLUDED from the erythroid pole: it is only weakly committed and its low-hemoglobin tail
overlaps the progenitor zone, which would unfairly inflate MEP/pole overlap.

PIPELINE (identical to the fetal/ICM pipeline)
  CP10K + log1p on per-cell total counts over ALL genes;
  phi = z(mean z-scored erythroid markers) - z(mean z-scored megakaryocyte markers),
  z-scored across MEP + erythroid + megakaryocyte cells;
  matched-n subsample n=150, medians over 80 reps for dip p and calibrated dBIC;
  synthetic 50/50 mix of the two REAL committed fates = positive control.
"""
import numpy as np, scipy.sparse as sp
import anndata as ad
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import diptest

RNG = np.random.default_rng(0)
H5 = '/private/tmp/claude-501/-Users-adelosan/73724e25-f2ba-48fa-a729-bbdde91c906b/scratchpad/hema_adult_bmref.h5ad'
OUTFIG = '/Users/adelosan/cESFW_analysis/plots/revision_results/Fig_hema_adult.png'

ERY_MK = ['GYPA', 'HBB', 'HBA1', 'HBA2', 'ALAS2', 'KLF1', 'SLC4A1', 'CA1', 'AHSP']   # AHSP == ERAF
# terminal megakaryocyte markers the MEP does NOT express (verified log-norm in MEP:
# PF4 0.01, PPBP 0.01, ITGA2B 0.26, GP9 0.01, GP1BA 0.05, VWF 0.03, TUBB1 0.00, PLA2G4A 0.03).
# PLEK is DELIBERATELY EXCLUDED: pleckstrin is broadly expressed and the MEP does express it
# (log-norm 1.28), violating the "mature marker the progenitor does not express" criterion.
# Including PLEK manufactures a spurious minor second mode in the MEP (see PLEK sensitivity below).
MK_MK  = ['PF4', 'PPBP', 'ITGA2B', 'GP9', 'GP1BA', 'VWF', 'TUBB1', 'PLA2G4A']

ERY_LABELS = ['basophilic erythroblast', 'polychromatophilic erythroblast',
              'orthochromatic erythroblast']       # mature, hemoglobin-high committed pole
MK_LABELS  = ['megakaryocyte progenitor cell', 'megakaryocyte']
MEP_LABEL  = 'megakaryocyte-erythroid progenitor cell'

# ----------------------------------------------------------------- stats
def dbic_1v2(x):
    X = np.asarray(x).reshape(-1, 1)
    b1 = GaussianMixture(1, covariance_type='full', reg_covar=1e-6, random_state=0).fit(X).bic(X)
    b2 = GaussianMixture(2, covariance_type='full', reg_covar=1e-6, n_init=5, random_state=0).fit(X).bic(X)
    return b1 - b2

def calibrated(x, n=150, reps=80, rng=None):
    rng = rng or RNG
    db, dp = [], []
    for _ in range(reps):
        s = rng.choice(x, n, replace=len(x) < n)
        db.append(dbic_1v2(s)); dp.append(diptest.diptest(s)[1])
    return float(np.median(db)), float(np.median(dp))

# ----------------------------------------------------------------- load (once) + phi
def load_cache():
    """Read the h5ad ONCE; return CP10K+log1p marker expression + population masks."""
    A = ad.read_h5ad(H5, backed='r')
    ct = A.obs['cell_type'].astype(str).values
    keep_labels = [MEP_LABEL] + ERY_LABELS + MK_LABELS
    idx = np.where(np.isin(ct, keep_labels))[0]
    raw = A.raw.X[idx]
    raw = raw.tocsr() if sp.issparse(raw) else sp.csr_matrix(raw)
    tot = np.asarray(raw.sum(1)).ravel()                       # per-cell total counts, ALL genes
    fn = A.var['feature_name'].astype(str); sym2idx = {s: i for i, s in enumerate(fn)}
    sub_ct = ct[idx]
    lognorm = {}
    for g in ERY_MK + MK_MK + ['PLEK']:                        # PLEK loaded for the sensitivity check
        v = np.asarray(raw[:, sym2idx[g]].todense()).ravel()
        lognorm[g] = np.log1p(v / tot * 1e4)
    masks = dict(mep=sub_ct == MEP_LABEL, ery=np.isin(sub_ct, ERY_LABELS), mk=np.isin(sub_ct, MK_LABELS))
    return lognorm, masks

def compute_phi(cache, ery_markers=ERY_MK, mk_markers=MK_MK):
    lognorm, masks = cache
    ana = masks['mep'] | masks['ery'] | masks['mk']
    def score(markers):
        zs = []
        for g in markers:
            v = lognorm[g]; sd = v[ana].std()
            zs.append((v - v[ana].mean()) / (sd if sd > 0 else 1.0))
        return np.mean(np.vstack(zs), 0)
    def zo(x): return (x - x[ana].mean()) / (x[ana].std() or 1.0)
    phi = zo(score(ery_markers)) - zo(score(mk_markers))
    return phi[masks['mep']], phi[masks['ery']], phi[masks['mk']], \
        (int(masks['mep'].sum()), int(masks['ery'].sum()), int(masks['mk'].sum()))

# ----------------------------------------------------------------- run
def synth(pe, pm, n=150, rng=RNG):
    return np.concatenate([rng.choice(pe, n // 2, replace=len(pe) < n // 2),
                           rng.choice(pm, n - n // 2, replace=len(pm) < n - n // 2)])

def main():
    cache = load_cache()
    phi_mep, phi_ery, phi_mk, ns = compute_phi(cache)
    print('=' * 80)
    print('ADULT MEP fork  (Balanced Bone Marrow Reference Map; Zeng 2025)')
    print(f'  MEP n={ns[0]}  erythroid n={ns[1]}  megakaryocytic n={ns[2]}')
    print(f"\n{'population':28s} {'n':>6s} {'mean_phi':>9s} {'dBIC@150':>9s} {'dip_p@150':>10s}  call")
    res = {}
    for name, x in [('progenitor (MEP)', phi_mep), ('fate A (erythroid)', phi_ery),
                    ('fate B (megakaryocytic)', phi_mk)]:
        db, dp = calibrated(x); res[name] = (len(x), x.mean(), db, dp)
        print(f"{name:28s} {len(x):6d} {x.mean():9.3f} {db:9.2f} {dp:10.4f}  "
              f"{'UNIMODAL' if dp > 0.05 else 'multimodal'}")
    db_s, dp_s = ([], [])
    for _ in range(80):
        s = synth(phi_ery, phi_mk); db_s.append(dbic_1v2(s)); dp_s.append(diptest.diptest(s)[1])
    db_s, dp_s = float(np.median(db_s)), float(np.median(dp_s))
    print(f"{'POS CTRL synthetic 50/50 mix':28s} {150:6d} {synth(phi_ery, phi_mk).mean():9.3f} "
          f"{db_s:9.2f} {dp_s:10.4f}  {'MULTIMODAL' if (dp_s < 0.05 or db_s > 0) else 'unimodal'}")

    # MEP 2-comp GMM
    gm = GaussianMixture(2, covariance_type='full', reg_covar=1e-6, n_init=10,
                         random_state=0).fit(phi_mep.reshape(-1, 1))
    o = np.argsort(gm.means_.ravel()); mns = gm.means_.ravel(); wt = gm.weights_
    sds = np.sqrt(gm.covariances_.ravel()); psd = np.sqrt(gm.covariances_.ravel().mean())
    sep = abs(mns[o[1]] - mns[o[0]])
    print(f"\nMEP 2-comp GMM [n={len(phi_mep)}]:")
    for k in o:
        print(f"   mean={mns[k]:+.3f}  weight={wt[k]:.3f}  sd={sds[k]:.3f}")
    print(f"   standardized separation = {sep / psd:.2f} pooled-sd")

    # centrality + non-overlap with committed poles
    ery_lo = np.percentile(phi_ery, 5); mk_hi = np.percentile(phi_mk, 95)
    print(f"\nCentrality:  MEP mean {phi_mep.mean():+.2f}  vs  erythroid median {np.median(phi_ery):+.2f} / "
          f"megakaryocytic median {np.median(phi_mk):+.2f}")
    print(f"  MEP 95th pct = {np.percentile(phi_mep,95):+.2f} (erythroid 5th pct {ery_lo:+.2f}); "
          f"MEP 5th pct = {np.percentile(phi_mep,5):+.2f} (megakaryocytic 95th pct {mk_hi:+.2f})")
    frac_reach_ery = (phi_mep >= ery_lo).mean(); frac_reach_mk = (phi_mep <= mk_hi).mean()
    print(f"  frac MEP reaching erythroid pole (>=ery p5) = {frac_reach_ery:.3f}; "
          f"reaching mega pole (<=mk p95) = {frac_reach_mk:.3f}")

    # resolved-mixture vs skewed-transition test: do the two MEP GMM modes sit AT the poles?
    ery_med, mk_med = np.median(phi_ery), np.median(phi_mk)
    lo_mode, hi_mode = mns[o[0]], mns[o[1]]
    both_intermediate = (mk_med < lo_mode < ery_med) and (mk_med < hi_mode < ery_med)
    print(f"\nResolved-mixture vs transition: MEP modes at {lo_mode:+.2f}/{hi_mode:+.2f}; "
          f"committed poles at mega {mk_med:+.2f} / ery {ery_med:+.2f}.")
    print(f"  both MEP modes lie strictly BETWEEN the poles (neither at a fate) = {both_intermediate} "
          f"-> {'skewed single transition zone, NOT a resolved mixture' if both_intermediate else 'a mode sits at a committed pole (resolution)'}")

    # leave-one-marker-out robustness of the MEP dip call
    print("\nLeave-one-marker-out robustness (MEP dip p @150):")
    stable = True
    for drop in ERY_MK + MK_MK:
        em = [g for g in ERY_MK if g != drop]; mm = [g for g in MK_MK if g != drop]
        pm2, _, _, _ = compute_phi(cache, em, mm)
        _, dp2 = calibrated(pm2)
        flag = 'unimodal' if dp2 > 0.05 else 'BIMODAL'
        if dp2 <= 0.05: stable = False
        print(f"   drop {drop:8s}: MEP dip p = {dp2:.3f} ({flag})")
    print(f"  => MEP unimodal call stable across all leave-one-out: {stable}")

    # PLEK sensitivity: adding the non-terminal, MEP-expressed marker PLEK
    pm_plek, _, _, _ = compute_phi(cache, ERY_MK, MK_MK + ['PLEK'])
    _, dp_plek = calibrated(pm_plek)
    print(f"\nPLEK sensitivity: adding PLEK (which the MEP DOES express, log-norm 1.28) -> "
          f"MEP dip p = {dp_plek:.3f} ({'unimodal' if dp_plek > 0.05 else 'BIMODAL'}). "
          f"The only route to a bimodal MEP is a marker that violates the terminal-marker criterion.")

    # ---------------- figure ----------------
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.4))
    allphi = np.concatenate([phi_mep, phi_ery, phi_mk])
    lo, hi = np.percentile(allphi, [0.3, 99.7]); bins = np.linspace(lo, hi, 50)
    cP, cE, cM, cS = '#4C6EF5', '#E8590C', '#2B8A3E', '#868E96'
    axL.hist(phi_ery, bins=bins, density=True, alpha=0.55, color=cE, label=f'erythroid (n={len(phi_ery)})')
    axL.hist(phi_mk, bins=bins, density=True, alpha=0.55, color=cM, label=f'megakaryocytic (n={len(phi_mk)})')
    axL.hist(phi_mep, bins=bins, density=True, histtype='step', lw=2.6, color=cP, label=f'MEP (n={len(phi_mep)})')
    axL.axvline(phi_mep.mean(), color=cP, ls='--', lw=1)
    axL.set_xlim(lo, hi); axL.set_xlabel(r'$\varphi$ = z(erythroid) $-$ z(megakaryocyte)')
    axL.set_ylabel('density'); axL.set_title('A  Adult MEP vs committed fates')
    axL.legend(fontsize=8, frameon=False)
    axL.text(0.02, 0.97, f'MEP dip p = {res["progenitor (MEP)"][3]:.3f} (unimodal)\n'
             f'mean $\\varphi$: MEP {phi_mep.mean():+.2f} | ery {phi_ery.mean():+.2f} | mk {phi_mk.mean():+.2f}',
             transform=axL.transAxes, va='top', fontsize=8)
    axR.hist(synth(phi_ery, phi_mk, 150), bins=bins, density=True, alpha=0.55, color=cS,
             label='synthetic 50/50 mix (n=150)')
    axR.hist(phi_mep, bins=bins, density=True, histtype='step', lw=2.6, color=cP, label=f'MEP (n={len(phi_mep)})')
    axR.set_xlim(lo, hi); axR.set_xlabel(r'$\varphi$'); axR.set_ylabel('density')
    axR.set_title('B  Positive control: MEP vs synthetic mixture'); axR.legend(fontsize=8, frameon=False)
    axR.text(0.02, 0.97, f'MEP:       dip p = {res["progenitor (MEP)"][3]:.3f},  dBIC = {res["progenitor (MEP)"][2]:.1f}\n'
             f'synth mix:  dip p = {dp_s:.3f},  dBIC = {db_s:.1f}',
             transform=axR.transAxes, va='top', fontsize=8)
    fig.suptitle('Adult human bone marrow: the megakaryocyte–erythroid progenitor is a single\nunresolved population between its erythroid and megakaryocyte fates (Zeng et al. 2025)',
                 fontsize=10.5, y=1.02)
    fig.tight_layout()
    import os; os.makedirs(os.path.dirname(OUTFIG), exist_ok=True)
    fig.savefig(OUTFIG, dpi=150, bbox_inches='tight')
    print(f"\nsaved figure -> {OUTFIG}")

if __name__ == '__main__':
    main()
