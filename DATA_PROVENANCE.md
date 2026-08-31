# Data provenance

Generated 2026-08-28.

**Read this first.** The scripts in this package were written for one workstation and
contain **absolute paths of the form `/Users/adelosan/cESFW_analysis/data/...`**. They
will not run elsewhere without repointing those paths. Most scripts obtain the embryo
reference through `_shared_data_loaders/_load_data_yanagida.py`, whose `DATA` constant
resolves to `<repo>/data` — set that one path and the majority of scripts work.

Only 11 accessions appear inside the scripts themselves, and only two scripts contain
anything resembling a download step. This file is therefore the authoritative record of
where the analysed files came from; it was compiled from the manuscript Methods, the
loader source, and the paths the scripts actually read.

## What each script reads

The complete per-script list of input files is in `data_inputs_by_figure.txt`.
The dominant input by far is the cESFW embryo atlas (used by ~43 of 65 scripts via
`load_embryo_data()`).

## Datasets

| Dataset | Accession | Files / role | Figures | Where to get it |
|---|---|---|---|---|
| Human embryo atlas (cESFW) | `Mendeley Data 34td4ds2r9` | Human_Embryo_Counts.csv, Human_Sample_Info.csv, Saved_cESFW_Genes.npy, UMAP_Human_Embryo_Embedding.npy — the reference used by ~43 of the 65 scripts | nearly all | https://data.mendeley.com/datasets/34td4ds2r9  (Radley et al. 2023; download and unpack into data/) |
| Yanagida 2021 blastocyst | `GSE171820` | GSE171820_Counts_raw.txt.gz + GSE171820_Sample_Assignments.csv.gz; appended to the base atlas by _load_data_yanagida.py (adds 228 blastocyst cells -> 2,384) | nearly all | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171820/suppl/ |
| Petropoulos 2016 | `E-MTAB-3929` | primary embryo dataset (inside the cESFW atlas) | 1,2,3,5,6,S19 | https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-3929 |
| Marmoset (Boroviak 2018) | `E-MTAB-7078` | counts + per-cell lineage labels (ICM_CS3/Epi_CS3/Hyp_CS3); labels from the Boroviak-Lab SpatialModelling GitHub repo | 6 | https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-7078 |
| Mouse (Nowotschin 2019) | `GSE123046` | 10x E3.5/E4.5 inner cells | 6 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123046/suppl/ |
| Mouse diapause | `GSE241462` | diapause single-cell | S5 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241462/suppl/ |
| Mouse paused ESC | `GSE81285` | Bulut-Karslioglu 2016 paused ESC | S5 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81285/suppl/ |
| Iyer 2024 blastoid | `GSE267302` | Ctrl vs RapaLink blastoid (mTOR inhibition) | 4 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267302/suppl/ |
| Human ICM ATAC (Wu 2018) | `GSE101571` | MACS peaks (hg19), 2PN/3PN ICM + naive/primed hESC; also refGene_hg19.txt.gz and primed_hESC_rpkm.wig.gz | 3,S16,S17 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE101nnn/GSE101571/suppl/ |
| Mouse ESC TF ChIP-seq (Chen 2008) | `GSE11431` | Esrrb GSM288355, Tcfcp2l1 GSM288350, Oct4/Sox2/Nanog/Klf4 + Myc-module control | 3,S18 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE11nnn/GSE11431/suppl/ |
| Cynomolgus (Nakamura 2016) | `GSE74767` | SC3seq_Cy_ProcessedData.txt.gz | S6 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74767/suppl/ |
| Rhesus | `GSE106481` | NHP conservation (exploratory) | S6 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106481/suppl/ |
| Messmer 2019 t2iLGo | `E-MTAB-6819` | naive PSC counts + sample info | 4 | https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6819 |
| Guo 2021 PXGL / Polo 2020 / Bi 2022 + Moya 2022 5iLAF / Mazid 2022 e4CL / hEPSC | `see Methods` | naive and alternative PSC datasets for reference mapping | 4 | per-dataset GEO/ArrayExpress accessions listed in the manuscript Methods |
| Dattani 2024 | `in-house processed` | naive->hypoblast timecourse; adata_n2b27_fax_fig5.h5ad plus derived CSVs in plots/ | 4 | Cell Stem Cell 31:1058-1071; processed object from the authors |
| Braun 2023 first-trimester brain | `see Methods` | braun_meta.tsv, braun_markers.tsv (thalamus subset) | 7 | CELLxGENE / published portal |
| Burclaff 2022 intestine | `UCSC Cell Browser "intestine-epithelium"` | gut_matrix.mtx.gz, gut_features.tsv.gz, gut_barcodes.tsv.gz, gut_meta.tsv | 7 | https://cells.ucsc.edu/?ds=intestine-epithelium |
| Arutyunyan 2023 placenta | `CELLxGENE Discover` | donor P13 trophoblast | 7 | https://cellxgene.cziscience.com/  (search Arutyunyan 2023 placenta) |
| Sikkema 2023 Human Lung Cell Atlas | `CELLxGENE` | healthy lung AT2/AT0/AT1 | S15 | https://cellxgene.cziscience.com/ |
| Habermann 2020 fibrotic lung | `GSE135893` | IPF lung AT2/transitional/AT1 | 7,S15 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/ |
| Zeng 2025 bone marrow | `see Methods` | balanced adult marrow reference (MEP) | S12 | published portal; see Methods |
| Gao 2018 ATAC peaks | `gao2018_atac` | icm_2pn_peaks.bed.gz | S16,S17 | see original publication |
| MSigDB Hallmark 2020 | `Enrichr GMT` | Enrichr.MSigDB_Hallmark_2020.gmt (cached at ~/.cache/gseapy/) | several | auto-downloaded by gseapy on first use |
| CollecTRI regulons | `decoupler / omnipath` | collectri_net.csv (fetched once, cached) | S2, TF activity | decoupler.get_collectri(organism="human") then .to_csv() |

## Caveats an auditor should know

1. **Some intermediate files are written to `/tmp`** (`/tmp/footprint_network.json`,
   `/tmp/chen2008_conservation.json`, `/tmp/icm_chromatin_priming.json`,
   `/tmp/ground_state_human.json`, `/tmp/progressive_assembly.json`,
   `/tmp/human_mod_magnitude.json`). These are produced by an upstream script in the
   same folder and consumed by the figure script — run the analysis script before the
   figure script. They are NOT inputs you need to download.
2. **Derived CSVs in `plots/`** (e.g. `ESRRB_embryo_ICM_positives.csv`,
   `ESRRB_ICMsig_dualcompetent.csv`, `dattani_icm_similarity.csv`,
   `ICM_signature_specificity_byStage.csv`, `ICM_NANOG_GATA6_balance.csv`) are outputs
   of scripts in this package, not primary data. They are the receipts behind the
   numbers quoted in the text.
3. **The bigwig inputs** (`corr_ICM/ICM_uncorrected.bw`, `corr_Naive/Naive_uncorrected.bw`)
   are produced by the TOBIAS pipeline (`ATACorrect`), not downloaded. See
   `Figure_S17/tobias_fork.py`.
4. Several PSC datasets are referenced by short local filenames
   (`Guo2021_sc_Counts.csv`, `Polo2020_Naive_QC_counts.csv`, `Moya2022_QC_counts.csv`,
   `Bi2022_5iLAF_nESC_QC_counts.csv`, `Mazid2022_sc_Counts.csv`, `HEPSC_Counts.csv`)
   that are **QC-filtered derivatives**, not raw downloads. Their accessions are in the
   manuscript Methods; the QC step is not included in this package.
