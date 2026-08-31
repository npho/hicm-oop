"""
Extended data loader that appends Yanagida 2021 blastocyst cells
(GSE171820, Origin='Blastocyst', n=228) to the base 2,156-cell atlas.

Drop-in replacement for _load_data.py: same return signature.
Combined atlas: up to 2,384 cells (2,156 base + Yanagida blastocyst cells
that are also present in Human_Sample_Info.csv).
"""
import numpy as np
import pandas as pd
import os

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

def load_embryo_data():
    # ── Base atlas ────────────────────────────────────────────────────────────
    Sample_Info_all = pd.read_csv(os.path.join(DATA, "Human_Sample_Info.csv"),
                                  header=0, index_col=0)
    Counts_all = pd.read_csv(os.path.join(DATA, "Human_Embryo_Counts.csv"),
                             header=0, index_col=0)
    cESFW_Genes = np.load(os.path.join(DATA, "Saved_cESFW_Genes.npy"),
                          allow_pickle=True)
    Embedding_all = np.load(os.path.join(DATA, "UMAP_Human_Embryo_Embedding.npy"))

    # ── Parse Yanagida 2021 blastocyst counts ─────────────────────────────────
    ya_path = os.path.join(DATA, "Blastoids_Yanagida2021",
                           "GSE171820_Counts_raw.txt.gz")
    sa_path = os.path.join(DATA, "Blastoids_Yanagida2021",
                           "GSE171820_Sample_Assignments.csv.gz")

    ya_raw = pd.read_csv(ya_path, sep=r'\s+', index_col=1,
                         quotechar='"', engine='python')
    ya_raw.index   = ya_raw.index.str.strip('"')
    ya_raw.columns = [c.strip('"') for c in ya_raw.columns]
    ya_raw = ya_raw.drop(columns=['gene_name', 'Length'], errors='ignore')
    # Convert column separators: file uses '-', Sample_Info uses '.'
    ya_raw.columns = [c.replace('-', '.') for c in ya_raw.columns]

    sa = pd.read_csv(sa_path, index_col=0)
    blasto_cells = sa[sa['Origin'] == 'Blastocyst'].index

    # Keep only blastocyst cells present in both count matrix AND Sample_Info
    ya_cells = [c for c in ya_raw.columns
                if c in blasto_cells and c in Sample_Info_all.index]
    ya_counts = ya_raw[ya_cells].T          # cells × genes

    # Align Yanagida genes to base atlas gene space
    shared_genes = [g for g in Counts_all.columns if g in ya_counts.columns]
    ya_aligned   = pd.DataFrame(0, index=ya_counts.index,
                                columns=Counts_all.columns, dtype=float)
    ya_aligned[shared_genes] = ya_counts[shared_genes].values

    # ── Combine ───────────────────────────────────────────────────────────────
    Counts_combined = pd.concat([Counts_all, ya_aligned])
    Counts_combined = Counts_combined.fillna(0)

    common_cells = Sample_Info_all.index[
        np.isin(Sample_Info_all.index, Counts_combined.index)
    ]
    si_inds = np.where(np.isin(Sample_Info_all.index, common_cells))[0]

    Sample_Info = Sample_Info_all.loc[common_cells]
    Counts      = Counts_combined.loc[common_cells]
    Embedding   = Embedding_all[si_inds]

    cESFW_Genes = cESFW_Genes[np.isin(cESFW_Genes, Counts.columns)]

    n_ya = len(ya_cells)
    print(f"Loaded: {Counts.shape[0]} cells × {Counts.shape[1]} genes | "
          f"{cESFW_Genes.shape[0]} cESFW genes | "
          f"embedding {Embedding.shape} | "
          f"full embedding {Embedding_all.shape} | "
          f"+{n_ya} Yanagida blastocyst cells")

    return (Sample_Info, Counts, cESFW_Genes, Embedding,
            Sample_Info_all, Embedding_all)
