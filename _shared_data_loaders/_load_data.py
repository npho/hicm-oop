"""Shared data loading helper for all figure scripts."""
import numpy as np
import pandas as pd
import os

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

def load_embryo_data():
    """Load counts, sample info, genes, and embedding.

    Returns two sets:
      - (Sample_Info, Counts, cESFW_Genes, Embedding): 2156 cells that have
        expression data — use for any analysis involving gene expression.
      - (Sample_Info_full, Embedding_full): all 2400 cells from the pre-computed
        embedding — use for UMAP scatter plots so the full topology is visible.
    """
    Sample_Info_all = pd.read_csv(os.path.join(DATA, "Human_Sample_Info.csv"),
                                  header=0, index_col=0)
    Counts_all = pd.read_csv(os.path.join(DATA, "Human_Embryo_Counts.csv"),
                             header=0, index_col=0)
    cESFW_Genes = np.load(os.path.join(DATA, "Saved_cESFW_Genes.npy"), allow_pickle=True)
    Embedding_all = np.load(os.path.join(DATA, "UMAP_Human_Embryo_Embedding.npy"))

    # The embedding rows correspond to Sample_Info row order (2400 cells).
    # Counts only has cells that were present in the merged matrix (2156 cells).
    # Subset everything to cells present in BOTH.
    common_cells = Sample_Info_all.index[
        np.isin(Sample_Info_all.index, Counts_all.index)
    ]
    si_inds = np.where(np.isin(Sample_Info_all.index, common_cells))[0]

    Sample_Info = Sample_Info_all.loc[common_cells]
    Counts      = Counts_all.loc[common_cells]
    Embedding   = Embedding_all[si_inds]

    Counts = Counts.fillna(0)
    cESFW_Genes = cESFW_Genes[np.isin(cESFW_Genes, Counts.columns)]

    print(f"Loaded: {Counts.shape[0]} cells × {Counts.shape[1]} genes | "
          f"{cESFW_Genes.shape[0]} cESFW genes | embedding {Embedding.shape} | "
          f"full embedding {Embedding_all.shape}")
    return (Sample_Info, Counts, cESFW_Genes, Embedding,
            Sample_Info_all, Embedding_all)
