# 04_feature_selection.py
import pandas as pd
import numpy as np

# ── Load normalized matrices ───────────────────────────────────────────────────
rna = pd.read_csv("data/processed/rna_vst.tsv",    sep="\t", index_col=0)
mut = pd.read_csv("data/processed/mut_binary.tsv", sep="\t", index_col=0)
meth = pd.read_csv("data/processed/meth_mval.tsv",  sep="\t", index_col=0)

# Fix RNA column names — dots back to dashes
rna.columns = rna.columns.str.replace(".", "-", regex=False)

print(f"RNA shape:         {rna.shape}")
print(f"Mutation shape:    {mut.shape}")
print(f"Methylation shape: {meth.shape}")

# ── Align samples ─────────────────────────────────────────────────────────────
common_samples = list(
    set(rna.columns) & set(mut.columns) & set(meth.columns)
)
print(f"\nCommon samples: {len(common_samples)}")

rna = rna[common_samples]
mut = mut[common_samples]
meth = meth[common_samples]

# ── Filter RNA by MAD (top 5000 genes) ────────────────────────────────────────
mad = rna.apply(lambda x: np.median(np.abs(x - x.median())), axis=1)
top_rna_genes = mad.nlargest(5000).index
rna_filtered = rna.loc[top_rna_genes]
print(f"RNA after MAD filter:          {rna_filtered.shape}")

# ── Filter mutations by recurrence (>=5% of samples) ─────────────────────────
recurrence = mut.mean(axis=1)
mut_filtered = mut[recurrence >= 0.05]
print(f"Mutations after recurrence filter: {mut_filtered.shape}")

# ── Filter methylation by variance (top 10000 CpGs) ──────────────────────────
meth_var = meth.var(axis=1)
top_cpgs = meth_var.nlargest(10000).index
meth_filtered = meth.loc[top_cpgs]
print(f"Methylation after variance filter: {meth_filtered.shape}")

# ── Save filtered matrices ────────────────────────────────────────────────────
rna_filtered.to_csv("data/processed/rna_filtered.tsv",  sep="\t")
mut_filtered.to_csv("data/processed/mut_filtered.tsv",  sep="\t")
meth_filtered.to_csv("data/processed/meth_filtered.tsv", sep="\t")

print("\nFeature selection complete. Files saved to data/processed/")
