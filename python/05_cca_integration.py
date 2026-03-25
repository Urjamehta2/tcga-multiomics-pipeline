# 05_cca_integration.py
from scipy import stats
import pandas as pd
import numpy as np
from sklearn.cross_decomposition import CCA
from sklearn.preprocessing import StandardScaler

# ── Load filtered matrices ────────────────────────────────────────────────────
rna = pd.read_csv("data/processed/rna_filtered.tsv",  sep="\t", index_col=0)
mut = pd.read_csv("data/processed/mut_filtered.tsv",  sep="\t", index_col=0)
meth = pd.read_csv("data/processed/meth_filtered.tsv", sep="\t", index_col=0)

# ── Align samples ─────────────────────────────────────────────────────────────
common_samples = list(
    set(rna.columns) & set(mut.columns) & set(meth.columns)
)
print(f"Samples going into CCA: {len(common_samples)}")

rna = rna[common_samples]
mut = mut[common_samples]
meth = meth[common_samples]

# ── Scale inputs ──────────────────────────────────────────────────────────────
X = StandardScaler().fit_transform(rna.T.values)
Y = StandardScaler().fit_transform(meth.T.values)

# ── Run CCA ───────────────────────────────────────────────────────────────────
print("Running CCA (this may take a few minutes)...")
cca = CCA(n_components=10, max_iter=1000)
cca.fit(X, Y)
X_c, Y_c = cca.transform(X, Y)
print("CCA complete.")

# ── Extract loadings ──────────────────────────────────────────────────────────
rna_loadings = pd.Series(cca.x_weights_[:, 0], index=rna.index)
meth_loadings = pd.Series(cca.y_weights_[:, 0], index=meth.index)

# ── Layer 1: top RNA genes by CCA loading ────────────────────────────────────
rna_top = set(rna_loadings.abs().nlargest(500).index)

# ── Layer 2: recurrently mutated genes ───────────────────────────────────────
mut_recurrent = set(mut.index)

# ── Layer 3: differential expression in mutated vs wildtype ──────────────────

de_results = []
for gene in mut_recurrent:
    if gene in rna.index:
        mutated = mut.loc[gene]
        mut_samples = mutated[mutated == 1].index.tolist()
        wt_samples = mutated[mutated == 0].index.tolist()
        if len(mut_samples) >= 5 and len(wt_samples) >= 5:
            expr_mut = rna.loc[gene, mut_samples].values
            expr_wt = rna.loc[gene, wt_samples].values
            t, p = stats.ttest_ind(expr_mut, expr_wt)
            log2fc = np.mean(expr_mut) - np.mean(expr_wt)
            de_results.append({
                "gene": gene,
                "log2FC": log2fc,
                "pval": p,
                "in_rna_top": gene in rna_top
            })

de_df = pd.DataFrame(de_results)
de_df["significant"] = (de_df["pval"] < 0.05) & (de_df["log2FC"].abs() > 0.5)

print("\n── Differentially expressed mutated genes ──")
print(de_df.sort_values("pval")[
      ["gene", "log2FC", "pval", "significant"]].to_string(index=False))

# ── Driver candidates: mutated + differentially expressed ────────────────────
driver_candidates = de_df[de_df["significant"]]["gene"].tolist()
print(f"\nDriver gene candidates: {len(driver_candidates)}")
print(sorted(driver_candidates))

# ── Save results ──────────────────────────────────────────────────────────────
scores = pd.DataFrame(X_c, index=common_samples,
                      columns=[f"CC{i+1}" for i in range(10)])
scores.to_csv("outputs/cca_sample_scores.tsv", sep="\t")

rna_loadings.abs().nlargest(500).to_csv(
    "outputs/top_rna_loadings.tsv", sep="\t", header=["loading"]
)

de_df.to_csv("outputs/de_mutated_genes.tsv", sep="\t", index=False)

pd.Series(sorted(driver_candidates)).to_csv(
    "outputs/driver_candidates.tsv", sep="\t", index=False, header=["gene"]
)

print("\nResults saved to outputs/")
