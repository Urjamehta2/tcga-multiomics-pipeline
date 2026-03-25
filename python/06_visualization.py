# 06_visualization.py
# Visualizations for multi-omic driver gene results

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats

# ── Load data ─────────────────────────────────────────────────────────────────
rna       = pd.read_csv("data/processed/rna_filtered.tsv",  sep="\t", index_col=0)
mut       = pd.read_csv("data/processed/mut_filtered.tsv",  sep="\t", index_col=0)
meth      = pd.read_csv("data/processed/meth_filtered.tsv", sep="\t", index_col=0)
scores    = pd.read_csv("outputs/cca_sample_scores.tsv",    sep="\t", index_col=0)
drivers   = pd.read_csv("outputs/driver_candidates.tsv",    sep="\t")["gene"].tolist()
rna_load  = pd.read_csv("outputs/top_rna_loadings.tsv",     sep="\t", index_col=0)

common_samples = list(
    set(rna.columns) & set(mut.columns) & set(meth.columns)
)

# ── Plot 1: CCA biplot (samples in CC1 vs CC2) ────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(scores["CC1"], scores["CC2"], alpha=0.6, s=20, color="#1D9E75")
ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
ax.axvline(0, color="gray", linewidth=0.5, linestyle="--")
ax.set_xlabel("Canonical Component 1")
ax.set_ylabel("Canonical Component 2")
ax.set_title("CCA biplot — TCGA-BRCA samples")
plt.tight_layout()
plt.savefig("outputs/figures/01_cca_biplot.png", dpi=150)
plt.close()
print("✓ CCA biplot saved")

# ── Plot 2: Top 30 RNA loadings (lollipop chart) ──────────────────────────────
top30 = rna_load.nlargest(30, "loading")
fig, ax = plt.subplots(figsize=(8, 8))
ax.hlines(y=top30.index, xmin=0, xmax=top30["loading"],
          color="#378ADD", linewidth=1.5)
ax.scatter(top30["loading"], top30.index, color="#378ADD", s=40, zorder=3)
ax.set_xlabel("Absolute CCA loading (CC1)")
ax.set_title("Top 30 RNA genes — CC1 loading")
ax.invert_yaxis()
plt.tight_layout()
plt.savefig("outputs/figures/02_rna_loadings_lollipop.png", dpi=150)
plt.close()
print("✓ Lollipop chart saved")

# ── Plot 3: Volcano plot (mutated vs wildtype expression) ─────────────────────
if len(drivers) > 0:
    results = []
    for gene in drivers:
        if gene not in rna.columns and gene in rna.index:
            mutated_samples = mut.columns[mut.loc[gene] == 1].tolist() \
                if gene in mut.index else []
            wt_samples = [s for s in common_samples if s not in mutated_samples]
            if len(mutated_samples) >= 5 and len(wt_samples) >= 5:
                expr_mut = rna.loc[gene, mutated_samples].values
                expr_wt  = rna.loc[gene, wt_samples].values
                t, p     = stats.ttest_ind(expr_mut, expr_wt)
                log2fc   = np.mean(expr_mut) - np.mean(expr_wt)
                results.append({"gene": gene, "log2FC": log2fc, "pval": p})

    if results:
        df_vol = pd.DataFrame(results)
        df_vol["-log10p"] = -np.log10(df_vol["pval"].clip(lower=1e-300))
        df_vol["significant"] = (
            (df_vol["pval"] < 0.05) & (df_vol["log2FC"].abs() > 1)
        )
        fig, ax = plt.subplots(figsize=(8, 6))
        colors = df_vol["significant"].map({True: "#D85A30", False: "#888780"})
        ax.scatter(df_vol["log2FC"], df_vol["-log10p"],
                   c=colors, alpha=0.7, s=30)
        for _, row in df_vol[df_vol["significant"]].iterrows():
            ax.annotate(row["gene"], (row["log2FC"], row["-log10p"]),
                        fontsize=7, alpha=0.8)
        ax.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.8)
        ax.axvline(-1, color="gray", linestyle="--", linewidth=0.8)
        ax.axvline(1,  color="gray", linestyle="--", linewidth=0.8)
        ax.set_xlabel("log2 Fold Change (mutated vs wildtype)")
        ax.set_ylabel("-log10(p-value)")
        ax.set_title("Volcano plot — driver gene candidates")
        plt.tight_layout()
        plt.savefig("outputs/figures/03_volcano.png", dpi=150)
        plt.close()
        print("✓ Volcano plot saved")

# ── Plot 4: Multi-omic heatmap (top 30 driver candidates) ─────────────────────
if len(drivers) > 0:
    top_drivers = drivers[:30]
    avail_rna  = [g for g in top_drivers if g in rna.index]
    avail_mut  = [g for g in top_drivers if g in mut.index]

    if avail_rna:
        fig = plt.figure(figsize=(14, 8))
        gs  = gridspec.GridSpec(2, 1, height_ratios=[4, 1], hspace=0.05)

        # RNA expression heatmap
        ax_rna = fig.add_subplot(gs[0])
        rna_sub = rna.loc[avail_rna, common_samples[:100]]  # first 100 samples
        sns.heatmap(rna_sub, ax=ax_rna, cmap="RdBu_r", center=0,
                    xticklabels=False, yticklabels=True,
                    cbar_kws={"label": "VST expression"})
        ax_rna.set_title("Multi-omic heatmap — driver gene candidates")
        ax_rna.set_ylabel("Gene")

        # Mutation bar
        if avail_mut:
            ax_mut = fig.add_subplot(gs[1])
            mut_sub = mut.loc[avail_mut, common_samples[:100]].mean(axis=0)
            ax_mut.bar(range(len(mut_sub)), mut_sub.values,
                       color="#D85A30", width=1.0)
            ax_mut.set_xlim(0, len(mut_sub))
            ax_mut.set_ylabel("Mut.\nfreq.")
            ax_mut.set_xlabel("Samples")

        plt.savefig("outputs/figures/04_multiomics_heatmap.png",
                    dpi=150, bbox_inches="tight")
        plt.close()
        print("✓ Multi-omic heatmap saved")

print("\nAll visualizations saved to outputs/figures/")
