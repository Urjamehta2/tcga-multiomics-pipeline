# 03_normalize.R
# Normalizes preprocessed matrices and exports TSV files for Python

library(DESeq2)
library(sva)

setwd("C:/Users/urjam/OneDrive/Desktop/PhD/OMIC")

# ── Load preprocessed data ────────────────────────────────────────────────────
rna_counts <- readRDS("data/processed/rna_counts.rds")
mut_matrix <- readRDS("data/processed/mut_matrix.rds")
meth_matrix <- readRDS("data/processed/meth_matrix.rds")

message("Preprocessed data loaded.")

# ── RNA-seq: VST normalization ────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
    countData = rna_counts,
    colData   = data.frame(row.names = colnames(rna_counts), condition = "BRCA"),
    design    = ~1
)
dds <- estimateSizeFactors(dds)
vst_mat <- assay(vst(dds, blind = TRUE))

message(paste("VST matrix dimensions:", nrow(vst_mat), "x", ncol(vst_mat)))

# ── RNA-seq: batch correction ─────────────────────────────────────────────────
# Extract plate ID from barcode (positions 22-25)
plate_id <- substr(colnames(vst_mat), 22, 25)
if (length(unique(plate_id)) > 1) {
    vst_mat <- ComBat(dat = vst_mat, batch = plate_id)
    message("Batch correction applied to RNA.")
}

# ── Methylation: beta to M-value conversion ───────────────────────────────────
# M = log2(beta / (1 - beta)), clip beta away from 0 and 1
meth_clipped <- pmin(pmax(meth_matrix, 0.001), 0.999)
meth_mval <- log2(meth_clipped / (1 - meth_clipped))

message(paste("M-value matrix dimensions:", nrow(meth_mval), "x", ncol(meth_mval)))

# ── Mutation: already binary, just align sample barcodes ─────────────────────
# Trim barcodes to 15 characters for consistency
colnames(mut_matrix) <- substr(colnames(mut_matrix), 1, 15)
colnames(vst_mat) <- substr(colnames(vst_mat), 1, 15)
colnames(meth_mval) <- substr(colnames(meth_mval), 1, 15)

# ── Export TSV files for Python ───────────────────────────────────────────────
write.table(vst_mat, "data/processed/rna_vst.tsv",
    sep = "\t", quote = FALSE
)
write.table(mut_matrix, "data/processed/mut_binary.tsv",
    sep = "\t", quote = FALSE
)
write.table(meth_mval, "data/processed/meth_mval.tsv",
    sep = "\t", quote = FALSE
)

message("Normalization complete. TSV files saved to data/processed/")
message("You can now run python/04_feature_selection.py")
