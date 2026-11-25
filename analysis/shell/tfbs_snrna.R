#!/usr/bin/env Rscript
# % of Gr28b+ ovipositor nuclei expressing each candidate TF (snRNA object is the
# B. dorsalis ovipositor SCT object, features are Drosophila-style gene names).
# Usage: Rscript shell/tfbs_snrna.R <out_dir>
suppressMessages({ library(Seurat) })
out <- commandArgs(trailingOnly = TRUE)[1]
o <- readRDS("/path/to/data/scRNA/raw_matrix/analysis/AjOvi-v1-SCT.rds")
DefaultAssay(o) <- "RNA"
cm <- GetAssayData(o, assay = "RNA", layer = "counts")
gp <- cm["Gr28b", ] > 0
rn <- rownames(cm); rl <- tolower(rn)
tfs <- readLines(file.path(out, "tf_names.txt"))
res <- t(sapply(tfs, function(tf) {
  i <- match(tolower(tf), rl)
  c(tf, if (is.na(i)) NA else mean(cm[rn[i], gp] > 0), sum(gp), is.na(i))
}))
write.table(data.frame(tf = res[, 1], pct_Gr28bplus = round(as.numeric(res[, 2]), 4),
                       n_Gr28b_nuclei = as.integer(res[, 3]), feature_missing = res[, 4] == "TRUE"),
            file.path(out, "snrna_detection.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
