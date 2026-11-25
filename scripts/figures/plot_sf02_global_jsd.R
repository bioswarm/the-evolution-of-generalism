#!/usr/bin/env Rscript
# Pure-source Supplementary Fig. 2 (Global 95-group JSD heatmap) in 100% TrueType Arial,
# using full standard species names (B. dorsalis, B. correcta, B. minax, B. cucurbitae, B. tau)
# with ZERO non-standard 4-letter abbreviations (Bdor/Bmin/Bcor/Bcuc/Btau).

REPO_ROOT   <- Sys.getenv("REPO_ROOT", ".")
REVIEW_ROOT <- Sys.getenv("REVIEW_ROOT", "/path/to/project")
OUT_DIR     <- Sys.getenv("OUT_DIR", file.path(REPO_ROOT, "outputs", "figures"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
REPRO_ROOT  <- Sys.getenv("REPRO_ROOT", file.path(REVIEW_ROOT, "reproduce"))
Sys.setenv(REPRO_ROOT = normalizePath(REPRO_ROOT))
source(file.path(REPRO_ROOT, "R", "00_common.R"))


suppressPackageStartupMessages({
  library(pheatmap)
  library(philentropy)
  library(grid)
})


ss <- sample_sheet(); setkey(ss, sample)
tp <- fig1c_tpm(ss)

P1 <- log2(tp + 2)
SP_FULL <- c(Bc = "B. correcta", Bd = "B. dorsalis", Bm = "B. minax",
             Zc = "B. cucurbitae", Zt = "B. tau")
grp <- with(ss, paste0(
  SP_FULL[sp], " (",
  ifelse(group == "LA", "La", c(FN = "nF", FM = "mF", MM = "mM")[group]),
  ", ", sub("larval_", "", tissue), ")"
))
G1 <- sapply(split(seq_along(grp), grp), function(k) rowMeans(P1[, k, drop = FALSE]))
G1 <- sweep(G1, 2, colSums(G1), "/")
D1 <- suppressMessages(philentropy::JSD(t(G1), unit = "log2"))
dimnames(D1) <- list(colnames(G1), colnames(G1))

ph1 <- pheatmap::pheatmap(
  D1,
  clustering_distance_rows = as.dist(D1),
  clustering_distance_cols = as.dist(D1),
  clustering_method = "complete",
  color = colorRampPalette(c("#4575b4", "#e0f3f8", "#fee090", "#d73027"))(100),
  fontsize = 6.2, fontfamily = "Arial", border_color = NA, silent = TRUE
)

out_pdf <- file.path(OUT_DIR, "Supplementary_Figure_2_global_JSD.pdf")
out_png <- file.path(OUT_DIR, "Supplementary_Figure_2_global_JSD.png")

cairo_pdf(out_pdf, width = 10.6, height = 10.6, family = "Arial")
grid.newpage(); grid.draw(ph1$gtable)
dev.off()

png(out_png, width = 10.6, height = 10.6, units = "in", res = 300, type = "cairo", family = "Arial")
grid.newpage(); grid.draw(ph1$gtable)
dev.off()

cat("Wrote Supplementary Fig. 2 (full species names, Arial) ->", out_pdf, "\n")
