#!/usr/bin/env Rscript
# Pure-source Supplementary Fig. 3 (Per-tissue PCA) in 100% TrueType Arial,
# using unified manuscript stage labels (mF, nF, mM, La — matching Fig. 1, Fig. 2, SF2, SF4, SF12).

REPO_ROOT   <- Sys.getenv("REPO_ROOT", ".")
REVIEW_ROOT <- Sys.getenv("REVIEW_ROOT", "/path/to/project")
OUT_DIR     <- Sys.getenv("OUT_DIR", file.path(REPO_ROOT, "outputs", "figures"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
REPRO_ROOT  <- Sys.getenv("REPRO_ROOT", file.path(REVIEW_ROOT, "reproduce"))
Sys.setenv(REPRO_ROOT = normalizePath(REPRO_ROOT))
source(file.path(REPRO_ROOT, "R", "00_common.R"))


suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
})


ss <- sample_sheet(); setkey(ss, sample)
tp <- fig1c_tpm(ss)

sf3_list <- lapply(c("ovipositor", "foreleg", "mouthpart", "antenna", "eye", "brain"), function(t) {
  i <- ss$tissue %in% c(t, paste0("larval_", t))
  x <- tp[, i]
  pp <- pca(x[rowSums(x >= 5) >= 1, ], scale = TRUE)
  p <- scatter(pp, ss[i], title = TIS_LAB[t], base = 8.5) +
    aes(shape = group) +
    scale_shape_manual(
      values = c(FM = 17, FN = 2, MM = 15, LA = 4),
      limits = c("FM", "FN", "MM", "LA"),
      labels = c(FM = "mF", FN = "nF", MM = "mM", LA = "La"),
      drop = FALSE, name = NULL
    ) +
    guides(
      colour = guide_legend(order = 1, override.aes = list(size = 2.4),
                            theme = theme(legend.text = element_text(size = 8.5, face = "italic", family = "Arial"))),
      shape  = guide_legend(order = 2, override.aes = list(size = 2.4),
                            theme = theme(legend.text = element_text(size = 8.5, face = "plain", family = "Arial")))
    )
  if (t != "mouthpart") p <- p + theme(legend.position = "none")
  p
})

p_sf3 <- wrap_plots(sf3_list, ncol = 3, guides = "collect")

out_pdf <- file.path(OUT_DIR, "Supplementary_Figure_3_per_tissue_PCA.pdf")
out_png <- file.path(OUT_DIR, "Supplementary_Figure_3_per_tissue_PCA.png")

cairo_pdf(out_pdf, width = 7.8, height = 5.2, family = "Arial")
print(p_sf3)
dev.off()

png(out_png, width = 7.8, height = 5.2, units = "in", res = 300, type = "cairo", family = "Arial")
print(p_sf3)
dev.off()

cat("Wrote Supplementary Fig. 3 (unified mF/nF/mM/La labels, Arial) ->", out_pdf, "\n")
