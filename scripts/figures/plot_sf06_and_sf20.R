#!/usr/bin/env Rscript
# Build Supplementary Fig. 6 (maturation PCA + divergence) and Supplementary Fig. 20 (Gr subfamily branch lengths) in 100% embedded TrueType Arial:
#   1. Combined Supplementary Fig. 6 (a: per-tissue maturation PCA + b: specialist divergence index)
#   2. Single-panel Gr subfamily total branch length (left panel only, highlighting Gr28b)

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

# ==============================================================================
# 1. Combined SF6: (a) Per-tissue maturation-effect PCA + (b) Divergence index
# ==============================================================================
ab <- c(correcta = "Bc", dorsalis = "Bd", minax = "Bm", cucurbitae = "Zc", tau = "Zt")
tt <- c(Ovi = "Ovipositor", Leg = "Foreleg", Mo = "Mouthpart", Ant = "Antenna", Eye = "Eye", Br = "Brain")
sf_pca_list <- lapply(c("Ovi", "Leg", "Mo", "Ant", "Eye", "Br"), function(t) {
  B <- as.matrix(read.table(file.path(CFG$published_patterns, t, "mashR", "mashR.beta"), header = TRUE, row.names = 1))
  q <- prcomp(t(B), center = TRUE, scale. = TRUE)
  v <- 100 * q$sdev^2 / sum(q$sdev^2)
  dd <- data.table(sp = ab[colnames(B)], PC1 = q$x[, 1], PC2 = q$x[, 2])
  p <- ggplot(dd, aes(PC1, PC2, colour = sp)) + geom_point(size = 2.8) +
    scale_colour_manual(values = SP_COL[c("Bc", "Bd", "Bm", "Zc", "Zt")],
                        labels = SP_LAB[c("Bc", "Bd", "Bm", "Zc", "Zt")], name = NULL) +
    labs(title = tt[t], x = sprintf("PC1 (%.1f%%)", v[1]), y = sprintf("PC2 (%.1f%%)", v[2])) +
    theme_pub(8.5) +
    theme(plot.title = element_text(size = 9, hjust = 0.5, face = "bold"),
          legend.text = element_text(size = 8.2, face = "italic"))
  p
})
pa <- wrap_plots(sf_pca_list, ncol = 3, guides = "collect") & theme(legend.position = "right")

dv14 <- fread(RES("divergence", "tissue_divergence.tsv"))
SPS <- c("Bc", "Bd", "Bm", "Zc", "Zt")
d14b <- melt(dv14, id.vars = c("measure", "tissue", "Bm_CI_low", "Bm_CI_high"), measure.vars = paste0("index_", SPS))
d14b[, ref_grp := ifelse(variable == "index_Bm", "Bm", "Gen")]
d14b[, tissue := factor(TIS_LAB[tissue], levels = c("Antenna", "Brain", "Eye", "Foreleg", "Mouthpart", "Ovipositor"))]
d14b[, mlab := ifelse(measure == "maturation_response", "Maturation response", "Mature female expression")]
pb <- ggplot(d14b, aes(value, tissue)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50", linewidth = 0.45) +
  geom_errorbarh(data = d14b[variable == "index_Bm"], aes(xmin = Bm_CI_low, xmax = Bm_CI_high),
                 height = 0, colour = "#eb6834", linewidth = 0.75) +
  geom_point(aes(colour = ref_grp, size = ref_grp)) +
  scale_size_manual(values = c(Bm = 3.1, Gen = 2.1), guide = "none") +
  scale_colour_manual(
    values = c(Gen = "#999999", Bm = "#eb6834"),
    breaks = c("Gen", "Bm"),
    labels = list("Generalist as reference", bquote(italic("B. minax") ~ "as reference")),
    name = NULL
  ) +
  facet_wrap(~ mlab, scales = "free_x") +
  labs(x = "Divergence index (reference / others)", y = NULL) +
  theme_pub(9.2) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9.0),
        strip.text = element_text(size = 9.5, face = "bold"),
        axis.title.x = element_text(size = 9.2),
        axis.text.y = element_text(size = 8.8),
        axis.text.x = element_text(size = 8.5))

p_comb <- (wrap_elements(full = pa) / pb) +
  plot_layout(heights = c(1.38, 0.92)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 11.5, face = "bold", family = "Arial"))

cairo_pdf(file.path(OUT_DIR, "Supplementary_Figure_6_maturation_PCA_and_divergence.pdf"), width = 8.0, height = 8.2, family = "Arial")
print(p_comb); dev.off()
png(file.path(OUT_DIR, "Supplementary_Figure_6_maturation_PCA_and_divergence.png"), width = 8.0, height = 8.2, units = "in", res = 300, type = "cairo", family = "Arial")
print(p_comb); dev.off()
log_msg("Wrote SuppFig6_maturation_pca_divergence.{pdf,png}")

# ==============================================================================
# 2. Single-panel (left-only) Gr subfamily total branch length
# ==============================================================================
gr_dt <- fread(RES("tables", "SuppData9_Gr_subfamily_branch_lengths.tsv"))
gr_dt <- gr_dt[order(clade_length)]
lbl <- ifelse(gr_dt$drosophila_sister == "-" | gr_dt$drosophila_sister == "" | is.na(gr_dt$drosophila_sister),
              paste0("no Drosophila ortholog (#", seq_len(nrow(gr_dt)), ")"),
              gr_dt$drosophila_sister)
lbl <- sub("^Gr22a,.*Gr22f$", "Gr22a\u2013f", lbl)
dup <- duplicated(lbl) | duplicated(lbl, fromLast = TRUE)
lbl[dup] <- paste0(lbl[dup], " (#", which(dup), ")")
gr_dt[, label := factor(lbl, levels = rev(lbl))]
gr_dt[, is_gr28b := as.logical(is_Gr28b)]

p_gr <- ggplot(gr_dt, aes(clade_length, label, fill = is_gr28b)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c(`FALSE` = "#a8a8a8", `TRUE` = "#e6007e"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(x = "Total branch length of Bactrocera clade", y = NULL) +
  theme_pub(9.0) +
  theme(axis.text.y = element_text(
          size = 8.2,
          colour = ifelse(rev(gr_dt$is_gr28b), "#e6007e", "#111111"),
          face = ifelse(rev(gr_dt$is_gr28b), "bold.italic", "plain")
        ),
        axis.title.x = element_text(size = 9.2, face = "bold"),
        plot.margin = margin(8, 12, 8, 8))

suppressWarnings({
  cairo_pdf(file.path(OUT_DIR, "Supplementary_Figure_20_Gr_subfamily_branch_lengths.pdf"), width = 5.4, height = 5.4, family = "Arial")
  print(p_gr); dev.off()
  png(file.path(OUT_DIR, "Supplementary_Figure_20_Gr_subfamily_branch_lengths.png"), width = 5.4, height = 5.4, units = "in", res = 300, type = "cairo", family = "Arial")
  print(p_gr); dev.off()
})
log_msg("Wrote SuppFig20_Gr_branch_length_left.{pdf,png}")
