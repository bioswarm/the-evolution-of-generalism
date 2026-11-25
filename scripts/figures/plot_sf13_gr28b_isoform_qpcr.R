#!/usr/bin/env Rscript
# Pure-source-code Supplementary Fig. 13 harmonized with the manuscript's exact palette (Fig. 1-5, SF14, SF16):
#   a: B. dorsalis (generalist; RefSeq variants X1 [Gr28b-X1], X3, X4, X5)
#   b: B. minax (specialist; RefSeq variants X1, X2 [Gr28b-X1 ortholog], X3, X5)
#   c: B. tsuneonis (specialist; RefSeq variants X1, X2 [Gr28b-X1 ortholog], X3, X4)
# Data loaded directly from raw_data/supplement.rawdata.xlsx (347 exact qRT-PCR replicates).
# Styling: soft academic pastel tissue palette matching SF16 & original SF13, species-colored titles matching Fig. 1/2/5,
#          clean white backgrounds (subtle ivory tint only on orthologous Gr28b-X1 facet), no gridlines, 100% TrueType Arial.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

REPO_ROOT <- Sys.getenv("REPO_ROOT", ".")
OUT_DIR   <- Sys.getenv("OUT_DIR", file.path(REPO_ROOT, "outputs", "figures"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dt <- fread(file.path(REPO_ROOT, "data", "SF13_qPCR_raw.tsv"))

tissue_ord <- c("Antenna", "Foreleg", "Mouthparts", "Ovipositor", "Brain", "Compound eye")
# Harmonized soft pastel palette matching original SF13, SF16, and Fig. 1/2/5 aesthetics
tissue_col <- c(
  "Antenna"      = "#F2B880",  # soft warm apricot (SF16 orange family)
  "Foreleg"      = "#9AD0A0",  # soft sage green (SF16 cucumber family)
  "Mouthparts"   = "#9BC0E6",  # soft sky blue (SF16 pumpkin family)
  "Ovipositor"   = "#DF8682",  # soft muted coral-crimson (B. dorsalis / ovipositor highlight)
  "Brain"        = "#C4B7DC",  # soft lavender (B. tau / SF14 purple family)
  "Compound eye" = "#E8A6C8"   # soft rose-pink (SF16 dragon fruit / orig SF13 eye)
)

dt[, tissue_fac := factor(tissue, levels = tissue_ord)]

# Deterministic horizontal offsets for replicate dots
dt <- dt[order(species, variant, tissue_fac, value)]
dt[, x_off := {
  n <- .N
  if (n == 1L) 0
  else if (n <= 3L) seq(-0.14, 0.14, length.out = n)
  else seq(-0.20, 0.20, length.out = n)
}, by = .(species, variant, tissue_fac)]

make_species_row <- function(sp_key, sp_label, sp_role, var_order, var_labels, orth_var, panel_letter) {
  d <- copy(dt[species == sp_key])
  d[, var_fac := factor(variant, levels = var_order, labels = var_labels)]
  d[, is_orth := (variant == orth_var)]
  d[, x_pos := as.numeric(tissue_fac) + x_off]

  sm <- d[, .(
    n = .N,
    mean_v = mean(value),
    sem_v  = sd(value) / sqrt(.N),
    lo     = pmax(0, mean(value) - sd(value) / sqrt(.N)),
    hi     = mean(value) + sd(value) / sqrt(.N)
  ), by = .(var_fac, is_orth, tissue_fac)]
  sm[, x_c := as.numeric(tissue_fac)]

  ymax <- 1.48
  ybrks <- c(0, 0.5, 1.0, 1.5)

  orth_df <- unique(sm[is_orth == TRUE, .(var_fac)])
  # Exact species colors from Fig. 1, Fig. 2, Fig. 5, SF19
  sp_col <- switch(sp_key,
                   "B. dorsalis"  = "#C62828",
                   "B. minax"     = "#1B7A3A",
                   "B. tsuneonis" = "#2B83BA")

  title_expr <- bquote(bold(.(panel_letter)) ~ "  " ~ bolditalic(.(sp_label)) ~ plain(.(paste0("(", sp_role, ")"))))

  ggplot() +
    geom_rect(data = orth_df, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
              fill = "#FAF6F0", inherit.aes = FALSE) +
    geom_col(data = sm, aes(x = x_c, y = mean_v, fill = tissue_fac),
             width = 0.62, alpha = 0.85, colour = "#333333", linewidth = 0.36, show.legend = FALSE) +
    geom_errorbar(data = sm, aes(x = x_c, ymin = lo, ymax = hi),
                  width = 0.22, linewidth = 0.48, colour = "#222222") +
    geom_point(data = d, aes(x = x_pos, y = value, fill = tissue_fac),
               shape = 21, size = 1.65, stroke = 0.36, colour = "#222222", show.legend = FALSE) +
    facet_wrap(~ var_fac, nrow = 1, labeller = label_parsed) +
    scale_fill_manual(values = tissue_col) +
    scale_x_continuous(breaks = seq_along(tissue_ord), labels = tissue_ord,
                       limits = c(0.4, length(tissue_ord) + 0.6), expand = c(0, 0)) +
    scale_y_continuous(breaks = ybrks, expand = c(0.015, 0)) +
    coord_cartesian(ylim = c(0, ymax), clip = "off") +
    labs(title = title_expr, x = NULL, y = "Relative expression") +
    theme_classic(base_size = 9.2, base_family = "Arial") +
    theme(
      plot.title = element_text(size = 10.2, hjust = 0, colour = sp_col, margin = margin(b = 5)),
      panel.background = element_rect(fill = "white", colour = NA),
      strip.background = element_rect(fill = "#F5F6F8", colour = "#CCCCCC", linewidth = 0.45),
      strip.text = element_text(size = 8.9, colour = "#111111", margin = margin(t = 3.5, b = 3.5)),
      panel.spacing = unit(9, "pt"),
      axis.line = element_line(colour = "#222222", linewidth = 0.50),
      axis.ticks = element_line(colour = "#222222", linewidth = 0.46),
      axis.ticks.length = unit(3, "pt"),
      axis.text.x = element_text(size = 8.2, colour = "#111111", angle = 36, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 8.5, colour = "#111111"),
      axis.title.y = element_text(size = 9.3, face = "bold", colour = "#111111", margin = margin(r = 6)),
      legend.position = "none",
      plot.margin = margin(6, 8, 6, 8)
    )
}

p_a <- make_species_row(
  "B. dorsalis", "B. dorsalis", "generalist",
  c("X1", "X3", "X4", "X5"),
  c("bold('X1 (')*bolditalic('Gr28b-X1')*bold(')')", "bold('X3')", "bold('X4')", "bold('X5')"),
  "X1", "a"
)
p_b <- make_species_row(
  "B. minax", "B. minax", "specialist",
  c("X1", "X2", "X3", "X5"),
  c("bold('X1')", "bold('X2 (')*bolditalic('Gr28b-X1')*bold(' ortholog)')", "bold('X3')", "bold('X5')"),
  "X2", "b"
)
p_c <- make_species_row(
  "B. tsuneonis", "B. tsuneonis", "specialist",
  c("X1", "X2", "X3", "X4"),
  c("bold('X1')", "bold('X2 (')*bolditalic('Gr28b-X1')*bold(' ortholog)')", "bold('X3')", "bold('X4')"),
  "X2", "c"
)

leg_df <- data.table(
  tissue = factor(tissue_ord, levels = tissue_ord),
  x = c(1.0, 2.05, 3.05, 4.30, 5.45, 6.35),
  y = 1
)
p_leg <- ggplot(leg_df, aes(x = x, y = y)) +
  annotate("text", x = 0.55, y = 1, label = "Tissue:", hjust = 1,
           size = 3.2, fontface = "bold", colour = "#111111", family = "Arial") +
  geom_point(aes(fill = tissue), shape = 22, size = 3.8, stroke = 0.45, colour = "#333333") +
  geom_text(aes(x = x + 0.14, label = tissue), hjust = 0, size = 3.1, colour = "#111111", family = "Arial") +
  annotate("rect", xmin = 7.75, xmax = 8.05, ymin = 0.80, ymax = 1.20,
           fill = "#FAF6F0", colour = "#BBBBBB", linewidth = 0.4) +
  annotate("text", x = 8.15, y = 1, label = "plain('Orthologous ') * italic('Gr28b-X1') * plain(' isoform')",
           parse = TRUE, hjust = 0, size = 3.05, colour = "#333333", family = "Arial") +
  scale_fill_manual(values = tissue_col, guide = "none") +
  coord_cartesian(xlim = c(0.0, 10.6), ylim = c(0.55, 1.45), expand = FALSE, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(2, 6, 4, 6))

p_all <- (p_a / p_b / p_c / p_leg) + plot_layout(heights = c(1, 1, 1, 0.16))

out_pdf <- file.path(OUT_DIR, "Supplementary_Figure_13_Gr28b_isoform_qPCR.pdf")
cairo_pdf(out_pdf, width = 10.2, height = 8.2, family = "Arial")
print(p_all)
dev.off()
cat("Wrote harmonized Supplementary Fig. 13 ->", out_pdf, "\n")
