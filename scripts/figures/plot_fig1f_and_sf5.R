#!/usr/bin/env Rscript
# Redraw Fig. 1f (upregulated smDEGs) and Supplementary Fig. 5 (downregulated smDEGs)
# from pure source data (published_patterns + SuppFig4_column_order.tsv),
# displaying total smDEG counts directly on each tissue panel.
# 100% embedded TrueType Arial via cairo_pdf.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

REPO_ROOT   <- Sys.getenv("REPO_ROOT", ".")
REVIEW_ROOT <- Sys.getenv("REVIEW_ROOT", "/path/to/project")
OUT_DIR     <- Sys.getenv("OUT_DIR", file.path(REPO_ROOT, "outputs", "figures"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
REPRO_ROOT  <- Sys.getenv("REPRO_ROOT", file.path(REVIEW_ROOT, "reproduce"))
Sys.setenv(REPRO_ROOT = normalizePath(REPRO_ROOT))
source(file.path(REPRO_ROOT, "R", "00_common.R"))


cols_ord <- fread(file.path(REPO_ROOT, "data", "SuppFig4_column_order.tsv"),
                  skip = 1, colClasses = "character")$pattern
stopifnot(length(cols_ord) == 31L)

pub <- published_patterns()
tissue_order <- c("brain", "eye", "antenna", "mouthpart", "foreleg", "ovipositor")
tissue_labels <- c(
  brain      = "Brain",
  eye        = "Eye",
  antenna    = "Antenna",
  mouthpart  = "Mouthpart",
  foreleg    = "Foreleg",
  ovipositor = "Ovipositor"
)

sp_keys <- c("tau", "cucurbitae", "minax", "dorsalis", "correcta")
sp_ital <- c(
  tau        = "B. tau",
  cucurbitae = "B. cucurbitae",
  minax      = "B. minax",
  dorsalis   = "B. dorsalis",
  correcta   = "B. correcta"
)
sp_cols <- c(
  tau        = "#7b52ab",
  cucurbitae = "#a83c8b",
  minax      = "#1b7a3a",
  dorsalis   = "#c62828",
  correcta   = "#e67e22"
)

build_upset <- function(direction = c("up", "down"),
                        show_species_labels = TRUE,
                        highlight_generalist = (direction == "up"),
                        base_size = 7.5,
                        dot_size = 1.55) {
  direction <- match.arg(direction)
  ch <- if (direction == "up") "2" else "0"

  # Compute exact counts from published_patterns()
  inter <- lapply(pub[tissue_order], function(p) {
    b <- vapply(strsplit(p, ""), function(v) paste(as.integer(v == ch), collapse = ""), "")
    table(b[b != "00000"])
  })

  df <- rbindlist(lapply(tissue_order, function(t) {
    cnts <- inter[[t]]
    n_vec <- as.integer(ifelse(is.na(cnts[cols_ord]), 0L, cnts[cols_ord]))
    data.table(tissue = t, pattern = cols_ord, idx = seq_along(cols_ord), n = n_vec)
  }))
  df[, tissue := factor(tissue, levels = tissue_order)]
  df[, bar_fill := ifelse(highlight_generalist & tissue == "ovipositor" & pattern == "11011",
                          "highlight", "normal")]

  # Total smDEG count (union across 5 species = sum of all 31 disjoint intersection bars)
  tot_df <- df[, .(total = sum(n)), by = tissue]
  tot_df[, label_expr := sprintf('bold("%s") ~ "(total = %s)"',
                                 tissue_labels[as.character(tissue)],
                                 format(total, big.mark = ",", trim = TRUE))]

  ymax <- 325
  y_lab_txt <- if (direction == "up") {
    "#Upregulated genes after sexual maturity"
  } else {
    "#Downregulated genes after sexual maturity"
  }

  # Bar panels
  p_bars <- ggplot(df, aes(x = idx, y = n, fill = bar_fill)) +
    annotate("rect", xmin = 0.45, xmax = 6.55, ymin = -Inf, ymax = Inf,
             fill = "#f7f5e6", colour = NA) +
    annotate("rect", xmin = 26.45, xmax = 31.55, ymin = -Inf, ymax = Inf,
             fill = "#e3e3e3", colour = NA) +
    geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.35, colour = "#222222") +
    geom_col(width = 0.78, colour = NA) +
    geom_text(data = tot_df,
              aes(x = 16, y = 272, label = label_expr),
              parse = TRUE, inherit.aes = FALSE,
              family = "Arial", size = base_size / 2.85, colour = "#111111") +
    scale_fill_manual(values = c(normal = "#348abd", highlight = "#b83227"), guide = "none") +
    scale_x_continuous(limits = c(0.3, 31.7), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, ymax), breaks = c(50, 150, 250), expand = c(0, 0)) +
    facet_grid(tissue ~ .) +
    labs(x = NULL, y = y_lab_txt) +
    theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "#888888", fill = NA, linewidth = 0.45),
      panel.spacing.y = unit(1.0, "pt"),
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = base_size - 0.5, colour = "#111111"),
      axis.title.y = element_text(size = base_size + 0.5, colour = "#111111",
                                  margin = margin(r = 3)),
      axis.ticks.y = element_line(linewidth = 0.35, colour = "#222222"),
      axis.ticks.length.y = unit(2.2, "pt"),
      plot.background = element_rect(fill = NA, colour = NA),
      plot.margin = margin(t = 2, r = 3, b = 0, l = 2)
    )

  # Annotate "88" cleanly above the red generalist-shared ovipositor bar
  if (highlight_generalist) {
    hi_df <- df[tissue == "ovipositor" & pattern == "11011"]
    p_bars <- p_bars +
      geom_text(data = hi_df, aes(x = idx, y = n + 12, label = n),
                inherit.aes = FALSE, family = "Arial", fontface = "bold",
                vjust = 0, size = (base_size - 0.8) / 2.85, colour = "#b83227")
  }

  # Dot matrix below
  dm <- CJ(idx = seq_along(cols_ord), k = 1:5)
  dm[, pattern := cols_ord[idx]]
  dm[, on := substr(pattern, k, k) == "1"]
  dm[, sp := sp_keys[k]]
  dm[, y := 6 - k]  # 5 = tau (top), 1 = correcta (bottom)

  dm_on <- dm[on == TRUE]
  seg <- dm_on[, .(ymin = min(y), ymax = max(y)), by = idx]
  seg <- seg[ymax > ymin]

  p_dots <- ggplot() +
    annotate("rect", xmin = 0.45, xmax = 6.55, ymin = 0.35, ymax = 5.65,
             fill = "#f7f5e6", colour = NA) +
    annotate("rect", xmin = 26.45, xmax = 31.55, ymin = 0.35, ymax = 5.65,
             fill = "#e3e3e3", colour = NA) +
    geom_segment(data = seg, aes(x = idx, xend = idx, y = ymin, yend = ymax),
                 linewidth = 0.42, colour = "#1a1a1a") +
    geom_point(data = dm_on, aes(x = idx, y = y, colour = sp), size = dot_size) +
    scale_colour_manual(values = sp_cols, guide = "none") +
    scale_x_continuous(limits = c(0.3, 31.7), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1.6, 5.75), breaks = 5:1,
                       labels = sp_ital[sp_keys], expand = c(0, 0)) +
    annotate("text", x = 3.5, y = -0.55, label = "Common DEGs",
             family = "Arial", size = base_size / 2.85, colour = "#111111", vjust = 1) +
    annotate("text", x = 29.0, y = -0.45, label = "Lineage-specific\nDEGs",
             family = "Arial", size = base_size / 2.85, colour = "#111111",
             vjust = 1, lineheight = 0.92) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = base_size, base_family = "Arial") +
    theme(
      axis.text.y = if (show_species_labels) {
        element_text(size = base_size - 0.3, face = "italic", colour = "#111111",
                     hjust = 1, margin = margin(r = 3))
      } else {
        element_blank()
      },
      plot.background = element_rect(fill = NA, colour = NA),
      plot.margin = margin(t = 1, r = 3, b = 3, l = 2)
    )

  (p_bars / p_dots + plot_layout(heights = c(6.1, 1.35))) &
    theme(plot.background = element_rect(fill = NA, colour = NA))
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 1. Standalone Fig. 1f (with species labels on left, publication vector Arial)
p_f_standalone <- build_upset("up", show_species_labels = TRUE,
                              highlight_generalist = TRUE, base_size = 8.0, dot_size = 1.75)
cairo_pdf(file.path(OUT_DIR, "Figure1f_upregulated_smDEGs.pdf"), width = 5.2, height = 5.4, family = "Arial")
print(p_f_standalone)
dev.off()

# 2. Exact A4 overlay PDF to splice seamlessly onto Figure1.pdf (595.28 x 841.89 pt)
# Exact verified coordinates (y_top from pdftotext -bbox):
#   Panel letter 'f':    x=246.2..249.5, y_top=254.9..271.3  (PRESERVED!)
#   Panel e legend:      xMax=265.07 (cucurbitae), y_top=460.4..560 (PRESERVED!)
#   Old panel f y-title: x=254.8..264.3, y_top=328..456
#   Old panel f ticks:   x=269.1..279.1, y_top=263..522
p_f_embed <- build_upset("up", show_species_labels = FALSE,
                         highlight_generalist = TRUE, base_size = 7.4, dot_size = 1.55)

W_pt <- 595.28
H_pt <- 841.89
cairo_pdf(file.path(OUT_DIR, "Figure1f_overlay_A4.pdf"),
          width = W_pt / 72, height = H_pt / 72,
          family = "Arial", bg = "transparent")
grid.newpage()
# Region 1 (above panel e legend, y_top 254..458.5): x = 252..582
grid.rect(x = unit(252 / W_pt, "npc"), y = unit((H_pt - 458.5) / H_pt, "npc"),
          width = unit((582 - 252) / W_pt, "npc"), height = unit((458.5 - 254) / H_pt, "npc"),
          just = c("left", "bottom"), gp = gpar(fill = "white", col = NA))
# Region 2 (alongside panel e legend, y_top 458..595): x = 267..582 (strictly right of cucurbitae 265.07!)
grid.rect(x = unit(267 / W_pt, "npc"), y = unit((H_pt - 595) / H_pt, "npc"),
          width = unit((582 - 267) / W_pt, "npc"), height = unit((595 - 458) / H_pt, "npc"),
          just = c("left", "bottom"), gp = gpar(fill = "white", col = NA))

# Place redrawn panel f into viewport x = 255..578 pt, y_top = 255..593 pt
vp_f <- viewport(
  x = unit(255 / W_pt, "npc"),
  y = unit((H_pt - 593) / H_pt, "npc"),
  width = unit((578 - 255) / W_pt, "npc"),
  height = unit((593 - 255) / H_pt, "npc"),
  just = c("left", "bottom")
)
pushViewport(vp_f)
print(p_f_embed, newpage = FALSE)
popViewport()
dev.off()

# 3. Matching Supplementary Fig. 5 (downregulated smDEGs with total counts on each tissue panel)
p_sf5 <- build_upset("down", show_species_labels = TRUE,
                     highlight_generalist = FALSE, base_size = 8.0, dot_size = 1.75)
cairo_pdf(file.path(OUT_DIR, "Supplementary_Figure_5_downregulated_smDEGs.pdf"), width = 5.4, height = 5.5, family = "Arial")
print(p_sf5)
dev.off()

cat("Successfully wrote Fig1f_upregulated_smDEGs.pdf, Fig1f_overlay_A4.pdf, SuppFig5_downregulated_smDEGs.pdf\n")
