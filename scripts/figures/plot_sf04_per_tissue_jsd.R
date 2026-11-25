#!/usr/bin/env Rscript
# Pure-source reproduction of Supplementary Fig. 3 (per-tissue pairwise JSD
# heatmaps, adult tissues only — larval brain & mouthpart excluded) in 100%
# embedded TrueType Arial, matching the original white-to-crimson sequential
# palette, row-shared JSD scales, top UPGMA dendrograms, species-colored labels,
# and primary cluster highlight boxes.

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
  library(philentropy)
  library(grid)
})

F_OUT <- OUT_DIR

# Load sample sheet & published Fig. 1c ortholog TPM matrix
# Note: gene_expression_wide_output.txt already carries the corrected Bc-MM-L/M
# labels (verified r = 0.966-0.969 vs corrected ortho_tpm); do NOT double-swap.
ss <- sample_sheet(); setkey(ss, sample)
tp <- as.matrix(fread(CFG$fig1c_tpm), rownames = 1)
tp <- tp[, !grepl("^Bt", colnames(tp))]
colnames(tp) <- gsub("-", ".", sub("-LA-?", "-", colnames(tp)))
stopifnot(setequal(colnames(tp), ss$sample))
tp <- tp[, ss$sample]

# Exclude larval samples (larval brain and larval mouthpart)
keep <- ss$group != "LA"
ss_a <- ss[keep]
tp_a <- tp[, keep]

sp_short <- c(Bc = "Bcor", Bd = "Bdor", Bm = "Bmin", Zc = "Bcuc", Zt = "Btau")
st_short <- c(FN = "nF", FM = "mF", MM = "mM")
grp_a <- paste0(sp_short[ss_a$sp], ".", st_short[ss_a$group], ".", ss_a$tissue)

P <- sweep(tp_a, 2, colSums(tp_a), "/")
G <- sapply(split(seq_along(grp_a), grp_a), function(k) rowMeans(P[, k, drop = FALSE]))

# Species label colors (consistent across all 6 panels; B. minax bold green)
LBL_COL <- c(
  Bmin = "#0F662F",
  Bcor = "#DE7700",
  Bdor = "#C62828",
  Bcuc = "#984EA3",
  Btau = "#6A51A3"
)

# White-to-crimson sequential palette matching original submitted SuppFig 3
JSD_PAL <- c("#FFFFFF", "#FADDDD", "#E99B9B", "#D14C4C", "#A8101C")

# Convert hclust object to dendrogram line segments, orienting smaller cluster left
build_dendro_segments <- function(hc, leaf_order) {
  n <- length(hc$labels)
  pos <- setNames(seq_len(n), leaf_order)
  m <- hc$merge
  h <- hc$height
  node_x <- numeric(n - 1)
  node_y <- h
  node_min <- logical(n - 1) # TRUE if all descendants are B. minax
  segs <- list()

  for (i in seq_len(n - 1)) {
    kids <- m[i, ]
    xs <- numeric(2); ys <- numeric(2); is_m <- logical(2)
    for (j in 1:2) {
      if (kids[j] < 0) {
        lbl <- hc$labels[-kids[j]]
        xs[j] <- pos[[lbl]]
        ys[j] <- 0
        is_m[j] <- startsWith(lbl, "Bmin.")
      } else {
        xs[j] <- node_x[kids[j]]
        ys[j] <- node_y[kids[j]]
        is_m[j] <- node_min[kids[j]]
      }
    }
    node_x[i] <- mean(xs)
    node_min[i] <- all(is_m)
    # Two vertical legs + one horizontal bridge
    segs[[length(segs) + 1]] <- data.frame(x = xs[1], xend = xs[1], y = ys[1], yend = h[i], is_min = is_m[1])
    segs[[length(segs) + 1]] <- data.frame(x = xs[2], xend = xs[2], y = ys[2], yend = h[i], is_min = is_m[2])
    segs[[length(segs) + 1]] <- data.frame(x = xs[1], xend = xs[2], y = h[i],  yend = h[i], is_min = node_min[i])
  }
  do.call(rbind, segs)
}

# Reorder hclust dendrogram so that the smaller top-level cluster (or mF in
# ovipositor, and B. minax when inside a cluster) sits on the left/top
orient_order <- function(hc, tissue) {
  dend <- as.dendrogram(hc)
  # Score each leaf so smaller/focal cluster goes first (smaller weight = left)
  ct <- cutree(hc, 2)
  sz <- table(ct)
  if (tissue == "ovipositor") {
    # Put mF first, Bmin.mF leftmost
    w <- ifelse(grepl("\\.mF$", hc$labels), 0, 100) +
         ifelse(startsWith(hc$labels, "Bmin."), 0, 10)
  } else {
    small_cl <- as.integer(names(which.min(sz)))
    w <- ifelse(ct[hc$labels] == small_cl, 0, 100) +
         ifelse(startsWith(hc$labels, "Bmin."), 0, 10)
  }
  names(w) <- hc$labels
  dend2 <- reorder(dend, w, agglo.FUN = mean)
  hc2 <- as.hclust(dend2)
  hc2
}

make_tissue_panel <- function(tissue, zmax, breaks, show_legend = FALSE) {
  M <- G[, grepl(paste0("\\.", tissue, "$"), colnames(G)), drop = FALSE]
  colnames(M) <- sub(paste0("\\.", tissue, "$"), "", colnames(M))
  D <- suppressMessages(philentropy::JSD(t(M), unit = "log2"))
  dimnames(D) <- list(colnames(M), colnames(M))

  # Save updated TSV table
  fwrite(data.table(group = rownames(D), D),
         RES("structure", paste0("SuppFig3_JSD_", tissue, ".tsv")), sep = "\t")

  hc_raw <- hclust(as.dist(D), method = "average")
  hc <- orient_order(hc_raw, tissue)
  ord <- hc$labels[hc$order]
  n <- length(ord)

  # Primary 2-cluster blocks along the ordered diagonal
  ct <- cutree(hc, 2)[ord]
  rle_ct <- rle(ct)
  ends <- cumsum(rle_ct$lengths)
  starts <- ends - rle_ct$lengths + 1
  boxes <- data.frame(
    xmin = starts - 0.5,
    xmax = ends + 0.5,
    ymin = (n - ends + 1) - 0.5,
    ymax = (n - starts + 1) + 0.5
  )

  # Heatmap long table (row 1 at TOP -> y = n..1)
  df <- expand.grid(col = ord, row = ord, stringsAsFactors = FALSE)
  df$val <- as.vector(t(D[ord, ord]))
  df$x <- match(df$col, ord)
  df$y <- n - match(df$row, ord) + 1

  # Label metadata for colored axis ticks (full standard species names, ZERO Bmin/Bcor abbreviations!)
  SP_FULL <- c(Bmin = "B. minax", Bcor = "B. correcta", Bdor = "B. dorsalis",
               Bcuc = "B. cucurbitae", Btau = "B. tau")
  fmt_lbl <- function(vec) {
    sp_k <- sub("\\..*$", "", vec)
    st_k <- sub("^[^.]+\\.", "", vec)
    paste0(SP_FULL[sp_k], " (", st_k, ")")
  }
  sp_ord <- sub("\\..*$", "", ord)
  lbl_cols_x <- unname(LBL_COL[sp_ord])
  lbl_face_x <- ifelse(sp_ord == "Bmin", "bold.italic", "italic")
  disp_x <- fmt_lbl(ord)
  # For y-axis, y=1 is bottom (rev(ord)) and y=n is top (ord[1])
  ord_y <- rev(ord)
  sp_y <- sub("\\..*$", "", ord_y)
  lbl_cols_y <- unname(LBL_COL[sp_y])
  lbl_face_y <- ifelse(sp_y == "Bmin", "bold.italic", "italic")
  disp_y <- fmt_lbl(ord_y)

  # Top dendrogram plot
  segs <- build_dendro_segments(hc, ord)
  hmax <- max(segs$yend)
  p_dend <- ggplot(segs) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend, colour = is_min),
                 linewidth = 0.48, lineend = "square") +
    scale_colour_manual(values = c(`TRUE` = "#0F662F", `FALSE` = "#222222"), guide = "none") +
    scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, hmax * 1.06), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(2, 0, 1, 0))

  # Main heatmap with TOP rotated labels, RIGHT horizontal labels, LEFT vertical tissue title
  fs_lbl <- if (n <= 10) 8.2 else 7.4
  p_hm <- ggplot(df, aes(x = x, y = y, fill = val)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = NA, colour = "#114890", linewidth = 1.05) +
    scale_fill_gradientn(
      colours = JSD_PAL,
      limits = c(0, zmax),
      oob = scales::squish,
      breaks = breaks,
      labels = sprintf("%.2f", breaks),
      name = "Jensen\u2013Shannon\nDivergence",
      guide = if (show_legend) guide_colourbar(
        barheight = unit(2.6, "cm"), barwidth = unit(0.36, "cm"),
        title.position = "left", title.hjust = 0.5,
        frame.colour = "grey75", frame.linewidth = 0.3,
        ticks.colour = "grey30", ticks.linewidth = 0.3
      ) else "none"
    ) +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = seq_len(n), labels = disp_x, position = "top"
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      breaks = seq_len(n), labels = disp_y, position = "right"
    ) +
    coord_fixed(xlim = c(0.5, n + 0.5), ylim = c(0.5, n + 0.5), clip = "off") +
    labs(x = NULL, y = NULL) +
    annotate("text", x = 0.5 - 0.065 * n, y = (n + 1) / 2, label = TIS_LAB[tissue],
             angle = 90, fontface = "bold", size = 3.55, family = "Arial",
             colour = "#111111", hjust = 0.5, vjust = 0.5) +
    theme_minimal(base_size = 8.5, base_family = "Arial") +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x.top = element_text(
        angle = 90, hjust = 0, vjust = 0.5,
        size = fs_lbl, colour = lbl_cols_x, face = lbl_face_x,
        margin = margin(b = 2)
      ),
      axis.text.y.right = element_text(
        hjust = 0, vjust = 0.5,
        size = fs_lbl, colour = lbl_cols_y, face = lbl_face_y,
        margin = margin(l = 2)
      ),
      legend.position = if (show_legend) "left" else "none",
      legend.title = element_text(angle = 90, size = 8.5, face = "plain", colour = "#111111",
                                  margin = margin(r = 5)),
      legend.text = element_text(size = 7.8, colour = "#111111"),
      legend.margin = margin(0, 2, 0, 0),
      plot.margin = margin(0, 14, 4, 24)
    )

  # Stack dendrogram directly above heatmap with matching x-axis width
  if (show_legend) {
    # Extract legend separately so dendrogram stays aligned strictly over the heatmap
    leg <- cowplot::get_legend(p_hm)
    p_hm_noleg <- p_hm + theme(legend.position = "none")
    core <- (p_dend / p_hm_noleg) + plot_layout(heights = c(0.24, 1))
    list(core = core, leg = wrap_elements(leg))
  } else {
    (p_dend / p_hm) + plot_layout(heights = c(0.24, 1))
  }
}

# Build Row 1 (Ovipositor, Foreleg, Mouthpart; shared scale 0 - 0.35)
r1_ovi <- make_tissue_panel("ovipositor", zmax = 0.35, breaks = c(0.05, 0.15, 0.25, 0.35), show_legend = TRUE)
r1_leg <- make_tissue_panel("foreleg",    zmax = 0.35, breaks = c(0.05, 0.15, 0.25, 0.35), show_legend = FALSE)
r1_mou <- make_tissue_panel("mouthpart",  zmax = 0.35, breaks = c(0.05, 0.15, 0.25, 0.35), show_legend = FALSE)

# Build Row 2 (Antenna, Eye, Brain; shared scale 0 - 0.14)
r2_ant <- make_tissue_panel("antenna",    zmax = 0.14, breaks = c(0.02, 0.06, 0.10, 0.14), show_legend = TRUE)
r2_eye <- make_tissue_panel("eye",        zmax = 0.14, breaks = c(0.02, 0.06, 0.10, 0.14), show_legend = FALSE)
r2_bra <- make_tissue_panel("brain",      zmax = 0.14, breaks = c(0.02, 0.06, 0.10, 0.14), show_legend = FALSE)

row1 <- wrap_plots(r1_ovi$leg, r1_ovi$core, r1_leg, r1_mou, ncol = 4, widths = c(0.26, 1, 1, 1))
row2 <- wrap_plots(r2_ant$leg, r2_ant$core, r2_eye, r2_bra, ncol = 4, widths = c(0.26, 1, 1, 1))
fig3 <- (row1 / row2) + plot_layout(heights = c(1, 1)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

pdf_path <- file.path(F_OUT, "Supplementary_Figure_4_per_tissue_JSD.pdf")
png_path <- file.path(F_OUT, "Supplementary_Figure_4_per_tissue_JSD.png")

W <- 13.6; H <- 8.6
suppressWarnings({
  cairo_pdf(pdf_path, width = W, height = H, family = "Arial")
  print(fig3)
  dev.off()

  png(png_path, width = W, height = H, units = "in", res = 300, type = "cairo", family = "Arial")
  print(fig3)
  dev.off()
})

log_msg("Saved ", pdf_path, " and ", png_path)
