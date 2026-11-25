#!/usr/bin/env Rscript
# Pure-source-code reproduction of the 5-panel Figure 2 (a-e):
#   a: 88 generalist-co-upregulated genes across 5 species x 6 tissues (13-gene blue clade box)
#   b: snRNA-seq UMAP of 17 cell populations in sexually mature B. dorsalis ovipositor (16,107 nuclei)
#   c: Vertical dot plot of the 11 detected focal genes across the 17 populations
#      (dot size = % expressed in decontX-cleaned RNA counts; colour = z-score of log1p mean SCT expression)
#   d: Expression of Gr28b (decontX-cleaned log-normalised counts) on UMAP with cluster 13, 9, 17 badges
#      and a magnified inset of cluster 17 (mechanosensory neurons: 11/34 = 32% Gr28b+)
#   e: Cross-species, multi-tissue Gr28b expression heatmap with mashr maturation significance
#
# Outputs:
#   reproduce/figures/final/Fig2_full.pdf   (pure vector, A4-proportioned)
#   reproduce/figures/final/Fig2_full.png   (300 dpi)
#   reproduce/figures/final/Fig2_be.png     (300 dpi, panels b-e only for response letter)

REPO_ROOT   <- Sys.getenv("REPO_ROOT", ".")
REVIEW_ROOT <- Sys.getenv("REVIEW_ROOT", "/path/to/project")
OUT_DIR     <- Sys.getenv("OUT_DIR", file.path(REPO_ROOT, "outputs", "figures"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
REPRO_ROOT  <- Sys.getenv("REPRO_ROOT", file.path(REVIEW_ROOT, "reproduce"))
Sys.setenv(REPRO_ROOT = normalizePath(REPRO_ROOT))
source(file.path(REPRO_ROOT, "R", "00_common.R"))


suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(grid)
  library(gridExtra)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(data.table)
})

# ------------------------------------------------------------------------------
# 0. Load data & verify exact manuscript counts
# ------------------------------------------------------------------------------
RDS_PATH <- Sys.getenv("SNRNA_RDS", file.path(REVIEW_ROOT, "improved_20261113/reproduce/results/snRNA/BdOV_decontX_reannot.rds"))
stopifnot(file.exists(RDS_PATH))
s <- readRDS(RDS_PATH)
ord <- SN_ORD
stopifnot(setequal(ord, unique(as.character(s$final_celltype))))
s$final_celltype <- factor(as.character(s$final_celltype), levels = ord)
Idents(s) <- "final_celltype"
DefaultAssay(s) <- "SCT"

ct_n <- table(s$final_celltype)[ord]
stopifnot(
  ncol(s) == 16107,
  ct_n["Tracheal epithelium"] == 551,
  ct_n["Epithelium (glob1+ ITP+)"] == 480,
  ct_n["Gr28b+ epithelium (Mdr49+ NaPi-T+)"] == 286,
  ct_n["Mechanosensory neurons (nSyb+ ChAT+ nompC+)"] == 34
)

# Log-normalised decontX-cleaned RNA assay for exact detection % and FeaturePlot
s <- NormalizeData(s, assay = "RNA", verbose = FALSE)

pal17 <- setNames(hue_pal()(length(ord)), seq_along(ord))
num   <- setNames(seq_along(ord), ord)

# ------------------------------------------------------------------------------
# 1. Panel a: 88-gene co-upregulated heatmap (exact 13-gene blue clade box)
# ------------------------------------------------------------------------------
pub   <- fread(file.path(CFG$published_patterns, "Ovi", "DEG.patterns"), header = FALSE, skip = 1, colClasses = "character")
genes88 <- pub[V2 %in% c("22122", "22022"), V1]
stopifnot(length(genes88) == 88)

ss_all <- sample_sheet()
ss_fm  <- ss_all[group == "FM"]
TP     <- fig1c_tpm(ss_all); rownames(TP) <- sub("_.*", "", rownames(TP))
TPM    <- ortho_tpm(ss_all)

TIS <- c("ovipositor", "foreleg", "mouthpart", "antenna", "eye", "brain")
SPO <- c(Bc = "B. correcta", Bd = "B. dorsalis", Bm = "B. minax", Zc = "B. cucurbitae", Zt = "B. tau")
SP_COLORS_A <- c(Bc = "#e8862a", Bd = "#c62828", Bm = "#1b7a3a", Zc = "#b84fa0", Zt = "#8c5fb8")

blocks_raw <- lapply(names(SPO), function(sp_code) {
  m <- sapply(TIS, function(t) rowMeans(TP[genes88, ss_fm[sp == sp_code & tissue == t, sample], drop = FALSE]))
  z <- t(scale(t(m))); z[is.na(z)] <- 0; t(z)
})
names(blocks_raw) <- names(SPO)
hc <- hclust(dist(t(blocks_raw[["Bc"]])), method = "complete")

blocks <- lapply(names(SPO), function(sp_code) {
  m <- sapply(TIS, function(t) rowMeans(log2(TPM[genes88, ss_fm[sp == sp_code & tissue == t, sample], drop = FALSE] + 1)))
  z <- t(scale(t(m))); z[is.na(z)] <- 0; t(z)
})
names(blocks) <- names(SPO)
Z88 <- do.call(rbind, blocks)

ovi_bc_raw <- blocks_raw[["Bc"]]["ovipositor", ] > 1.7
clades <- lapply(seq_len(nrow(hc$merge)), function(i) {
  f <- function(k) if (k < 0) -k else c(f(hc$merge[k, 1]), f(hc$merge[k, 2]))
  f(i)
})
score <- vapply(clades, function(m) {
  if (length(m) >= 10 && length(m) <= 20) sum(ovi_bc_raw[m]) - 0.5 * sum(!ovi_bc_raw[m]) else -Inf
}, 0)
box13 <- colnames(Z88)[clades[[which.max(score)]]]
stopifnot(length(box13) == 13)

hl <- lapply(names(SPO), function(sp_code) {
  rng <- if (sp_code == "Bm") c(-3, 3) else c(-2, 2)
  Heatmap(
    blocks[[sp_code]],
    name = SPO[[sp_code]],
    col = colorRamp2(c(rng[1], 0, rng[2]),
                     c("white", colorRampPalette(c("white", SP_COLORS_A[[sp_code]]))(5)[2], SP_COLORS_A[[sp_code]])),
    cluster_rows = FALSE,
    cluster_columns = if (sp_code == "Bc") hc else FALSE,
    column_order = if (sp_code == "Bc") NULL else hc$order,
    show_column_names = FALSE,
    row_labels = tools::toTitleCase(TIS),
    row_names_gp = gpar(fontsize = 7.3, fontfamily = "Arial"),
    rect_gp = gpar(col = "grey65", lwd = 0.25),
    row_title = SPO[[sp_code]],
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 8.2, fontface = "bold.italic", col = SP_COLORS_A[[sp_code]], fontfamily = "Arial"),
    column_dend_height = unit(8, "mm"),
    height = unit(6 * 2.55, "mm"),
    heatmap_legend_param = list(
      title = NULL, at = c(rng[1], 0, rng[2]),
      legend_height = unit(10, "mm"),
      grid_width = unit(2.6, "mm"),
      labels_gp = gpar(fontsize = 6.3, fontfamily = "Arial")
    )
  )
})
ht_a <- Reduce(`%v%`, hl)

draw_panel_a <- function() {
  hd <- draw(ht_a, gap = unit(1.3, "mm"), merge_legends = FALSE, newpage = FALSE,
             padding = unit(c(1, 2, 1, 2), "mm"))
  ord_cols <- unlist(column_order(hd))
  pos <- range(match(box13, colnames(Z88)[ord_cols]))
  stopifnot(diff(pos) + 1 == length(box13))
  nc <- ncol(Z88)
  for (sp_name in SPO) {
    decorate_heatmap_body(sp_name, {
      grid.rect(x = unit((pos[1] - 1) / nc, "npc"),
                width = unit((diff(pos) + 1) / nc, "npc"),
                just = "left",
                gp = gpar(col = "#1f6fb5", fill = NA, lwd = 1.25))
    })
  }
}

# ------------------------------------------------------------------------------
# 2. Panel b: UMAP of 17 cell populations + clean numbered circle legend
# ------------------------------------------------------------------------------
um <- as.data.frame(Embeddings(s, "umap"))
um$ct <- s$final_celltype
um$k  <- as.integer(num[as.character(um$ct)])
um$k_fac <- factor(um$k, levels = seq_along(ord))

# Hand-tuned non-overlapping label positions (data coordinates)
cen <- um %>%
  group_by(k, ct) %>%
  summarise(x = median(umap_1), y = median(umap_2), n = n(), .groups = "drop")

# Small offsets to prevent overlapping adjacent clusters
offsets <- tibble::tribble(
  ~k,  ~dx,   ~dy,
   1,  0.0,   0.0,
   2,  0.0,   0.0,
   3,  0.0,   0.0,
   4,  1.1,   0.4,
   5,  0.0,   1.4,
   6,  0.0,   0.0,
   7, -1.2,  -1.1,
   8,  0.4,  -0.3,
   9, -0.3,   1.2,
  10, -0.6,  -1.4,
  11, -1.2,  -0.2,
  12,  0.2,   1.3,
  13, -0.7,   1.4,
  14,  1.4,  -0.2,
  15,  1.2,   0.9,
  16,  1.2,  -1.1,
  17, -1.7,  -0.6
)
cen <- cen %>% left_join(offsets, by = "k") %>% mutate(lx = x + dx, ly = y + dy)

# Plot cluster 17 last so it is never buried
um_ord <- um %>% arrange(k == 17)

pb_umap <- ggplot(um_ord, aes(umap_1, umap_2, colour = k_fac)) +
  geom_point(size = 0.28, stroke = 0, alpha = 0.92) +
  scale_colour_manual(values = pal17, guide = "none") +
  geom_text(data = cen %>% filter(k != 17), aes(lx, ly, label = k),
            inherit.aes = FALSE, size = 2.15, colour = "grey10", fontface = "plain", family = "Arial") +
  # Make cluster 17 label bold and unmistakable
  geom_point(data = cen %>% filter(k == 17), aes(lx, ly), inherit.aes = FALSE,
             size = 3.6, colour = pal17["17"], alpha = 0.95) +
  geom_text(data = cen %>% filter(k == 17), aes(lx, ly, label = k),
            inherit.aes = FALSE, size = 2.0, colour = "white", fontface = "bold", family = "Arial") +
  coord_cartesian(xlim = c(-17, 15.5), ylim = c(-13.5, 16)) +
  labs(x = "umap_1", y = "umap_2") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35),
    axis.text = element_text(size = 6.5, colour = "black"),
    axis.title = element_text(size = 7.2, colour = "black"),
    plot.margin = margin(3, 2, 2, 2)
  )

# Display labels (rename cluster 13 consistently with other Epithelium clusters)
disp_ord <- ifelse(ord == "Gr28b+ epithelium (Mdr49+ NaPi-T+)",
                   "Epithelium (Gr28b+ Mdr49+ NaPi-T+)", ord)

# Custom numbered colored-circle legend (exact aesthetic of Figure2_V1.pdf)
leg_df <- tibble::tibble(
  k = seq_along(ord),
  y = rev(seq_along(ord)),
  label = paste0(disp_ord, " (", as.integer(ct_n[ord]), ")"),
  col = pal17[as.character(seq_along(ord))],
  is_hi = k %in% c(13, 17)
)

pb_leg <- ggplot(leg_df, aes(x = 0, y = y)) +
  geom_point(aes(colour = factor(k)), size = 3.0) +
  scale_colour_manual(values = pal17, guide = "none") +
  geom_text(aes(label = k), size = 1.85, colour = "grey10", fontface = "bold", family = "Arial") +
  geom_text(aes(x = 0.16, label = label, fontface = ifelse(is_hi, "bold", "plain")),
            hjust = 0, size = 1.98, colour = "black", family = "Arial") +
  scale_x_continuous(limits = c(-0.14, 4.65), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, 17.6), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(2, 0, 2, 0))

# ------------------------------------------------------------------------------
# 3. Panel c: Vertical DotPlot of 11 focal genes across 17 populations
#    Dot size  = % Expressed (decontX-cleaned RNA counts > 0, matches text 35% / 32% / 20%)
#    Dot color = z-score of log1p(AverageExpression SCT)
# ------------------------------------------------------------------------------
focal <- c("Gr28b", "NKCC", "Bd1G02800", "Esp", "Yp3", "Bd4G00648", "Obp56c", "Cyp318a1", "yellow-g", "Bd4G00733", "Pde8")
alias <- c(Bd4G00733 = "CG3191", Bd1G02800 = "CG6675")
feat  <- sapply(focal, function(g) {
  q <- if (g %in% names(alias)) alias[[g]] else g
  grep(paste0("^", q, "(\\.1)?$"), rownames(s[["SCT"]]), value = TRUE)[1]
})

# % expressed from RNA/decontX counts > 0
cnt_rna <- GetAssayData(s, assay = "RNA", layer = "counts")[feat, ]
pct_mat <- sapply(ord, function(ct) rowMeans(cnt_rna[, s$final_celltype == ct, drop = FALSE] > 0) * 100)
rownames(pct_mat) <- names(feat); colnames(pct_mat) <- as.character(seq_along(ord))

# z-score from SCT log1p(AverageExpression)
avg_sct <- as.matrix(AverageExpression(s, features = feat, assays = "SCT", group.by = "final_celltype")$SCT)
avg_sct <- avg_sct[feat, ord, drop = FALSE]
rownames(avg_sct) <- names(feat); colnames(avg_sct) <- as.character(seq_along(ord))
z_mat <- t(scale(t(log1p(avg_sct)))); z_mat[is.na(z_mat)] <- 0

# Sanity check exact Gr28b numbers
stopifnot(
  round(pct_mat["Gr28b", "13"], 1) == 34.6,
  round(pct_mat["Gr28b", "17"], 1) == 32.4,
  round(pct_mat["Gr28b", "9"],  1) == 19.6
)

dot_df <- expand.grid(gene = focal, k = seq_along(ord), stringsAsFactors = FALSE) %>%
  mutate(
    pct = mapply(function(g, kk) pct_mat[g, as.character(kk)], gene, k),
    z   = mapply(function(g, kk) z_mat[g, as.character(kk)], gene, k),
    z_clip = pmin(pmax(z, -1), 3.2),
    gene = factor(gene, levels = focal),
    k_fac = factor(k, levels = rev(seq_along(ord)))
  )

# Coloured gene x-axis labels (Gr28b bold red italic, others italic)
gene_labels <- lapply(focal, function(g) {
  if (g == "Gr28b") bquote(bolditalic(.(g))) else bquote(italic(.(g)))
})
gene_cols <- ifelse(focal == "Gr28b", "#b2182b", "black")

# Subtle row highlight boxes for cluster 13 and cluster 17 (optional ultra-light tint)
row_hi <- data.frame(
  ymin = c(17 - 13 + 0.5, 17 - 17 + 0.5),
  ymax = c(17 - 13 + 1.5, 17 - 17 + 1.5),
  xmin = 0.5, xmax = length(focal) + 0.5
)

pc_dot <- ggplot(dot_df, aes(x = gene, y = k_fac)) +
  geom_rect(data = row_hi, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "#fff3e0", alpha = 0.55) +
  geom_tile(fill = NA, colour = "grey92", linewidth = 0.2) +
  geom_point(aes(size = pct, fill = z_clip), shape = 21, colour = "grey30", stroke = 0.22) +
  scale_size_area(
    max_size = 4.0,
    limits = c(0, 100),
    breaks = c(10, 30, 60, 90),
    name = "% Expressed"
  ) +
  scale_fill_gradientn(
    colours = c("#f5f2ef", "#fdd0a2", "#fd8d3c", "#d94801", "#7f2704"),
    values  = rescale(c(-1, -0.2, 0.5, 1.8, 3.2)),
    limits  = c(-1, 3.2),
    breaks  = c(0, 1, 2, 3),
    name    = "Expression\n(z-score)"
  ) +
  scale_x_discrete(position = "top", labels = gene_labels) +
  scale_y_discrete(labels = function(x) x) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 6.5, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text.x.top = element_text(angle = 42, hjust = 0, vjust = 0, size = 6.2, colour = gene_cols),
    axis.text.y = element_text(size = 6.2, colour = ifelse(rev(seq_along(ord)) %in% c(13, 17), "black", "grey25"),
                               face = ifelse(rev(seq_along(ord)) %in% c(13, 17), "bold", "plain")),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 5.5, face = "plain"),
    legend.text  = element_text(size = 5.2),
    legend.key.height = unit(2.8, "mm"),
    legend.key.width  = unit(2.4, "mm"),
    legend.spacing.y  = unit(1.2, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 1, 2, 1)
  ) +
  guides(
    fill = guide_colorbar(order = 1, barwidth = unit(2.2, "mm"), barheight = unit(13, "mm")),
    size = guide_legend(order = 2, override.aes = list(fill = "#d94801", colour = "grey30"))
  )

# ------------------------------------------------------------------------------
# 4. Panel d: Expression of Gr28b (decontX log-normalised RNA) + badges 13, 9, 17
#    + magnified zoom-in inset of cluster 17 (mechanosensory neurons, 11/34 = 32%)
# ------------------------------------------------------------------------------
um$gr28b <- as.numeric(GetAssayData(s, assay = "RNA", layer = "data")["Gr28b", ])
um_gr <- um %>% arrange(gr28b)   # plot positive cells on top
gr_cap <- quantile(um$gr28b[um$gr28b > 0], 0.98)  # clean color dynamic range (~2.5)
um_gr$gr_plot <- pmin(um_gr$gr28b, gr_cap)

# Badge coordinates for clusters 13, 9, 17 (placed cleanly beside each cluster)
badges_d <- tibble::tibble(
  k    = c(13, 9, 17),
  bx   = c(-8.6, -10.9, -5.8),
  by   = c(8.6,   6.0,  -6.2),
  col  = pal17[c("13", "9", "17")],
  tcol = c("grey10", "grey10", "white")
)

# Zoom box around cluster 17 (x in [-8.45, -7.50], y in [-6.80, -5.85])
c17 <- um_gr %>% filter(k == 17)
zbox <- list(xmin = -8.45, xmax = -7.50, ymin = -6.80, ymax = -5.85)
# Inset placement in open bottom-left space of Panel d (wide enough for crisp title + 34 cells)
ins  <- list(xmin = -16.8, xmax = -8.2, ymin = -13.4, ymax = -7.7)

# Connector lines from zoom box to inset box
conn_df <- tibble::tibble(
  x    = c(zbox$xmin, zbox$xmax),
  y    = c(zbox$ymax, zbox$ymin),
  xend = c(ins$xmin,  ins$xmax),
  yend = c(ins$ymax,  ins$ymax)
)

pd_main <- ggplot(um_gr, aes(umap_1, umap_2, colour = gr_plot)) +
  geom_point(size = 0.42, stroke = 0) +
  scale_colour_gradientn(
    colours = c("grey87", "#d7b5a6", "#a63603", "#5c1404"),
    values  = rescale(c(0, 0.35, 1.3, gr_cap)),
    limits  = c(0, gr_cap),
    breaks  = c(0, 1, 2),
    name    = NULL
  ) +
  # Dashed box around cluster 17 + connector lines to inset
  annotate("rect", xmin = zbox$xmin, xmax = zbox$xmax, ymin = zbox$ymin, ymax = zbox$ymax,
           fill = NA, colour = pal17["17"], linewidth = 0.48, linetype = "dashed") +
  geom_segment(data = conn_df, aes(x = x, y = y, xend = xend, yend = yend),
               inherit.aes = FALSE, colour = pal17["17"], linewidth = 0.35, linetype = "dashed") +
  # Cluster badges 13, 9, 17
  geom_point(data = badges_d, aes(bx, by), inherit.aes = FALSE,
             size = 3.8, colour = badges_d$col) +
  geom_text(data = badges_d, aes(bx, by, label = k), inherit.aes = FALSE,
            size = 2.0, colour = badges_d$tcol, fontface = "bold", family = "Arial") +
  coord_cartesian(xlim = c(-17, 15.5), ylim = c(-13.5, 16)) +
  labs(title = expression("Expression of " * italic("Gr28b")), x = "umap_1", y = "umap_2") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 8.2, hjust = 0.5, margin = margin(b = 4)),
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35),
    axis.text = element_text(size = 6.5, colour = "black"),
    axis.title = element_text(size = 7.2, colour = "black"),
    legend.position = c(0.93, 0.87),
    legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.width = unit(2.5, "mm"),
    legend.key.height = unit(4.0, "mm"),
    legend.text = element_text(size = 5.8),
    plot.margin = margin(3, 4, 2, 2)
  )

# Magnified inset plot of cluster 17 (all 34 nuclei, 11 Gr28b+ = 32%)
xr <- range(c17$umap_1); yr <- range(c17$umap_2)
xspan <- diff(xr); yspan <- diff(yr)
pd_inset <- ggplot(c17, aes(umap_1, umap_2, fill = gr_plot)) +
  geom_point(shape = 21, size = 1.65, colour = "grey25", stroke = 0.22) +
  scale_fill_gradientn(
    colours = c("grey87", "#d7b5a6", "#a63603", "#5c1404"),
    values  = rescale(c(0, 0.35, 1.3, gr_cap)),
    limits  = c(0, gr_cap),
    guide   = "none"
  ) +
  annotate("text", x = mean(xr), y = yr[2] + yspan * 0.42,
           label = "17: 32% (11/34)", hjust = 0.5, vjust = 0.5,
           size = 1.60, fontface = "bold", colour = "#b2182b", family = "Arial") +
  coord_cartesian(
    xlim = c(xr[1] - xspan * 0.32, xr[2] + xspan * 0.32),
    ylim = c(yr[1] - yspan * 0.18, yr[2] + yspan * 0.65)
  ) +
  theme_void(base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "white", colour = pal17["17"], linewidth = 0.55),
    plot.margin = margin(0, 0, 0, 0)
  )

pd_full <- pd_main +
  annotation_custom(
    grob = ggplotGrob(pd_inset),
    xmin = ins$xmin, xmax = ins$xmax, ymin = ins$ymin, ymax = ins$ymax
  )

# ------------------------------------------------------------------------------
# 5. Panel e: Cross-species, multi-tissue Gr28b expression
#    Uses the Fig. 1c gene-level TPM table (gene_expression_wide_output.txt) with
#    verified species labels (fixing the original tau/cucurbitae label swap).
#    Generates BOTH:
#      (a) pe_bar  -> original-style 5-species TPM replicate dot + mean-bar plot (Main Fig. 2e)
#      (b) pe_heat -> log2(TPM + 1) heatmap with mashr stats (Response to Reviewers only)
# ------------------------------------------------------------------------------
G_ID <- "Bd1G02958"
M_F1C <- fig1c_tpm(ss_all)
G_F1C <- grep("^Bd1G02958_", rownames(M_F1C), value = TRUE)[1]
stopifnot(!is.na(G_F1C))
TPM_E <- M_F1C[G_F1C, ]

SPP_E <- c(Bc = "B. correcta", Bd = "B. dorsalis", Zc = "B. cucurbitae", Zt = "B. tau", Bm = "B. minax")
SPP_BAR_ORD <- c("Bc", "Bd", "Bm", "Zc", "Zt")  # exact same top-to-bottom order as Panel a
SP_COL_E <- SP_COLORS_A[names(SPP_E)]
cols_e <- rbind(
  CJ(tissue = TIS, group = c("FN", "FM"), sorted = FALSE),
  data.table(tissue = c("larval_brain", "larval_mouthpart"), group = "LA")
)

d_e <- rbindlist(lapply(names(SPP_E), function(sp_code) {
  cols_e[, {
    k <- ss_all[sp == sp_code & tissue == .BY$tissue & group == .BY$group, sample]
    v <- as.numeric(TPM_E[k])
    .(sp = sp_code, mean_tpm = mean(v), sd_tpm = sd(v), log2 = mean(log2(v + 1)), n = length(k))
  }, by = .(tissue, group)]
}))
d_e[, `:=`(mark = "", log2FC = NA_real_, lfsr = NA_real_)]
for (ti in TIS) {
  pm <- fread(RES("mashr", "maturation", sprintf("patterns_%s.tsv", ti)))
  lf <- fread(RES("mashr", "maturation", sprintf("lfsr_%s.tsv", ti)))
  ms <- fread(RES("mashr", "maturation", sprintf("missing_%s.tsv", ti)))
  for (sp_code in names(SPP_E)) {
    idx <- d_e$sp == sp_code & d_e$tissue == ti & d_e$group == "FM"
    if (!G_ID %in% pm$gene || ms[gene == G_ID][[sp_code]] == 1) {
      d_e[idx, mark := "n.t."]
      next
    }
    b <- pm[gene == G_ID][[sp_code]]; q <- lf[gene == G_ID][[sp_code]]
    d_e[idx, `:=`(log2FC = b, lfsr = q)]
    if (q < 0.05 && abs(b) > 1) d_e[idx, mark := if (q < 0.001) "***" else if (q < 0.01) "**" else "*"]
  }
}

# Sync SuppData7 (and legacy reproduce/results/tables/SuppData2) with exact plotted TPM values
fwrite(d_e[, .(tissue, group, sp, mean_tpm, sd_tpm, log2, n, mark, log2FC, lfsr)],
       file.path(REPO_ROOT, "data", "supplementary_data", "SuppData7_Gr28b_expression_and_tests.tsv"),
       sep = "\t")
fwrite(d_e[, .(tissue, group, sp, mean_tpm, sd_tpm, log2, n, mark, log2FC, lfsr)],
       RES("tables", "SuppData2_Gr28b_expression_and_tests.tsv"),
       sep = "\t")

d_e[, x := factor(paste(tissue, group), levels = paste(cols_e$tissue, cols_e$group))]
d_e[, y := factor(SPP_E[sp], levels = rev(SPP_E))]
xl_e <- ifelse(cols_e$group == "FN", "nF",
        ifelse(cols_e$group == "FM", "mF",
        ifelse(cols_e$tissue == "larval_brain", "La-Br", "La-Mo")))

# 5a. Heatmap version (Response to Reviewers comparison figure)
pe_heat <- ggplot(d_e, aes(x, y, fill = log2)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = mark, colour = log2 > 3.2), size = 2.35, fontface = "plain") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_fill_gradientn(
    colours = c("#faf2ee", "#d7b5a6", "#a86c66", "#671919"),
    limits  = c(0, 5.5),
    breaks  = c(0, 2, 4),
    name    = expression(log[2]*"(TPM + 1)")
  ) +
  scale_x_discrete(labels = xl_e, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  geom_vline(xintercept = seq(2.5, 12.5, 2), colour = "white", linewidth = 1.4) +
  annotate("text", x = seq(1.5, 11.5, 2), y = 0.16,
           label = tools::toTitleCase(TIS), size = 2.35, vjust = 1, colour = "black") +
  coord_cartesian(ylim = c(0.5, 5.5), clip = "off") +
  labs(title = expression("Expression of " * italic("Gr28b")), x = NULL, y = NULL) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 8.2, hjust = 0.5, margin = margin(b = 4)),
    axis.line  = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 5.6, colour = "black", margin = margin(t = 2)),
    axis.text.y = element_text(size = 7.2, face = "bold.italic",
                               colour = rev(unname(SP_COL_E)), margin = margin(r = 4)),
    legend.title = element_text(size = 6.0),
    legend.text  = element_text(size = 5.8),
    legend.key.width  = unit(2.8, "mm"),
    legend.key.height = unit(4.0, "mm"),
    plot.margin = margin(3, 4, 16, 4)
  )

# 5b. Original-style 5-species replicate dot + horizontal mean-bar plot (Main Fig. 2e)
# Paired x positions: 6 adult tissue pairs (gap between tissues) + 2 larval tissues
x_pos_map <- data.table(
  tissue = cols_e$tissue,
  group  = cols_e$group,
  xpos   = c(1.0, 2.0,   3.5, 4.5,   6.0, 7.0,   8.5, 9.5,   11.0, 12.0,   13.5, 14.5,   16.0, 17.2),
  xlab   = xl_e
)

d_rep <- rbindlist(lapply(SPP_BAR_ORD, function(sp_code) {
  x_pos_map[, {
    k <- ss_all[sp == sp_code & tissue == .BY$tissue & group == .BY$group, sample]
    v <- as.numeric(TPM_E[k])
    .(sp = sp_code, tpm = v, rep_idx = seq_along(v))
  }, by = .(tissue, group, xpos, xlab)]
}))
d_rep[, sp_fac := factor(SPP_E[sp], levels = SPP_E[SPP_BAR_ORD])]
# Slight deterministic horizontal offset for the 3 biological replicates so none overlap
d_rep[, x_jit := xpos + c(-0.14, 0, 0.14)[rep_idx]]

d_mean <- d_rep[, .(mean_tpm = mean(tpm)), by = .(sp, sp_fac, tissue, group, xpos)]
d_sp_lbl <- d_rep[, .(ymax = max(tpm) * 1.14), by = .(sp, sp_fac)]
d_sp_lbl[, col := unname(SP_COLORS_A[sp])]

# Clean per-species y-axis breaks (generous clearance from upper/lower facet borders)
clean_breaks <- function(lims) {
  hi <- lims[2]
  if (hi >= 35) c(0, 20, 40)
  else if (hi >= 18) c(0, 10, 20)
  else if (hi >= 10) c(0, 4, 8, 12)
  else if (hi >= 5.5) c(0, 2, 4, 6)
  else c(0, 1, 2, 3)
}

tis_centers <- data.table(
  x = c(1.5, 4.0, 6.5, 9.0, 11.5, 14.0, 16.6),
  label = c(tools::toTitleCase(TIS), "Larva")
)

# Subtle tissue group separator x-positions
tis_seps <- c(2.75, 5.25, 7.75, 10.25, 12.75, 15.25)

pe_bar_core <- ggplot() +
  geom_vline(xintercept = tis_seps, colour = "#ececec", linewidth = 0.35) +
  geom_point(data = d_rep, aes(x = x_jit, y = tpm),
             colour = "#8c8c8c", size = 1.05, alpha = 0.9, stroke = 0) +
  geom_segment(data = d_mean,
               aes(x = xpos - 0.36, xend = xpos + 0.36, y = mean_tpm, yend = mean_tpm),
               colour = "#111111", linewidth = 0.52) +
  geom_text(data = d_sp_lbl,
            aes(x = 17.65, y = ymax * 0.95, label = sp_fac, colour = sp),
            fontface = "bold.italic", family = "Arial", size = 2.55, hjust = 1, vjust = 1) +
  scale_colour_manual(values = SP_COLORS_A, guide = "none") +
  scale_x_continuous(
    limits = c(0.3, 17.9),
    breaks = x_pos_map$xpos,
    labels = x_pos_map$xlab,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = clean_breaks,
    expand = expansion(mult = c(0.06, 0.22))
  ) +
  facet_grid(sp_fac ~ ., scales = "free_y") +
  labs(title = expression("Expression of " * italic("Gr28b")), x = NULL, y = "Expression (TPM)") +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(
    plot.title       = element_text(size = 8.2, family = "Arial", hjust = 0.5, margin = margin(b = 4)),
    panel.grid       = element_blank(),
    panel.border     = element_rect(colour = "#333333", fill = NA, linewidth = 0.4),
    panel.spacing.y  = unit(0, "pt"),
    strip.background = element_blank(),
    strip.text       = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(linewidth = 0.35, colour = "#222222"),
    axis.ticks.length.y = unit(2.0, "pt"),
    axis.text.x      = element_text(size = 6.0, family = "Arial", colour = "#111111", margin = margin(t = 2.5)),
    axis.text.y      = element_text(size = 6.4, family = "Arial", colour = "#111111"),
    axis.title.y     = element_text(size = 7.2, family = "Arial", colour = "#111111", margin = margin(r = 3)),
    plot.margin      = margin(3, 4, 1, 4)
  )

pe_bar_tis <- ggplot(tis_centers, aes(x = x, y = 0.55, label = label)) +
  geom_text(family = "Arial", size = 2.38, colour = "#111111", vjust = 1) +
  scale_x_continuous(limits = c(0.3, 17.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 7, base_family = "Arial") +
  theme(
    panel.grid   = element_blank(),
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    plot.margin  = margin(0, 4, 2, 4)
  )

pe_bar <- (pe_bar_core / pe_bar_tis) + plot_layout(heights = c(1, 0.09))

# ------------------------------------------------------------------------------
# 6. Assemble pure-vector 5-panel Figure 2 (A4 portrait, 8.27 x 11.2 in)
#    Main Figure 2 uses pe_bar (original dot + mean-bar plot) in Panel e;
#    Standalone pe_bar and pe_heat are BOTH exported to reproduce/figures/.
# ------------------------------------------------------------------------------
draw_panel_letter <- function(letter, x_npc, y_npc) {
  grid.text(letter, x = unit(x_npc, "npc"), y = unit(y_npc, "npc"),
            just = c("left", "top"),
            gp = gpar(fontsize = 11.5, fontface = "bold", fontfamily = "Arial"))
}

render_full <- function(panel_e_plot = pe_bar) {
  grid.newpage()
  pushViewport(viewport(x = 0.5, y = 0.5, width = 0.96, height = 0.96))

  # Row 1: Panel a (top 40.5% of height)
  pushViewport(viewport(x = 0.02, y = 0.995, width = 0.96, height = 0.405, just = c("left", "top")))
  draw_panel_a()
  popViewport()
  draw_panel_letter("a", 0.002, 0.998)

  # Row 2: Panel b (UMAP + legend) and Panel c (DotPlot)
  pushViewport(viewport(x = 0.008, y = 0.578, width = 0.340, height = 0.275, just = c("left", "top")))
  grid.draw(ggplotGrob(pb_umap))
  popViewport()

  pushViewport(viewport(x = 0.348, y = 0.578, width = 0.285, height = 0.275, just = c("left", "top")))
  grid.draw(ggplotGrob(pb_leg))
  popViewport()

  pushViewport(viewport(x = 0.638, y = 0.586, width = 0.360, height = 0.292, just = c("left", "top")))
  grid.draw(ggplotGrob(pc_dot))
  popViewport()

  draw_panel_letter("b", 0.002, 0.588)
  draw_panel_letter("c", 0.630, 0.588)

  # Row 3: Panel d (Gr28b UMAP + cluster 17 inset) and Panel e (Gr28b expression across species)
  pushViewport(viewport(x = 0.008, y = 0.282, width = 0.375, height = 0.278, just = c("left", "top")))
  grid.draw(ggplotGrob(pd_full))
  popViewport()

  pushViewport(viewport(x = 0.390, y = 0.286, width = 0.605, height = 0.284, just = c("left", "top")))
  grid.draw(patchwork::patchworkGrob(panel_e_plot))
  popViewport()

  draw_panel_letter("d", 0.002, 0.290)
  draw_panel_letter("e", 0.388, 0.290)

  popViewport()
}

render_be <- function(panel_e_plot = pe_bar) {
  grid.newpage()
  pushViewport(viewport(x = 0.5, y = 0.5, width = 0.96, height = 0.96))

  # Top half: b + c
  pushViewport(viewport(x = 0.008, y = 0.985, width = 0.340, height = 0.475, just = c("left", "top")))
  grid.draw(ggplotGrob(pb_umap))
  popViewport()

  pushViewport(viewport(x = 0.348, y = 0.985, width = 0.285, height = 0.475, just = c("left", "top")))
  grid.draw(ggplotGrob(pb_leg))
  popViewport()

  pushViewport(viewport(x = 0.638, y = 0.995, width = 0.360, height = 0.495, just = c("left", "top")))
  grid.draw(ggplotGrob(pc_dot))
  popViewport()

  draw_panel_letter("b", 0.002, 0.998)
  draw_panel_letter("c", 0.630, 0.998)

  # Bottom half: d + e
  pushViewport(viewport(x = 0.008, y = 0.485, width = 0.375, height = 0.475, just = c("left", "top")))
  grid.draw(ggplotGrob(pd_full))
  popViewport()

  pushViewport(viewport(x = 0.390, y = 0.490, width = 0.605, height = 0.485, just = c("left", "top")))
  grid.draw(patchwork::patchworkGrob(panel_e_plot))
  popViewport()

  draw_panel_letter("d", 0.002, 0.495)
  draw_panel_letter("e", 0.388, 0.495)

  popViewport()
}

FIG_DIR <- OUT_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
PDF_OUT <- file.path(OUT_DIR, "Figure2.pdf")
PNG_OUT <- file.path(OUT_DIR, "Figure2.png")
BE_PNG  <- file.path(OUT_DIR, "Figure2_be.png")

# Main Figure 2 (with original dot + mean-bar plot in panel e)
cairo_pdf(PDF_OUT, width = 8.27, height = 11.2, family = "Arial")
render_full(pe_bar)
dev.off()

png(PNG_OUT, width = 8.27, height = 11.2, units = "in", res = 300, type = "cairo", family = "Arial")
render_full(pe_bar)
dev.off()

png(BE_PNG, width = 8.27, height = 6.5, units = "in", res = 300, type = "cairo", family = "Arial")
render_be(pe_bar)
dev.off()

# Export BOTH standalone Panel e versions (barplot + heatmap) into reproduce/figures/
for (stem in c("Fig2e_Gr28b_barplot")) {
  cairo_pdf(file.path(FIG_DIR, paste0(stem, ".pdf")), width = 5.4, height = 3.2, family = "Arial")
  print(pe_bar)
  dev.off()
  png(file.path(FIG_DIR, paste0(stem, ".png")), width = 5.4, height = 3.2, units = "in", res = 300, type = "cairo", family = "Arial")
  print(pe_bar)
  dev.off()
}
for (stem in c("Fig2e_Gr28b_heatmap", "Fig2d_Gr28b_heatmap")) {
  cairo_pdf(file.path(FIG_DIR, paste0(stem, ".pdf")), width = 5.6, height = 2.5, family = "Arial")
  print(pe_heat)
  dev.off()
  png(file.path(FIG_DIR, paste0(stem, ".png")), width = 5.6, height = 2.5, units = "in", res = 300, type = "cairo", family = "Arial")
  print(pe_heat)
  dev.off()
}

cat("Wrote:\n  ", PDF_OUT, "\n  ", PNG_OUT, "\n  ", BE_PNG,
    "\n  Fig2e_Gr28b_barplot.{pdf,png}\n  Fig2e_Gr28b_heatmap.{pdf,png}\n")
