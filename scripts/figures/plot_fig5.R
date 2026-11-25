#!/usr/bin/env Rscript
# Pure-source-code 4-panel Main Fig. 5 (= reproduce/figures/Fig4):
#   Streamlined, enlarged-font publication layout (verbose details moved to Fig. 5 legend).
#   a, b : Hi-C + TADs, Isoforms (Gr28b-X1 highlighted in orange), ATAC-seq MO peaks,
#          Ovipositor RNA-seq coverage (M1-M3, N1-N3), and 5-species Synteny.
#   c    : B. dorsalis MO ATAC peaks (P1-P10), Gr28b gene structure, 5-species % identity
#          heatmap + Locus median + B. minax ATAC accessibility row.
#   d    : High-resolution Jalview-style MSA of disrupted generalist-vs-specialist TFBS motifs
#          in P2 (pdm2 / nub) and P3 (tin; HLH4C / Fer2 / wor).
# 100% embedded TrueType Arial + Liberation Mono.

REPO_ROOT   <- Sys.getenv("REPO_ROOT", ".")
REVIEW_ROOT <- Sys.getenv("REVIEW_ROOT", "/path/to/project")
OUT_DIR     <- Sys.getenv("OUT_DIR", file.path(REPO_ROOT, "outputs", "figures"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
REPRO_ROOT  <- Sys.getenv("REPRO_ROOT", file.path(REVIEW_ROOT, "reproduce"))
Sys.setenv(REPRO_ROOT = normalizePath(REPRO_ROOT))
source(file.path(REPRO_ROOT, "R", "00_common.R"))


suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(rhdf5)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(grid)
})

H    <- CFG$hic_root
ODIR <- OUT_DIR

SP <- list(
  Bmin = list(name = "B. minax", col = "#1b7a3a", chrom = "BM1", gr = c(9477063, 9508221), wide = c(8290000, 9820000),
              zoom = c(9452500, 9532500), cool = "5k/hic_corrected.5k.cool", tad = "10k/TAD_10k_domains.bed",
              gtf = "Bmin.gtf", gene = "Bm1G00667",
              peaks = "Bm_MO.pval0.01.300K.narrowPeak",
              rna = c(paste0("Bm-FM-O", 1:3, ".bigwig"), paste0("Bm-FN-O", 1:3, ".bigwig")), synteny = NULL),
  Bdor = list(name = "B. dorsalis", col = "#c62828", chrom = "NC_064303.1", gr = c(108169689, 108194489), wide = c(106510000, 108510000),
              zoom = c(108142000, 108222000), cool = "10k/hic_corrected.10k.cool", tad = "10k/TAD_10k_domains.bed",
              gtf = "Bdor.gtf", gene = "Bd1G02958",
              peaks = "Bd_MO.pval0.01.300K.narrowPeak",
              rna = c(paste0("Bd-FM-O", 1:3, ".bigwig"), paste0("Bd-FN-O", 1:3, ".bigwig")),
              synteny = c(correcta = "Bdor_vs_Bcor", minax = "Bdor_vs_Bmin", cucurbitae = "Bdor_vs_Bcuc", tau = "Bdor_vs_Btau"))
)
SYN_COL <- c(correcta = "#e8862a", minax = "#1b7a3a", cucurbitae = "#b84fa0", tau = "#8c5fb8")
th <- theme_classic(base_size = 8.4, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35),
        axis.ticks = element_line(linewidth = 0.35),
        axis.text = element_text(colour = "#111111"),
        plot.margin = margin(1, 3, 1, 3))
blank_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.line.x = element_blank(), axis.title.x = element_blank())
mb <- function(x) x / 1e6

# Bd1G02958.1 (Gr28b-X1) coordinates
LD       <- c(108181355, 108181424)
EX1B     <- 108181492

make_peak_bells <- function(pk_dt, x_origin = 0, x_scale = 1, width_adjust = 1.0, min_half_bp = 0) {
  if (!nrow(pk_dt)) return(data.table(x = numeric(0), y = numeric(0), id = integer(0)))
  rbindlist(lapply(seq_len(nrow(pk_dt)), function(i) {
    a0 <- pk_dt$V2[i] + 1; b0 <- pk_dt$V3[i]
    sx <- a0 + pk_dt$V10[i]
    h  <- pk_dt$V7[i]
    hl <- max((sx - a0) * width_adjust, min_half_bp)
    hr <- max((b0 - sx) * width_adjust, min_half_bp)
    a  <- sx - hl; b <- sx + hr
    xs <- seq(a, b, length.out = 70)
    y  <- ifelse(xs <= sx,
                 h * (0.5 - 0.5 * cos(pi * (xs - a) / max(sx - a, 1))),
                 h * (0.5 + 0.5 * cos(pi * (xs - sx) / max(b - sx, 1))))
    data.table(x = (c(a, xs, b) - x_origin) / x_scale, y = c(0, y, 0), id = i)
  }))
}

# ------------------------------------------------------------------------------
# 1. Panels a & b
# ------------------------------------------------------------------------------
read_cool <- function(f, chrom, lo, hi) {
  names <- h5read(f, "chroms/name"); ci <- match(chrom, names)
  coff <- h5read(f, "indexes/chrom_offset"); b0 <- coff[ci]; b1 <- coff[ci + 1]
  st <- h5read(f, "bins/start", index = list((b0 + 1):b1)); w <- h5read(f, "bins/weight", index = list((b0 + 1):b1))
  bs <- st[2] - st[1]; keep <- which(st + bs > lo & st < hi); i0 <- b0 + min(keep) - 1; i1 <- b0 + max(keep) - 1
  boff <- h5read(f, "indexes/bin1_offset", index = list((i0 + 1):(i1 + 2))); p0 <- boff[1]; p1 <- boff[length(boff)]
  b1i <- h5read(f, "pixels/bin1_id", index = list((p0 + 1):p1)); b2i <- h5read(f, "pixels/bin2_id", index = list((p0 + 1):p1))
  cnt <- h5read(f, "pixels/count", index = list((p0 + 1):p1)); h5closeAll()
  ok <- b2i >= i0 & b2i <= i1
  d <- data.table(i = b1i[ok] - b0 + 1, j = b2i[ok] - b0 + 1, v = as.numeric(cnt[ok]))
  d[, v := v * w[i] * w[j]]; d <- d[is.finite(v)]
  d[, `:=`(x1 = st[i], x2 = st[j])][, .(x1, x2, v, bs = bs)]
}

hic_panel <- function(d, tads, lo, hi, depth) {
  bs <- d$bs[1]; d <- d[x2 - x1 <= depth]
  d[, `:=`(cx = (x1 + x2) / 2 + bs / 2, cy = (x2 - x1) / 2)]
  d[, id := .I]
  poly <- d[, .(x = c(cx - bs / 2, cx, cx + bs / 2, cx), y = c(cy, cy + bs / 2, cy, cy - bs / 2)), by = .(id, v)]
  poly[, y := pmax(y, 0)]
  tads[, tid := .I]
  tri <- tads[, .(x = c(start, (start + end) / 2, end), y = c(0, (end - start) / 2, 0)), by = tid]
  ggplot() +
    geom_polygon(data = poly, aes(x, y, group = id, fill = log1p(v)), colour = NA) +
    geom_path(data = tri, aes(x, y, group = tid), linewidth = 0.38, colour = "black") +
    scale_fill_gradientn(
      colours = rev(c("#a50026", "#f46d43", "#fee090", "#e0f3f8", "#74add1", "#313695")),
      limits = c(0, log1p(100)), oob = scales::squish,
      breaks = log1p(c(1, 5, 20, 100)), labels = c(1, 5, 20, 100), name = "Contacts"
    ) +
    coord_cartesian(xlim = c(lo, hi), ylim = c(0, depth / 2), expand = FALSE) + th +
    scale_x_continuous(labels = function(x) sprintf("%.1f", mb(x)), name = "Chr 1 (Mb)") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(size = 7.8),
          axis.title.x = element_text(size = 8.2),
          legend.key.height = unit(0.28, "cm"), legend.key.width = unit(0.18, "cm"),
          legend.title = element_text(size = 7.0), legend.text = element_text(size = 6.8))
}

SP_rna_max <- 1.1 * max(sapply(names(SP), function(sp) {
  s <- SP[[sp]]; R <- file.path(H, sp, "results_hicexplorer")
  max(sapply(s$rna, function(f) max(as.numeric(import(BigWigFile(file.path(R, f)),
    which = GRanges(s$chrom, IRanges(s$gr[1], s$gr[2])), as = "NumericList")[[1]]), na.rm = TRUE)))
}))

ATAC_YMAX <- 4.30

panel_ab <- function(sp) {
  s <- SP[[sp]]; R <- file.path(H, sp, "results_hicexplorer"); lo <- s$wide[1]; hi <- s$wide[2]
  d <- read_cool(file.path(R, s$cool), s$chrom, lo, hi)
  d <- d[x1 >= lo & x2 <= hi]
  tads <- fread(file.path(R, s$tad), header = FALSE)[V1 == s$chrom & V2 >= lo & V3 <= hi, .(start = V2, end = V3)]
  p_hic <- hic_panel(d, tads, lo, hi, depth = hi - lo) +
    labs(title = bquote(bolditalic(.(s$name)) ~ bold("(Chr 1)"))) +
    theme(plot.title = element_text(size = 9.5, hjust = 0.5))

  zl <- s$zoom[1]; zh <- s$zoom[2]
  hz_light <- annotate("rect", xmin = s$gr[1], xmax = s$gr[2], ymin = -Inf, ymax = Inf, fill = "#6fcf97", alpha = 0.12)
  ex <- import(file.path(R, s$gtf))
  ex <- ex[as.character(seqnames(ex)) == s$chrom & ex$type == "exon" & grepl(s$gene, ex$transcript_id)]
  ed <- data.table(start = start(ex), end = end(ex), tx = ex$transcript_id)
  foc_tx <- switch(sp, Bdor = "Bd1G02958.1", Bmin = "Bm1G00667.3")
  foc_lb <- switch(sp, Bdor = "bolditalic('Gr28b-X1')", Bmin = "bold('X1 ortholog')")
  tx <- ed[, .(s = min(start), e = max(end)), by = tx]
  tx[, foc := tx == foc_tx]
  tx[, y := { o <- order(foc); r <- integer(.N); r[o] <- seq_len(.N); r }]
  tx[, lx := s - 0.018 * (zh - zl)]
  ed <- merge(ed, tx[, .(tx, y, foc)], by = "tx")

  p_gene <- ggplot() + hz_light +
    geom_segment(data = tx[foc == FALSE], aes(x = s, xend = e, y = y, yend = y), linewidth = 0.35, colour = "#8c96a0") +
    geom_segment(data = tx[foc == TRUE],  aes(x = s, xend = e, y = y, yend = y), linewidth = 0.65, colour = "#e6550d") +
    geom_rect(data = ed[foc == FALSE], aes(xmin = start, xmax = end, ymin = y - 0.35, ymax = y + 0.35), fill = "#8c96a0") +
    geom_rect(data = ed[foc == TRUE],  aes(xmin = start, xmax = end, ymin = y - 0.35, ymax = y + 0.35), fill = "#e6550d") +
    geom_text(data = tx[foc == TRUE],  aes(x = lx, y = y), label = foc_lb, parse = TRUE, colour = "#d94801", hjust = 1, size = 2.65, family = "Arial") +
    annotate("text", x = mean(s$gr), y = nrow(tx) + 1.15, label = "Gr28b", fontface = "bold.italic", size = 3.0, family = "Arial") +
    coord_cartesian(xlim = c(zl, zh), ylim = c(0.4, nrow(tx) + 1.75), expand = FALSE, clip = "off") +
    labs(y = "Isoforms") + th + blank_x +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          axis.title.y = element_text(size = 8.2))

  qw <- c(s$gr[1] - 5000, s$gr[2] + 5000)
  pk_uniq <- unique(fread(file.path(R, s$peaks), header = FALSE, select = 1:3))[
    V1 == s$chrom & V3 > qw[1] & V2 < qw[2], .(start = V2 + 1, end = V3)]
  n_gene <- nrow(pk_uniq[end > s$gr[1] & start < s$gr[2]])

  pkz <- fread(file.path(R, s$peaks), header = FALSE)[V1 == s$chrom & V3 > zl & V2 < zh]
  pk_poly <- make_peak_bells(pkz, width_adjust = 1.65, min_half_bp = 650)
  if (nrow(pk_poly)) pk_poly[, y := pmin(y, ATAC_YMAX)]

  p_atac <- ggplot() +
    geom_vline(xintercept = s$gr, linetype = "dashed", colour = "#66b288", linewidth = 0.3) +
    (if (sp == "Bdor") annotate("rect", xmin = EX1B - 2400, xmax = EX1B + 1400, ymin = 0, ymax = ATAC_YMAX * 0.84,
                                fill = "#9e6bc4", alpha = 0.18)) +
    (if (nrow(pk_poly)) geom_polygon(data = pk_poly, aes(x, y, group = id), fill = "#3b82c4", alpha = 0.62, colour = NA)) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey25") +
    (if (sp == "Bdor") annotate("text", x = EX1B - 500, y = ATAC_YMAX * 0.92, label = "P2–P4",
                                size = 2.65, fontface = "bold", colour = "#6a3d9a", family = "Arial")) +
    (if (sp == "Bdor") {
      annotate("text", x = zl + 1600, y = ATAC_YMAX * 0.94, hjust = 0, vjust = 1, size = 2.65,
               fontface = "bold", colour = "#1f4e8c", family = "Arial",
               label = sprintf("%d peaks", n_gene))
    } else {
      annotate("text", x = mean(s$gr), y = ATAC_YMAX * 0.88, hjust = 0.5, vjust = 1, size = 2.65,
               fontface = "bold", colour = "#1f4e8c", family = "Arial",
               label = sprintf("%d peaks", n_gene))
    }) +
    scale_y_continuous(breaks = c(0, 2, 4)) +
    coord_cartesian(xlim = c(zl, zh), ylim = c(0, ATAC_YMAX), expand = FALSE) +
    labs(y = "ATAC-seq") + th + blank_x +
    theme(axis.title.y = element_text(size = 8.2), axis.text.y = element_text(size = 7.4))

  bin <- function(bwf, chr, n = 600) {
    gr <- GRanges(chr, IRanges(zl, zh)); v <- import(BigWigFile(bwf), which = gr, as = "NumericList")[[1]]
    data.table(x = seq(zl, zh, length.out = length(v)), v = as.numeric(v))[
      , .(v = mean(v, na.rm = TRUE)), by = .(x = cut(x, n, labels = FALSE))][
      , x := zl + (x - 0.5) * (zh - zl) / n]
  }
  rn <- rbindlist(lapply(seq_along(s$rna), function(k) cbind(bin(file.path(R, s$rna[k]), s$chrom),
                  lab = c(paste0("Mature ", 1:3), paste0("New ", 1:3))[k])))
  rn[, lab := factor(lab, levels = c(paste0("Mature ", 1:3), paste0("New ", 1:3)))]
  ymax <- max(SP_rna_max, na.rm = TRUE)
  p_rna <- ggplot(rn) + hz_light + geom_area(aes(x, pmin(v, ymax)), fill = "#c0504d") +
    scale_y_continuous(breaks = round(ymax * 0.85, -2), expand = c(0, 0)) +
    facet_grid(lab ~ ., switch = "y") +
    coord_cartesian(xlim = c(zl, zh), ylim = c(0, ymax), expand = FALSE) +
    labs(y = "Ovipositor RNA-seq") + th +
    scale_x_continuous(labels = function(x) sprintf("%.3f", mb(x)),
                       name = "Chr 1, 80-kb zoom (Mb)", n.breaks = 4) +
    theme(strip.background = element_blank(), strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, size = 7.4, colour = "#111111"),
          panel.spacing = unit(0.8, "pt"),
          axis.text.y = element_text(size = 6.8),
          axis.text.x = element_text(size = 7.6),
          axis.title.x = element_text(size = 8.2),
          axis.title.y = element_text(size = 8.2))

  if (!is.null(s$synteny)) {
    sy <- rbindlist(lapply(names(s$synteny), function(k) {
      x <- fread(file.path(R, paste0(s$synteny[[k]], ".gui-gapped1000.genomic.bed")), header = FALSE, skip = 1)
      x[V1 == s$chrom & V3 > zl & V2 < zh, .(start = V2, end = V3, sp = k)]
    }))
    sy[, sp := factor(sp, levels = rev(names(s$synteny)))]
    yt_s <- c(Bdor = 5, Bcor = 4, Bmin = 3, Zcuc = 2, Ztau = 1)
    yA <- mean(yt_s[c("Bdor", "Bcor")]); yB <- yt_s["Bmin"]
    yD <- mean(yt_s[c("Zcuc", "Ztau")]); yC <- mean(c(yA, yB))
    ts_s <- rbind(
      data.frame(x = 0, xe = 0, y = min(yC, yD), ye = max(yC, yD)),
      data.frame(x = 0, xe = 1, y = yC, ye = yC),
      data.frame(x = 0, xe = 1, y = yD, ye = yD),
      data.frame(x = 1, xe = 1, y = min(yA, yB), ye = max(yA, yB)),
      data.frame(x = 1, xe = 2, y = yA, ye = yA),
      data.frame(x = 1, xe = 3, y = yB, ye = yB),
      data.frame(x = 1, xe = 1, y = min(yt_s["Zcuc"], yt_s["Ztau"]), ye = max(yt_s["Zcuc"], yt_s["Ztau"])),
      data.frame(x = 1, xe = 3, y = yt_s["Zcuc"], ye = yt_s["Zcuc"]),
      data.frame(x = 1, xe = 3, y = yt_s["Ztau"], ye = yt_s["Ztau"]),
      data.frame(x = 2, xe = 2, y = min(yt_s["Bdor"], yt_s["Bcor"]), ye = max(yt_s["Bdor"], yt_s["Bcor"])),
      data.frame(x = 2, xe = 3, y = yt_s["Bdor"], ye = yt_s["Bdor"]),
      data.frame(x = 2, xe = 3, y = yt_s["Bcor"], ye = yt_s["Bcor"])
    )
    tx0 <- zl - 11500; tstep <- 2200
    ts_s$x  <- tx0 + ts_s$x * tstep
    ts_s$xe <- tx0 + ts_s$xe * tstep
    lab_s <- data.frame(
      x = tx0 + 3 * tstep + 900,
      y = c(5, 4, 3, 2, 1),
      lab = c("B. dorsalis", "B. correcta", "B. minax", "B. cucurbitae", "B. tau"),
      ff  = c("bold.italic", "italic", "italic", "italic", "italic")
    )
    p_bot <- ggplot(sy) +
      geom_segment(data = ts_s, aes(x = x, y = y, xend = xe, yend = ye),
                   linewidth = 0.35, colour = "grey25", inherit.aes = FALSE) +
      geom_text(data = lab_s, aes(x = x, y = y, label = lab, fontface = ff),
                hjust = 0, size = 2.45, colour = "#111111", family = "Arial", inherit.aes = FALSE) +
      geom_rect(aes(xmin = pmax(start, zl), xmax = pmin(end, zh),
                    ymin = as.numeric(sp) - 0.35, ymax = as.numeric(sp) + 0.35, fill = sp)) +
      scale_fill_manual(values = SYN_COL, guide = "none") +
      coord_cartesian(xlim = c(zl, zh), ylim = c(0.4, 5.6), expand = FALSE, clip = "off") +
      labs(y = "Synteny") + th + blank_x +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
            axis.title.y = element_text(size = 8.0),
            plot.margin = margin(2, 3, 1, 3))
  } else {
    p_bot <- plot_spacer()
  }

  heights <- c(0, 0.85, 0.84, 2.25, 0.92)
  heights[1] <- 0.21 * sum(heights[-1]) / (1 - 0.21)
  wrap_plots(list(p_hic, p_gene, p_atac, p_rna, p_bot), ncol = 1, heights = heights)
}

# ------------------------------------------------------------------------------
# 2. Shared data for Panels c and d
# ------------------------------------------------------------------------------
TDIR  <- file.path(RES("locus", "tfbs"), "peak_anchored")
LATC  <- c(Bdor = "B. dorsalis", Bcor = "B. correcta", Bmin = "B. minax",
           Zcuc = "B. cucurbitae", Ztau = "B. tau")
SPORD <- c("Bdor", "Bcor", "Bmin", "Zcuc", "Ztau")
GTAG  <- c(Bdor = "G", Bcor = "G", Bmin = "S", Zcuc = "G", Ztau = "G")
GCOL  <- c(G = "#1b9e77", S = "#b34710")

hom  <- fread(file.path(TDIR, "homolog_table.tsv"))
bmat <- fread(file.path(TDIR, "bmin_atac_homologs.tsv"))
flid <- fread(file.path(TDIR, "fulllength_identity.tsv"))
bgm  <- fread(file.path(TDIR, "locus_background_summary.tsv"))
qbed <- fread(file.path(TDIR, "queries_flank100.bed"), header = FALSE)

tree_segs <- function(yt) {
  yA <- mean(yt[c("Bdor", "Bcor")]); yB <- yt["Bmin"]
  yD <- mean(yt[c("Zcuc", "Ztau")]); yC <- mean(c(yA, yB))
  rbind(
    data.frame(x = 0, xe = 0, y = min(yC, yD), ye = max(yC, yD)),
    data.frame(x = 0, xe = 1, y = yC, ye = yC),
    data.frame(x = 0, xe = 1, y = yD, ye = yD),
    data.frame(x = 1, xe = 1, y = min(yA, yB), ye = max(yA, yB)),
    data.frame(x = 1, xe = 2, y = yA, ye = yA),
    data.frame(x = 1, xe = 3, y = yB, ye = yB),
    data.frame(x = 1, xe = 1, y = min(yt["Zcuc"], yt["Ztau"]), ye = max(yt["Zcuc"], yt["Ztau"])),
    data.frame(x = 1, xe = 3, y = yt["Zcuc"], ye = yt["Zcuc"]),
    data.frame(x = 1, xe = 3, y = yt["Ztau"], ye = yt["Ztau"]),
    data.frame(x = 2, xe = 2, y = min(yt["Bdor"], yt["Bcor"]), ye = max(yt["Bdor"], yt["Bcor"])),
    data.frame(x = 2, xe = 3, y = yt["Bdor"], ye = yt["Bdor"]),
    data.frame(x = 2, xe = 3, y = yt["Bcor"], ye = yt["Bcor"])
  )
}

peaks <- qbed[grepl("^peak_", V4)][order(V7)]
peaks[, pid := paste0("P", .I)]
peaks[, mid_bp := (V7 + V8) / 2]
peaks[, mid_kb := (mid_bp - EX1B) / 1000]
peaks[, xlab := sprintf("%s\n%.1f", pid, mid_kb)]

# ------------------------------------------------------------------------------
# 3. Panel c: Streamlined peak + gene + conservation heatmap + B. minax ATAC
# ------------------------------------------------------------------------------
GX_LO <- EX1B - 6800
GX_HI <- EX1B + 18700
g2x <- function(bp) 0.55 + (bp - GX_LO) / (GX_HI - GX_LO) * (10.45 - 0.55)

R_bd <- file.path(H, "Bdor/results_hicexplorer")
pk_c_raw <- fread(file.path(R_bd, SP$Bdor$peaks), header = FALSE)[
  V1 == SP$Bdor$chrom & V3 > GX_LO & V2 < GX_HI]
pk_c_poly <- make_peak_bells(pk_c_raw, x_origin = 0, x_scale = 1, width_adjust = 1.05, min_half_bp = 155)
pk_c_poly[, gx := g2x(x)]
pk_c_raw[, id := .I]
pk_c_raw[, is_utr := (V2 + V3) / 2 >= 108179500 & (V2 + V3) / 2 <= 108182100]
pk_c_poly <- merge(pk_c_poly, pk_c_raw[, .(id, is_utr)], by = "id")

Y_ATAC0 <- 9.75
Y_ATAC_H <- 2.20
ATAC_C_MAX <- 4.10
pk_c_poly[, py := Y_ATAC0 + pmin(y, ATAC_C_MAX) / ATAC_C_MAX * Y_ATAC_H]

pk_lab <- copy(peaks)[, col_x := .I]
pk_lab[, gx := g2x(mid_bp)]
pk_lab[, max_v7 := sapply(seq_len(.N), function(i) {
  sub <- pk_c_raw[V3 > pk_lab$V7[i] & V2 < pk_lab$V8[i]]
  if (nrow(sub)) max(sub$V7) else 2.0
})]
pk_lab[, py_top := Y_ATAC0 + max_v7 / ATAC_C_MAX * Y_ATAC_H]
pk_lab[, lx := gx]
pk_lab[, ly := py_top + 0.16]
pk_lab[pid == "P3", `:=`(lx = gx - 0.08, ly = py_top + 0.42)]
pk_lab[pid == "P4", `:=`(lx = gx + 0.22, ly = py_top + 0.14)]
pk_lab[pid == "P5", `:=`(lx = gx - 0.08, ly = py_top + 0.16)]
pk_lab[pid == "P6", `:=`(lx = gx + 0.04, ly = py_top + 0.56)]
pk_lab[pid == "P7", `:=`(lx = gx + 0.24, ly = py_top + 0.16)]
pk_lab[, txt := pid]
pk_lab[, is_utr := pid %in% c("P2", "P3", "P4")]

gtf_bd <- import(file.path(R_bd, SP$Bdor$gtf))
g_sub  <- gtf_bd[as.character(seqnames(gtf_bd)) == SP$Bdor$chrom & end(gtf_bd) > GX_LO & start(gtf_bd) < GX_HI]
gdt    <- data.table(type = as.character(g_sub$type), tx = g_sub$transcript_id,
                     s = pmax(start(g_sub), GX_LO), e = pmin(end(g_sub), GX_HI))
grc_ex  <- gdt[tx == "Bd1G02958.1" & type == "exon"]
grc_cds <- gdt[tx == "Bd1G02958.1" & type == "CDS"]
alt_ex  <- gdt[tx %in% c("Bd1G02958.5", "Bd1G02958.2", "Bd1G02958.4") & type == "exon"]
grd_utr <- data.table(s = 108187719, e = 108187912)
nbr_ex  <- gdt[tx %in% c("Bd1G02960.1", "Bd1G02962.1") & type == "exon"]

conn_dt <- rbindlist(lapply(seq_len(nrow(pk_lab)), function(i) {
  x0 <- pk_lab$col_x[i]; x1 <- pk_lab$gx[i]
  t  <- seq(0, 1, length.out = 25)
  s  <- 3 * t^2 - 2 * t^3
  data.table(id = i, x = x0 + (x1 - x0) * s, y = 6.46 + (7.64 - 6.46) * t,
             pid = pk_lab$pid[i], is_utr = pk_lab$is_utr[i])
}))

kb_ticks <- data.table(kb = c(-5, 0, 5, 10, 15))
kb_ticks[, gx := g2x(EX1B + kb * 1000)]
kb_ticks[, lab := ifelse(kb == 0, "0", sprintf("%+d kb", kb))]

hm <- merge(flid[grepl("^peak_", query) & species %in% SPORD, .(query, species, fl_ident)],
            hom[grepl("^peak_", query) & species %in% SPORD, .(query, species, homolog)],
            by = c("query", "species"), all.x = TRUE)
hm[, pident := ifelse(homolog == "Y" & fl_ident != "", as.numeric(fl_ident), NA_real_)]
hm <- rbind(hm, data.table(query = peaks$V4, species = "Bdor", pident = 100), fill = TRUE)
hm[, x := match(query, peaks$V4)][, y := 7 - match(species, SPORD)]

BGX <- nrow(peaks) + 1.35
bgd <- rbind(data.table(species = "Bdor", med = 100),
             bgm[species %in% SPORD, .(species, med = as.numeric(median))])
bgd[, y := 7 - match(species, SPORD)]

bmr <- bmat[grepl("^peak_", query), .(query, ov = peak_overlap, has = bmin_homolog)]
bmr[, x := match(query, peaks$V4)][, y := 1]

yt_h <- setNames(7 - seq_along(SPORD), SPORD)
ts_h <- tree_segs(yt_h)
ts_h$x  <- -4.35 + ts_h$x * 0.36
ts_h$xe <- -4.35 + ts_h$xe * 0.36
tg_h <- data.frame(
  x = -2.92,
  y = yt_h[SPORD],
  g = GTAG[SPORD],
  gcol = GCOL[GTAG[SPORD]],
  lab = LATC[SPORD]
)

utr_gx1 <- g2x(108179450); utr_gx2 <- g2x(108182150)

panelC <- ggplot(hm) +
  annotate("rect", xmin = utr_gx1, xmax = utr_gx2, ymin = 7.66, ymax = Y_ATAC0 + Y_ATAC_H + 0.15,
           fill = "#b08cc9", alpha = 0.16) +
  annotate("rect", xmin = 1.50, xmax = 4.50, ymin = 0.48, ymax = 6.48,
           fill = "#f3eaf8", alpha = 0.45, colour = "#7b3294", linewidth = 0.55, linetype = "dashed") +
  # Layer 1: ATAC peak track
  geom_polygon(data = pk_c_poly[is_utr == FALSE], aes(gx, py, group = id), fill = "#3b82c4", alpha = 0.65, colour = NA) +
  geom_polygon(data = pk_c_poly[is_utr == TRUE],  aes(gx, py, group = id), fill = "#7b3294", alpha = 0.68, colour = NA) +
  annotate("segment", x = 0.55, xend = 10.45, y = Y_ATAC0, yend = Y_ATAC0, linewidth = 0.32, colour = "grey25") +
  annotate("segment", x = 0.55, xend = 0.55, y = Y_ATAC0, yend = Y_ATAC0 + 4 / ATAC_C_MAX * Y_ATAC_H, linewidth = 0.3, colour = "grey25") +
  annotate("segment", x = 0.46, xend = 0.55, y = Y_ATAC0 + c(0, 2, 4) / ATAC_C_MAX * Y_ATAC_H,
           yend = Y_ATAC0 + c(0, 2, 4) / ATAC_C_MAX * Y_ATAC_H, linewidth = 0.25, colour = "grey25") +
  annotate("text", x = 0.40, y = Y_ATAC0 + c(0, 2, 4) / ATAC_C_MAX * Y_ATAC_H, label = c("0", "2", "4"),
           hjust = 1, size = 2.4, colour = "grey25", family = "Arial") +
  annotate("text", x = -0.25, y = Y_ATAC0 + Y_ATAC_H * 0.50, label = "ATAC peaks",
           hjust = 1, size = 2.8, fontface = "bold", colour = "#111111", family = "Arial") +
  geom_text(data = pk_lab, aes(x = lx, y = ly, label = txt),
            size = 2.55, fontface = ifelse(pk_lab$is_utr, "bold", "plain"),
            colour = ifelse(pk_lab$is_utr, "#6a3d9a", "#1f4e8c"), family = "Arial", vjust = 0) +
  annotate("segment", x = pk_lab[pid == "P3"]$gx, xend = pk_lab[pid == "P3"]$gx,
           y = pk_lab[pid == "P3"]$py_top + 0.04, yend = pk_lab[pid == "P3"]$ly - 0.04,
           linewidth = 0.25, colour = "#6a3d9a") +
  annotate("segment", x = pk_lab[pid == "P6"]$gx, xend = pk_lab[pid == "P6"]$gx,
           y = pk_lab[pid == "P6"]$py_top + 0.04, yend = pk_lab[pid == "P6"]$ly - 0.04,
           linewidth = 0.25, colour = "#1f4e8c") +
  annotate("text", x = (utr_gx1 + utr_gx2) / 2, y = 12.70,
           label = "P2–P4 (promoter / 5'UTR)",
           size = 2.75, fontface = "bold", colour = "#6a3d9a", family = "Arial") +
  # Layer 2: Gene structure track
  annotate("segment", x = g2x(EX1B), xend = g2x(108194489), y = 8.95, yend = 8.95, linewidth = 0.45, colour = "#e6550d") +
  annotate("rect", xmin = g2x(LD[1]), xmax = g2x(LD[2]), ymin = 8.76, ymax = 9.14, fill = "#7f8c8d") +
  geom_rect(data = grc_ex, aes(xmin = g2x(s), xmax = g2x(e), ymin = 8.74, ymax = 9.16),
            fill = "#fdae6b", colour = "#d94801", linewidth = 0.22, inherit.aes = FALSE) +
  geom_rect(data = grc_cds, aes(xmin = g2x(s), xmax = g2x(e), ymin = 8.68, ymax = 9.22),
            fill = "#e6550d", colour = NA, inherit.aes = FALSE) +
  annotate("text", x = g2x(EX1B) - 0.14, y = 8.95, label = "bolditalic('Gr28b-X1')",
           parse = TRUE, hjust = 1, size = 2.65, colour = "#b34710", family = "Arial") +
  annotate("segment", x = g2x(GX_LO), xend = g2x(108194489), y = 8.12, yend = 8.12, linewidth = 0.30, colour = "#8c96a0") +
  geom_rect(data = alt_ex, aes(xmin = g2x(s), xmax = g2x(e), ymin = 7.95, ymax = 8.29),
            fill = "#8c96a0", colour = NA, inherit.aes = FALSE) +
  geom_rect(data = grd_utr, aes(xmin = g2x(s), xmax = g2x(e), ymin = 7.93, ymax = 8.31),
            fill = "#fdae6b", colour = "#e66101", linewidth = 0.25, inherit.aes = FALSE) +
  annotate("text", x = g2x(108187820), y = 8.53, label = "X5", size = 2.55, fontface = "bold", colour = "#d95f02", family = "Arial") +
  annotate("segment", x = g2x(108194907), xend = g2x(108197440), y = 8.12, yend = 8.12, linewidth = 0.28, colour = "grey55") +
  annotate("segment", x = g2x(108198490), xend = g2x(GX_HI),     y = 8.12, yend = 8.12, linewidth = 0.28, colour = "grey55") +
  geom_rect(data = nbr_ex, aes(xmin = g2x(s), xmax = g2x(e), ymin = 7.96, ymax = 8.28),
            fill = "#b8b8b8", colour = "grey50", linewidth = 0.18, inherit.aes = FALSE) +
  annotate("text", x = g2x(108197800), y = 8.55, label = "Neighbors", size = 2.45, colour = "grey40", family = "Arial") +
  annotate("text", x = -0.25, y = 8.55, label = "Gene",
           hjust = 1, size = 2.8, fontface = "bold", colour = "#111111", family = "Arial") +
  # Layer 3: S-curve connectors
  geom_path(data = conn_dt[is_utr == FALSE], aes(x, y, group = id), linewidth = 0.30, colour = "grey65", inherit.aes = FALSE) +
  geom_path(data = conn_dt[is_utr == TRUE],  aes(x, y, group = id), linewidth = 0.48, colour = "#7b3294", inherit.aes = FALSE) +
  annotate("segment", x = 0.55, xend = 10.45, y = 7.66, yend = 7.66, linewidth = 0.26, colour = "grey45") +
  geom_segment(data = kb_ticks, aes(x = gx, xend = gx, y = 7.66, yend = 7.52), linewidth = 0.26, colour = "grey45", inherit.aes = FALSE) +
  geom_label(data = kb_ticks, aes(x = gx, y = 7.45, label = lab), size = 2.35, colour = "grey30", family = "Arial",
             vjust = 1, fill = "white", linewidth = 0, label.padding = unit(0.06, "lines"), inherit.aes = FALSE) +
  # Layer 4: Heatmap + Locus median + B. minax ATAC
  geom_tile(aes(x, y, fill = pident), width = 0.90, height = 0.84, colour = "white", linewidth = 0.32) +
  geom_text(data = hm[!is.na(pident)], aes(x, y, label = sprintf("%.0f", pident), colour = pident >= 58),
            size = 2.50, fontface = "bold", family = "Arial", show.legend = FALSE) +
  geom_tile(data = bgd, aes(x = BGX, y, fill = med), width = 0.90, height = 0.84, colour = "white", linewidth = 0.32) +
  geom_text(data = bgd, aes(x = BGX, y, label = sprintf("%.0f", med), colour = med >= 58),
            size = 2.50, fontface = "bold", family = "Arial", show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#1a1a1a"), guide = "none") +
  geom_text(data = hm[is.na(pident)], aes(x, y), label = "×", size = 3.4, colour = "grey45", fontface = "bold", family = "Arial") +
  geom_point(data = bmr[has == "Y"], aes(x, y, shape = ov), size = 2.6, fill = "#1f4e8c", colour = "#1f4e8c", stroke = 0.7) +
  geom_text(data = bmr[has == "N"], aes(x, y), label = "×", size = 3.4, colour = "grey45", fontface = "bold", family = "Arial") +
  annotate("segment", x = 0.48, xend = nrow(peaks) + 0.52, y = 1.50, yend = 1.50, linewidth = 0.25, colour = "grey70") +
  # Cladogram -> G/S -> Species
  geom_segment(data = ts_h, aes(x, y, xend = xe, yend = ye), linewidth = 0.42, colour = "grey25") +
  geom_rect(data = tg_h, aes(xmin = x - 0.27, xmax = x + 0.27, ymin = y - 0.34, ymax = y + 0.34), fill = tg_h$gcol) +
  geom_text(data = tg_h, aes(x, y, label = g), size = 2.45, colour = "white", fontface = "bold", family = "Arial") +
  geom_text(data = tg_h, aes(x = 0.32, y = y, label = lab), hjust = 1, size = 2.75, fontface = "italic", colour = "#111111", family = "Arial") +
  annotate("text", x = 0.32, y = 1.0, label = "italic('B. minax')~'ATAC'", parse = TRUE, hjust = 1, size = 2.75, colour = "#111111", family = "Arial") +
  # Compact inline B. minax ATAC key (top right)
  annotate("point", x = 7.65, y = 12.70, shape = 16, size = 2.2, colour = "#1f4e8c") +
  annotate("text",  x = 7.90, y = 12.70, label = "peak", hjust = 0, size = 2.45, colour = "grey20", family = "Arial") +
  annotate("point", x = 9.05, y = 12.70, shape = 1, size = 2.2, colour = "#1f4e8c", stroke = 0.65) +
  annotate("text",  x = 9.30, y = 12.70, label = "no peak", hjust = 0, size = 2.45, colour = "grey20", family = "Arial") +
  annotate("text",  x = 10.95, y = 12.70, label = "×", hjust = 0.5, size = 3.0, fontface = "bold", colour = "grey45", family = "Arial") +
  annotate("text",  x = 11.18, y = 12.70, label = "absent", hjust = 0, size = 2.45, colour = "grey20", family = "Arial") +
  scale_shape_manual(values = c(Y = 16, N = 1), guide = "none") +
  scale_fill_gradientn(
    colours = c("#f7fbff", "#9ecae1", "#3182bd", "#084594"), limits = c(0, 100), na.value = "#eaeaea",
    name = "Identity (%)",
    guide = guide_colourbar(barheight = unit(36, "pt"), barwidth = unit(7, "pt"), title.position = "top")
  ) +
  scale_x_continuous(
    breaks = c(seq_len(nrow(peaks)), BGX),
    labels = c(peaks$xlab, "Locus\nmedian"),
    limits = c(-4.55, BGX + 1.35), expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0.35, 13.05), expand = c(0, 0)) +
  labs(x = "Distance to X1 5'UTR start (kb)", y = NULL) +
  th +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.3),
        axis.text.x = element_text(size = 7.4, lineheight = 0.92, colour = "#111111"),
        axis.title.x = element_text(size = 8.2, face = "bold", hjust = 0.62, margin = margin(t = 3)),
        legend.text = element_text(size = 7.2), legend.title = element_text(size = 7.6, face = "bold"),
        legend.position = "right", legend.margin = margin(0, 0, 0, -4),
        plot.margin = margin(2, 4, 2, 2))

# ------------------------------------------------------------------------------
# 4. Panel d: Clean, large-font MSA (coordinates & FIMO details moved to legend!)
# ------------------------------------------------------------------------------
load_mafft <- function(fa_path) {
  aln <- scan(fa_path, what = "", quiet = TRUE)
  hd  <- grepl("^>", aln)
  sq  <- setNames(toupper(tapply(aln[!hd], cumsum(hd)[!hd], paste0, collapse = "")), gsub(">", "", aln[hd]))
  sch <- lapply(sq, function(s) strsplit(s, "")[[1]])
  list(sch = sch, bdp = cumsum(sch[["Bdor"]] != "-"))
}

aln_p2 <- load_mafft(file.path(TDIR, "P2_chained.mafft.fa"))
aln_p3 <- load_mafft(file.path(TDIR, "E_homologs.mafft.fa"))

make_motif_panel <- function(aln_obj, cols, motif_cols, title_label, hit_clades, show_tree = TRUE) {
  sch <- aln_obj$sch
  ref <- sch[["Bdor"]][cols]
  n_sp <- length(SPORD)
  dt <- rbindlist(lapply(seq_along(SPORD), function(i) {
    sp <- SPORD[i]
    v  <- sch[[sp]][cols]
    st <- ifelse(v == "-", "gap", ifelse(sp == "Bdor" | v == ref, "same", "sub"))
    data.table(x = seq_along(cols), y = n_sp + 1 - i, sp = sp, ch = v, st = st)
  }))
  mx1 <- which(cols == motif_cols[1]) - 0.5
  mx2 <- which(cols == motif_cols[2]) + 0.5
  bx <- rbindlist(lapply(hit_clades, function(cl) {
    ys <- n_sp + 1 - match(cl, SPORD)
    data.table(x1 = mx1, x2 = mx2, y1 = min(ys) - 0.48, y2 = max(ys) + 0.48)
  }))

  x_min <- if (show_tree) -12.6 else 0.3
  p <- ggplot(dt) +
    geom_tile(aes(x, y, fill = st), width = 0.94, height = 0.86, colour = "white", linewidth = 0.25) +
    scale_fill_manual(values = c(same = "#dceaf7", sub = "#fddbc7", gap = "#f2f2f2"), guide = "none") +
    geom_text(aes(x, y, label = ch, colour = st, fontface = ifelse(st == "sub", "bold", "plain")),
              family = "Liberation Mono", size = 2.55, show.legend = FALSE) +
    scale_colour_manual(values = c(same = "#1a1a1a", sub = "#b2182b", gap = "grey55"), guide = "none") +
    geom_rect(data = bx, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
              fill = NA, colour = "#08519c", linewidth = 0.72) +
    annotate("segment", x = (mx1 + mx2) / 2, xend = (mx1 + mx2) / 2, y = n_sp + 0.52, yend = n_sp + 0.92,
             linewidth = 0.45, colour = "#08519c") +
    annotate("text", x = (mx1 + mx2) / 2, y = n_sp + 1.25, label = title_label,
             size = 2.85, fontface = "bold", colour = "#08519c", family = "Arial")

  if (show_tree) {
    yt_e <- setNames(n_sp + 1 - seq_along(SPORD), SPORD)
    ts_e <- tree_segs(yt_e)
    ts_e$x  <- -12.2 + ts_e$x * 0.95
    ts_e$xe <- -12.2 + ts_e$xe * 0.95
    tg_e <- data.frame(x = -8.5, y = yt_e[SPORD], g = GTAG[SPORD], gcol = GCOL[GTAG[SPORD]], lab = LATC[SPORD])
    p <- p +
      geom_segment(data = ts_e, aes(x, y, xend = xe, yend = ye), linewidth = 0.42, colour = "grey25") +
      geom_rect(data = tg_e, aes(xmin = x - 0.68, xmax = x + 0.68, ymin = y - 0.36, ymax = y + 0.36), fill = tg_e$gcol) +
      geom_text(data = tg_e, aes(x, y, label = g), size = 2.35, colour = "white", fontface = "bold", family = "Arial") +
      geom_text(data = tg_e, aes(x = 0.15, y = y, label = lab), hjust = 1, size = 2.75, fontface = "italic", colour = "#111111", family = "Arial")
  }

  p + coord_cartesian(xlim = c(x_min, length(cols) + 0.6), ylim = c(0.42, n_sp + 1.75), clip = "off", expand = FALSE) +
    th +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(), plot.margin = margin(2, 4, 1, 2))
}

pd_p2 <- make_motif_panel(
  aln_p2, 2712:2735, c(2718, 2729),
  "P2: pdm2 / nub",
  list(c("Bdor", "Bcor"), c("Zcuc", "Ztau")),
  show_tree = TRUE
)
pd_p3_tin <- make_motif_panel(
  aln_p3, 550:573, c(556, 567),
  "P3: tin",
  list(c("Bdor", "Bcor")),
  show_tree = FALSE
)
pd_p3_hlh <- make_motif_panel(
  aln_p3, 633:658, c(639, 652),
  "P3: HLH4C / Fer2 / wor",
  list(c("Bdor", "Bcor"), c("Zcuc", "Ztau")),
  show_tree = FALSE
)

pd_leg <- ggplot() +
  annotate("rect", xmin = 0.12, xmax = 0.15, ymin = 0.20, ymax = 0.80, fill = "#dceaf7", colour = "grey65", linewidth = 0.25) +
  annotate("text", x = 0.16, y = 0.50, label = "Identical", hjust = 0, size = 2.75, colour = "#111111", family = "Arial") +
  annotate("rect", xmin = 0.34, xmax = 0.37, ymin = 0.20, ymax = 0.80, fill = "#fddbc7", colour = "grey65", linewidth = 0.25) +
  annotate("text", x = 0.38, y = 0.50, label = "Substitution", hjust = 0, size = 2.75, colour = "#b2182b", fontface = "bold", family = "Arial") +
  annotate("rect", xmin = 0.58, xmax = 0.62, ymin = 0.20, ymax = 0.80, fill = NA, colour = "#08519c", linewidth = 0.80) +
  annotate("text", x = 0.63, y = 0.50, label = "Intact predicted TFBS", hjust = 0, size = 2.75, colour = "#08519c", fontface = "bold", family = "Arial") +
  coord_cartesian(xlim = c(0, 1.0), ylim = c(0, 1), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(2, 2, 1, 2))

panelD <- ((pd_p2 | pd_p3_tin | pd_p3_hlh) + plot_layout(widths = c(37.2, 24.3, 26.3))) /
          pd_leg + plot_layout(heights = c(1.0, 0.17))

# ------------------------------------------------------------------------------
# 5. Assemble clean 4-panel Figure 5 (larger 8.6 x 10.4 in canvas -> crisp 8-9.5 pt!)
# ------------------------------------------------------------------------------
row1 <- wrap_elements(panel_ab("Bmin")) + wrap_elements(panel_ab("Bdor"))
row2 <- wrap_elements(full = panelC)
row3 <- wrap_elements(full = panelD)

fig <- (row1 / row2 / row3) +
  plot_layout(heights = c(4.65, 3.05, 1.48)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 12.5, face = "bold", family = "Arial"))

PDF_OUT <- file.path(ODIR, "Figure5.pdf")
PNG_OUT <- file.path(ODIR, "Figure5.png")

cairo_pdf(PDF_OUT, width = 8.6, height = 10.4, family = "Arial")
print(fig)
dev.off()

png(PNG_OUT, width = 8.6, height = 10.4, units = "in", res = 300, type = "cairo", family = "Arial")
print(fig)
dev.off()

cat("Wrote streamlined Figure 5 ->", PDF_OUT, "\n")
