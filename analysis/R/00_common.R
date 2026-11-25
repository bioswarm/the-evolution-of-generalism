# Shared helpers sourced by every analysis script: configuration, paths, sample sheet, ortholog mapping, plotting theme.
suppressPackageStartupMessages({ library(yaml); library(data.table) })

.this_dir <- function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
}
.find_root <- function() {
  env <- Sys.getenv("REPRO_ROOT"); if (nzchar(env)) return(normalizePath(env))
  d <- .this_dir()
  for (i in 1:5) { if (file.exists(file.path(d, "config.yaml"))) return(normalizePath(d)); d <- dirname(d) }
  stop("config.yaml not found: run from the reproduce/ directory or set REPRO_ROOT")
}
ROOT <- .find_root()
CFG  <- yaml::read_yaml(file.path(ROOT, "config.yaml"))
OUT  <- normalizePath(if (nzchar(Sys.getenv("REPRO_OUT"))) Sys.getenv("REPRO_OUT") else file.path(ROOT, CFG$out), mustWork = FALSE)
RES  <- function(...) { p <- file.path(OUT, "results", ...); dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE); p }
FIG  <- function(...) { p <- file.path(OUT, "figures", ...); dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE); p }
set.seed(2026)

SPECIES <- unlist(CFG$species); AB <- unlist(CFG$abbrev); LATIN <- unlist(CFG$latin)
TISSUE_CODE <- c(A = "antenna", B = "brain", E = "eye", L = "foreleg", M = "mouthpart", O = "ovipositor")
SP_COL <- c(Bc = "#e8862a", Bd = "#c62828", Bm = "#1b7a3a", Zc = "#b84fa0", Zt = "#8c5fb8", Bt = "#2a78d6")

# ---- sample sheet: one row per kallisto run, with the corrected identity ----
sample_sheet <- function() {
  relabel <- unlist(CFG$relabel)
  rows <- lapply(SPECIES, function(sp) {
    dirs <- list.dirs(file.path(CFG$kallisto_gene, sp), recursive = FALSE, full.names = FALSE)
    dirs <- dirs[file.exists(file.path(CFG$kallisto_gene, sp, dirs, "abundance.h5"))]
    true <- ifelse(dirs %in% names(relabel), relabel[dirs], dirs)
    m <- regmatches(true, regexec("^(\\w\\w)-(?:(FN|FM|MM)-)?([ABELMO])(\\d)$", true, perl = TRUE))
    data.table(raw_dir = dirs, id = true, species = sp, sp = sapply(m, `[`, 2),
               group = ifelse(sapply(m, `[`, 3) == "", "LA", sapply(m, `[`, 3)),
               tcode = sapply(m, `[`, 4), rep = as.integer(sapply(m, `[`, 5)),
               path = file.path(CFG$kallisto_gene, sp, dirs))
  })
  s <- rbindlist(rows)
  s[, tissue := ifelse(group == "LA", paste0("larval_", TISSUE_CODE[tcode]), TISSUE_CODE[tcode])]
  s[, sample := gsub("-", ".", id)][, relabelled := raw_dir != id][]
}

# kallisto target id -> gene id used in the ortholog table
to_gene <- function(sp, tx) if (sp == "Bcor") sub("^Bc_(\\d+)-P\\d+$", "Bc\\1", tx) else sub("\\.\\d+$", "", tx)
orthologues <- function() fread(CFG$orthologues)
fly_names <- function() { x <- fread(CFG$fly_names, header = FALSE); setNames(x$V2, x$V1) }

# ---- ortholog-level TPM matrix (rows = Bdor gene id of each orthogroup; columns = samples) ----
ortho_tpm <- function(ss = sample_sheet(), renormalise = TRUE) {
  f <- RES("tables", "ortho_tpm.tsv")
  if (file.exists(f)) { m <- as.matrix(fread(f), rownames = 1); return(m[, ss$sample, drop = FALSE]) }
  o <- orthologues()
  m <- sapply(seq_len(nrow(ss)), function(i) {
    a <- fread(file.path(ss$path[i], "abundance.tsv"), select = c("target_id", "tpm"))
    g <- to_gene(ss$species[i], a$target_id); a$tpm[match(o[[ss$species[i]]], g)]
  })
  dimnames(m) <- list(o$Bdor, ss$sample)
  if (renormalise) m <- sweep(m, 2, colSums(m, na.rm = TRUE), "/") * 1e6
  fwrite(data.table(gene = rownames(m), m), f, sep = "\t"); m
}

theme_pub <- function(base = 8) ggplot2::theme_classic(base_size = base, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3), axis.ticks = ggplot2::element_line(linewidth = 0.3),
                 axis.text = ggplot2::element_text(size = base * 0.88, colour = "grey15"),
                 strip.background = ggplot2::element_blank(),
                 strip.text = ggplot2::element_text(size = base, face = "bold"),
                 legend.text = ggplot2::element_text(size = base * 0.88))

save_arial_pdf <- function(path, width, height, plot_obj = NULL, expr = NULL) {
  cairo_pdf(path, width = width, height = height, family = "Arial")
  on.exit(dev.off(), add = TRUE)
  if (!is.null(plot_obj)) print(plot_obj) else eval(expr)
}
# ---- Fig. 1c / 1e style: open symbols, species in italics (tau..correcta), bold group labels with leader lines ----
FIG1_SHAPE <- c(Antenna = 1, Brain = 4, Eye = 0, Foreleg = 6, Mouthpart = 5, Ovipositor = 2, `Larval brain` = 3, `Larval mouthpart` = 4)
# PC signs are arbitrary: orient so that the named tissue lies on the given side (as in the published panels)
orient_pc <- function(d, pc, tissue, sign) { if (sign(mean(d[[pc]][d$tissue == tissue])) != sign) d[[pc]] <- -d[[pc]]; d }
FIG1_SP <- c(Zt = "tau", Zc = "cucurbitae", Bm = "minax", Bd = "dorsalis", Bc = "correcta")
fig1_scatter <- function(d, ve, labels, size = 1.4, tissues = names(FIG1_SHAPE)) {
  # d: data.table with PC1, PC2, sp (Bc/Bd/...), tissue (names of FIG1_SHAPE); labels: data.table(label, x, y, gx, gy)
  #    (x, y = label position; gx, gy = point the leader line points to, NA for no line)
  d <- copy(d)[, st := ifelse(grepl("^Larval", tissue), 1.3, 0.45)][, sp := factor(sp, levels = names(FIG1_SP))][, tissue := factor(tissue, levels = tissues)]
  ggplot2::ggplot(d, ggplot2::aes(PC1, PC2, colour = sp, shape = tissue)) +
    # leader line from just outside the label (12% of the way) to 92% of the way to the group centre
    ggplot2::geom_segment(data = labels[!is.na(gx)], ggplot2::aes(x = x + 0.12 * (gx - x), y = y + 0.12 * (gy - y),
                          xend = x + 0.92 * (gx - x), yend = y + 0.92 * (gy - y)), inherit.aes = FALSE, colour = "grey45", linewidth = 0.25) +
    ggplot2::geom_point(ggplot2::aes(stroke = st), size = size) +
    ggplot2::geom_text(data = labels, ggplot2::aes(x, y, label = label), inherit.aes = FALSE, size = 2.2, fontface = "bold", colour = "grey20",
                       lineheight = 0.9) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::scale_colour_manual(values = SP_COL[names(FIG1_SP)], labels = FIG1_SP, name = NULL) +
    ggplot2::scale_shape_manual(values = FIG1_SHAPE[tissues], drop = FALSE, name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 1, override.aes = list(shape = 16, size = 2.2)),
                    shape = ggplot2::guide_legend(order = 2, override.aes = list(colour = "grey35", size = 1.8,
                                                  stroke = ifelse(grepl("^Larval", tissues), 1.3, 0.45)))) +
    ggplot2::labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) + theme_pub(7) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 6), legend.key.size = ggplot2::unit(0.32, "cm"),
                   legend.spacing.y = ggplot2::unit(0, "mm"), legend.margin = ggplot2::margin(0, 0, 0, 0)) }
fig1_italic_species <- function(p) p + ggplot2::guides(colour = ggplot2::guide_legend(order = 1, override.aes = list(shape = 16, size = 2.2),
                                                       theme = ggplot2::theme(legend.text = ggplot2::element_text(face = "italic", size = 6))))
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# ---- notebook helpers ----
# In Jupyter, each cell runs on a ragg recording device opened by the evaluate package. Code that measures text before
# opening its own device (ComplexHeatmap, grid) would then use ragg font metrics instead of the pdf metrics of a script
# run (default device pdf("Rplots.pdf")). Use evaluate's pdf(NULL) fallback so notebook and script runs give identical files.
if (isTRUE(getOption("jupyter.in_kernel"))) {
  utils::assignInNamespace("has_ragg", function() FALSE, "evaluate")
  # ComplexHeatmap::draw, when it sees IRkernel loaded, reroutes the plot through grid.grabExpr + replay, which
  # rescales it slightly (Supp. Fig. 8 heatmap). We save and show figures through our own device / show_fig(), so
  # force the plain, script-identical path. Patch now if loaded, otherwise when the package loads.
  .no_jup <- function() try(utils::assignInNamespace("is_under_jupyter", function() FALSE, "ComplexHeatmap"), silent = TRUE)
  if ("ComplexHeatmap" %in% loadedNamespaces()) .no_jup() else setHook(packageEvent("ComplexHeatmap", "onLoad"), function(...) .no_jup())
}
# show a saved figure file inline (PDF rasterised at 100 dpi, PNG scaled to <= width_px); no-op outside Jupyter
show_fig <- function(path, width_px = 1100) {
  if (!isTRUE(getOption("jupyter.in_kernel"))) return(invisible(path))
  im <- if (grepl("\\.pdf$", path)) magick::image_read(path, density = 100) else magick::image_read(path)
  for (k in seq_along(im)) { p <- magick::image_background(im[k], "white")
    if (magick::image_info(p)$width > width_px) p <- magick::image_resize(p, sprintf("%dx", width_px))
    IRdisplay::display_png(magick::image_write(p, format = "png")) }
  invisible(path) }
# stop early when an upstream notebook has not been run
need <- function(..., by) { f <- c(...); miss <- f[!file.exists(f)]
  if (length(miss)) stop("missing ", paste(miss, collapse = ", "), "\nrun notebooks/", by, ".ipynb first", call. = FALSE)
  invisible(TRUE) }
# compare a value with the number quoted in the manuscript
check_num <- function(label, value, expected, tol = 0) {
  ok <- isTRUE(abs(value - expected) <= tol + 1e-12)
  cat(sprintf("%-60s %s (text: %s) %s\n", label, format(value), format(expected), if (ok) "OK" else "MISMATCH"))
  if (!ok) stop("value differs from the manuscript: ", label, call. = FALSE)
  invisible(ok) }

# ---- helpers shared by several figure notebooks (ggplot2 / patchwork must be loaded by the caller) ----
TIS_LAB <- c(antenna = "Antenna", brain = "Brain", eye = "Eye", foreleg = "Foreleg", mouthpart = "Mouthpart", ovipositor = "Ovipositor",
      larval_brain = "Larval brain", larval_mouthpart = "Larval mouthpart")
SP_LAB <- c(Bc = "B. correcta", Bd = "B. dorsalis", Bm = "B. minax", Zc = "B. cucurbitae", Zt = "B. tau", Bt = "B. tsuneonis")   # genus as in the text (R2.7)
SHAPE <- c(Antenna = 1, Brain = 4, Eye = 0, Foreleg = 6, Mouthpart = 5, Ovipositor = 2, `Larval brain` = 3, `Larval mouthpart` = 8)
pca <- function(M, scale = TRUE, n = 20) { M <- M[apply(M, 1, var) > 0, , drop = FALSE]
  p <- prcomp(t(M), center = TRUE, scale. = scale, rank. = n); list(x = p$x, ve = 100 * p$sdev^2 / sum(p$sdev^2)) }
scatter <- function(p, meta, title = NULL, flip2 = FALSE, base = 8) {
  d <- data.table(PC1 = p$x[, 1], PC2 = p$x[, 2] * ifelse(flip2, -1, 1), meta)
  ggplot(d, aes(PC1, PC2, colour = sp, shape = TIS_LAB[tissue])) + geom_point(size = 1.4, stroke = 0.55) +
    scale_colour_manual(values = SP_COL, labels = SP_LAB, name = NULL) + scale_shape_manual(values = SHAPE, name = NULL) +
    labs(x = sprintf("PC1 (%.1f%%)", p$ve[1]), y = sprintf("PC2 (%.1f%%)", p$ve[2]), title = title) + theme_pub(base) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 2.0)),
           shape  = guide_legend(order = 2, override.aes = list(size = 2.0))) +
    theme(plot.title = element_text(size = base + 0.5, hjust = 0.5, face = "bold"),
          legend.text = element_text(size = base - 0.5),
          legend.key.size = unit(0.32, "cm"))
}
# published Fig. 1c input (TPM of 9,801 1:1 orthologs, 285 libraries); columns already carry corrected Bc-MM-L/M labels
fig1c_tpm <- function(ss) {
  tp <- as.matrix(fread(CFG$fig1c_tpm), rownames = 1); tp <- tp[, !grepl("^Bt", colnames(tp))]
  colnames(tp) <- gsub("-", ".", sub("-LA-?", "-", colnames(tp)))     # Bc-LA-B1 -> Bc.B1 (sample sheet naming)
  stopifnot(setequal(colnames(tp), ss$sample)); tp[, ss$sample] }

# published per-tissue mashr pattern codes (5 characters: tau, cucurbitae, minax, dorsalis, correcta; 2 up, 1 ns, 0 down)
PT <- c(antenna = "Ant", brain = "Br", eye = "Eye", foreleg = "Leg", mouthpart = "Mo", ovipositor = "Ovi")
published_patterns <- function() lapply(PT, function(d) {
  x <- fread(file.path(CFG$published_patterns, d, "DEG.patterns"), header = FALSE, skip = 1, colClasses = "character"); setNames(x$V2, x$V1) })
# UpSet-style intersection panels (Fig. 1f up; Supp. Fig. 6 both directions, all intersections)
upset_panel <- function(direction, cols = NULL, pub = published_patterns()) {
  ch <- if (direction == "up") "2" else "0"
  inter <- lapply(pub, function(p) { b <- vapply(strsplit(p, ""), function(v) paste(as.integer(v == ch), collapse = ""), ""); table(b[b != "00000"]) })
  if (is.null(cols)) { tot <- Reduce(function(a, b) { k <- union(names(a), names(b)); setNames(ifelse(is.na(a[k]), 0, a[k]) + ifelse(is.na(b[k]), 0, b[k]), k) }, inter)
    cols <- names(sort(tot, decreasing = TRUE)) }
  d <- rbindlist(lapply(names(inter), function(t) data.table(tissue = t, pattern = cols, n = as.integer(ifelse(is.na(inter[[t]][cols]), 0, inter[[t]][cols])))))
  d[, pattern := factor(pattern, levels = cols)][, tissue := factor(tissue, levels = c("brain", "eye", "antenna", "mouthpart", "foreleg", "ovipositor"))]
  bars <- ggplot(d, aes(pattern, n)) + geom_col(fill = "#1f9fe0", width = 0.8) + facet_grid(tissue ~ ., switch = "y") +
    geom_hline(yintercept = 50, linetype = 2, linewidth = 0.3) + labs(x = NULL, y = sprintf("#%sregulated genes after sexual maturity", if (direction == "up") "Up" else "Down")) +
    theme_pub() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.placement = "outside", strip.text.y.left = element_text(angle = 0))
  dm <- CJ(pattern = factor(cols, levels = cols), k = 1:5)[, on := substr(as.character(pattern), k, k) == "1"]
  dm[, row := factor(c("tau", "cucurbitae", "minax", "dorsalis", "correcta")[k], levels = rev(c("tau", "cucurbitae", "minax", "dorsalis", "correcta")))]
  seg <- dm[on == TRUE, .(lo = min(as.integer(row)), hi = max(as.integer(row))), by = pattern]
  dots <- ggplot(dm, aes(pattern, row)) + geom_point(aes(colour = ifelse(on, as.character(row), "off")), size = 1.3) +
    geom_segment(data = seg, aes(x = pattern, xend = pattern, y = lo, yend = hi), linewidth = 0.4, inherit.aes = FALSE) +
    scale_colour_manual(values = c(setNames(SP_COL[c("Zt", "Zc", "Bm", "Bd", "Bc")], c("tau", "cucurbitae", "minax", "dorsalis", "correcta")), off = "grey88"), guide = "none") +
    labs(x = NULL, y = NULL) + theme_pub() + theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
    axis.text.y = element_text(face = "italic"))
  list(plot = bars / dots + plot_layout(heights = c(6, 1.2)), cols = cols, data = d)
}

# re-annotated snRNA object (notebook 00c) with the final labels as ordered identities (Seurat must be loaded by the caller)
SN_ORD <- c("Muscle", "Fat body", "Hemocytes", "Immune-activated hemocytes (AMP+)", "Tracheal epithelium",
     "Cuticular epithelium", "Epithelium (aspr+ Eaat1+)", "Epithelium (Ance-3+ Ddc+)", "Epithelium (glob1+ ITP+)",
     "Epithelium (FucTA+ KCNQ+)", "Epithelium (cad+ Poxm+)", "Ion-transporting epithelium (CAH1+ Nha1+)",
     "Gr28b+ epithelium (Mdr49+ NaPi-T+)", "Sensory organ support cells (sv+ Su(H)+)", "Sensory organ support cells (sv+ dsb+)",
            "Glia (repo+)", "Mechanosensory neurons (nSyb+ ChAT+ nompC+)")
snrna_load <- function() {
  s <- readRDS(RES("snRNA", "BdOV_decontX_reannot.rds")); DefaultAssay(s) <- "SCT"
  stopifnot(setequal(SN_ORD, unique(as.character(s$final_celltype))))
  s$final_celltype <- factor(as.character(s$final_celltype), levels = SN_ORD); Idents(s) <- "final_celltype"; s }

# ---- composite figures (300 dpi, text width 5.77 in, bold 9-pt panel letters); magick must be loaded by the caller ----
DPI <- 300; TW <- 5.77; PAD <- round(0.06 * DPI); LSIZE <- round(9 / 72 * DPI)
FINAL <- function(f) FIG("final", f)
load_panel <- function(f, width_in = NULL) {
  im <- if (grepl("\\.pdf$", f)) image_read(f, density = DPI) else image_read(f)
  im <- image_background(image_flatten(im), "white")
  if (!is.null(width_in)) im <- image_resize(im, sprintf("%dx", round(width_in * DPI)), filter = "Lanczos")
  im }
# rows: list of rows, each a list of list(file, width_in); letters assigned left -> right, top -> bottom
stack_panels <- function(rows, letters, width_in = TW) {
  W <- round(width_in * DPI); k <- 0; row_ims <- list()
  for (r in rows) {
    ims <- lapply(r, function(p) load_panel(p[[1]], p[[2]])); h <- max(vapply(ims, function(i) image_info(i)$height, 0))
    row <- image_blank(W, h + PAD, "white"); x <- PAD
    for (im in ims) {
      k <- k + 1; row <- image_composite(row, im, offset = sprintf("+%d+%d", x, PAD))
      if (nzchar(letters[k])) row <- image_annotate(row, letters[k], font = "Arial", weight = 700, size = LSIZE, location = sprintf("+%d+%d", x - PAD, 0))
      x <- x + image_info(im)$width + PAD }
    row_ims[[length(row_ims) + 1]] <- row }
  out <- image_append(do.call(c, row_ims), stack = TRUE)
  image_append(c(out, image_blank(W, PAD, "white")), stack = TRUE) }
save_final <- function(im, name) { dir.create(dirname(FINAL("x")), showWarnings = FALSE, recursive = TRUE)
  image_write(im, FINAL(name), format = "png", density = DPI)
  i <- image_info(im); log_msg(sprintf("%-14s %d x %d px (%.2f x %.2f in)", name, i$width, i$height, i$width / DPI, i$height / DPI)) }
