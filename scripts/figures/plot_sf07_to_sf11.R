#!/usr/bin/env Rscript
# Generate Supplementary Figs. 7a, 7b, 8, 9, 10, and 11 in 100% TrueType Arial:
#   - Supplementary Fig. 7a/7b: Tissue UpSet plots of upregulated / downregulated DEGs
#   - Supplementary Fig. 8: GO and KEGG enrichment of ovipositor maturation-associated genes
#   - Supplementary Fig. 9: Ambient RNA decontamination (decontX) across ovipositor cell populations
#   - Supplementary Fig. 10: Canonical marker gene dot plot across all 17 snRNA-seq cell populations
#   - Supplementary Fig. 11: Single-nucleus UMAP feature plots of the 11 detected ovipositor candidate genes

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
  library(grid)
  library(ComplexHeatmap)
  library(dplyr)
  library(tidyr)
  library(Seurat)
})

save_both <- function(pdf_path, png_path, w, h, obj = NULL, expr = NULL) {
  cairo_pdf(pdf_path, width = w, height = h, family = "Arial")
  if (!is.null(obj)) print(obj) else eval(expr)
  dev.off()
  png(png_path, width = w, height = h, units = "in", res = 300, type = "cairo", family = "Arial")
  if (!is.null(obj)) print(obj) else eval(expr)
  dev.off()
  cat(sprintf("Wrote %s\n", pdf_path))
}
copy_renum <- function(...) invisible(NULL)
ROOT_REV <- REVIEW_ROOT
# ==============================================================================
# 5. SuppFig 6a & 6b: UpSet plots of up/downregulated DEGs -> Arial 7.5-10.5 pt
# ==============================================================================
log_msg("--- SuppFig 6a & 6b ---")
pub <- published_patterns()
tis_up_down <- function(code) {
  sets <- lapply(pub, function(p) names(p)[p == code])
  names(sets) <- TIS_LAB[names(sets)]
  sets[c("Ovipositor", "Eye", "Foreleg", "Antenna", "Mouthpart", "Brain")]
}
upset_tissue <- function(code, colour, title) {
  sets <- tis_up_down(code)
  m <- make_comb_mat(sets)
  m <- m[comb_size(m) > 0]
  o <- order(-comb_size(m))
  UpSet(m, comb_order = o, set_order = rev(names(sets)), comb_col = colour, bg_col = "#F0F0F0",
        row_names_gp = gpar(fontsize = 8.5, fontfamily = "Arial"),
        top_annotation = upset_top_annotation(
          m, add_numbers = TRUE, numbers_gp = gpar(fontsize = 7.0, fontfamily = "Arial"),
          gp = gpar(fill = colour, col = NA), height = unit(4.2, "cm"),
          annotation_name_rot = 90, annotation_name_gp = gpar(fontsize = 8.5, fontfamily = "Arial"),
          axis_param = list(gp = gpar(fontsize = 8, fontfamily = "Arial"))
        ),
        left_annotation = upset_left_annotation(
          m, gp = gpar(fill = colour, col = NA),
          annotation_name_gp = gpar(fontsize = 8.5, fontfamily = "Arial"),
          axis_param = list(direction = "reverse", gp = gpar(fontsize = 8, fontfamily = "Arial"))
        ),
        right_annotation = NULL, row_names_side = "right",
        column_title = title, column_title_gp = gpar(fontsize = 10.5, fontface = "bold", fontfamily = "Arial"))
}
ht6a <- upset_tissue("22222", "#E41A1C", "Upregulated DEGs")
ht6b <- upset_tissue("00000", "#3778B4", "Downregulated DEGs")
save_both(file.path(OUT_DIR, "Supplementary_Figure_7a_upregulated_DEGs.pdf"), file.path(OUT_DIR, "Supplementary_Figure_7a_upregulated_DEGs.png"), 7.5, 4.6, expr = quote(draw(ht6a)))
save_both(file.path(OUT_DIR, "Supplementary_Figure_7b_downregulated_DEGs.pdf"), file.path(OUT_DIR, "Supplementary_Figure_7b_downregulated_DEGs.png"), 7.5, 4.6, expr = quote(draw(ht6b)))
copy_renum(file.path(OUT_DIR, "Supplementary_Figure_7a_upregulated_DEGs.pdf"), "SuppFig6a_upregulated_DEGs.pdf")
copy_renum(file.path(OUT_DIR, "Supplementary_Figure_7b_downregulated_DEGs.pdf"), "SuppFig6b_downregulated_DEGs.pdf")

# ==============================================================================
# 6. SuppFig 7: GO & KEGG enrichment -> Arial 8.5 pt
# ==============================================================================
log_msg("--- SuppFig 7 ---")
res7 <- fread(RES("tables", "SuppData_enrichment_ovipositor.tsv"))
bars7 <- function(d, fill_lo, fill_hi, ylab) {
  d <- d[order(qvalue)][1:min(.N, 16)]
  d[, Description := factor(Description, levels = rev(Description))]
  ggplot(d, aes(Count, Description, fill = qvalue)) + geom_col(width = 0.78) +
    scale_fill_gradient(low = fill_lo, high = fill_hi, name = expression(italic(q))) +
    labs(x = "#Genes", y = ylab) + theme_pub(8.5) +
    theme(axis.title.y = element_text(size = 9, face = "bold"),
          axis.text.y = element_text(size = 7.8),
          legend.title = element_text(size = 8.5), legend.text = element_text(size = 7.5))
}
sig7 <- res7[qvalue < 0.05]
p_sf7 <- (bars7(sig7[grepl("^GO", db)], "#D6604D", "#4393C3", "Gene ontology terms") /
          bars7(sig7[db == "KEGG"], "#D6604D", "#4393C3", "KEGG pathways")) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 11, face = "bold", family = "Arial"))
save_both(file.path(OUT_DIR, "Supplementary_Figure_8_ovipositor_enrichment.pdf"), file.path(OUT_DIR, "Supplementary_Figure_8_ovipositor_enrichment.png"), 7.5, 8.2, obj = p_sf7)
copy_renum(file.path(OUT_DIR, "Supplementary_Figure_8_ovipositor_enrichment.pdf"), "SuppFig7_enrichment.pdf")

# ==============================================================================
# 7. SuppFig 8 (markers DotPlot), SuppFig 9 (focal genes FeaturePlot), SuppFig 15 (ambient RNA)
# ==============================================================================
log_msg("--- SuppFig 8, 9, 15 (snRNA) ---")
sn_rds <- Sys.getenv("SNRNA_RDS", file.path(ROOT_REV, "improved_20261113/reproduce/results/snRNA/BdOV_decontX_reannot.rds"))
s <- readRDS(sn_rds); DefaultAssay(s) <- "SCT"
s$final_celltype <- factor(as.character(s$final_celltype), levels = SN_ORD)
Idents(s) <- "final_celltype"
num <- setNames(seq_along(SN_ORD), SN_ORD)

# SuppFig 8: Marker DotPlot (true vector Arial 7.5-8.5 pt!)
mk <- c("Mhc", "sls", "Zasp66", "TpnC73F", "Lsp2", "Fbp1", "apolpp", "Hml", "Pxn", "Def", "AttA", "DptB", "trh", "btl",
        "kkv", "grh", "Cht7", "Cpr49Ae", "aspr", "Eaat1", "Ance-3", "Ddc", "glob1", "ITP", "FucTA", "KCNQ", "cad", "Poxm",
        "CAH1", "Nha1", "Mdr49", "NaPi-T", "sv", "Su(H)", "dsb", "repo", "nrv2", "nSyb", "brp", "Syt1", "ChAT", "nompC")
mk <- intersect(mk, rownames(s))
DISP_SN <- ifelse(SN_ORD == "Gr28b+ epithelium (Mdr49+ NaPi-T+)",
                  "Epithelium (Gr28b+ Mdr49+ NaPi-T+)", SN_ORD)
s$ct_num <- factor(paste0(num[as.character(s$final_celltype)], "  ", DISP_SN[num[as.character(s$final_celltype)]]),
                   levels = paste0(num, "  ", DISP_SN))
pd8 <- DotPlot(s, features = rev(mk), group.by = "ct_num", dot.scale = 3.2) + coord_flip() + theme_pub(8.5) +
  theme(axis.text.y = element_text(size = 7.6, face = "italic"),
        axis.text.x = element_text(size = 7.4, angle = 45, hjust = 1),
        legend.text = element_text(size = 7.5), legend.title = element_text(size = 8.0),
        legend.key.size = unit(0.36, "cm"), plot.margin = margin(4, 6, 4, 4)) +
  xlab(NULL) + ylab(NULL)
save_both(file.path(OUT_DIR, "Supplementary_Figure_10_snRNA_markers.pdf"), file.path(OUT_DIR, "Supplementary_Figure_10_snRNA_markers.png"), 7.5, 9.2, obj = pd8)
copy_renum(file.path(OUT_DIR, "Supplementary_Figure_10_snRNA_markers.pdf"), "SuppFig8_markers.pdf")

# SuppFig 9: Focal genes UMAP FeaturePlots — ONE unified color scale across all 11 panels!
focal <- c("Pde8", "Bd4G00733", "yellow-g", "Lmpt", "Cyp318a1", "Obp56c", "Bd4G00648", "Yp3", "Esp", "Bd1G02800", "Bd3G00354", "Gr28b", "NKCC")
alias <- c(Bd4G00733 = "CG3191", Bd1G02800 = "CG6675", Bd3G00354 = "CG34290")
feat <- sapply(focal, function(g) { q <- if (g %in% names(alias)) alias[[g]] else g
  h <- grep(paste0("^", q, "(\\.1)?$"), rownames(s), value = TRUE); if (length(h)) h[1] else NA })
feat <- feat[!is.na(feat)]
hot <- scale_colour_gradient(
  low = "grey88", high = "#7f2704",
  limits = c(0, 4), oob = scales::squish,
  breaks = 0:4, labels = c("0", "1", "2", "3", "\u22654"),
  name = "Expression\n(SCT log-norm)",
  guide = guide_colourbar(barheight = unit(2.4, "cm"), barwidth = unit(0.32, "cm"), title.position = "top")
)
ttl <- setNames(names(feat), feat)
ttl[feat %in% alias] <- paste0(names(feat)[feat %in% alias], " (", feat[feat %in% alias], ")")
fp9 <- function(g, title = g) FeaturePlot(s, features = g, order = TRUE, pt.size = 0.06, raster = FALSE) +
  hot + theme_pub(8.0) + ggtitle(title) +
  theme(plot.title = element_text(size = 8.5, face = "italic", hjust = 0.5),
        legend.title = element_text(size = 8.0, face = "bold"),
        legend.text = element_text(size = 7.5), axis.title = element_blank(),
        axis.text = element_text(size = 7.0))
p_sf9 <- (wrap_plots(c(lapply(feat, function(g) fp9(g, ttl[[g]])), list(guide_area())), ncol = 4) +
          plot_layout(guides = "collect"))
save_both(file.path(OUT_DIR, "Supplementary_Figure_11_focal_genes_atlas.pdf"), file.path(OUT_DIR, "Supplementary_Figure_11_focal_genes_atlas.png"), 7.5, 5.8, obj = p_sf9)
copy_renum(file.path(OUT_DIR, "Supplementary_Figure_11_focal_genes_atlas.pdf"), "SuppFig9_focal_genes.pdf")

# SuppFig 15: Ambient RNA removal (true vector Arial 7.5-8.5 pt!)
dd15 <- fread(RES("snRNA", "SuppFig17_decontX_detection.tsv"))
dd15[, celltype := factor(ifelse(celltype == "Gr28b+ epithelium (Mdr49+ NaPi-T+)",
                                 "Epithelium (Gr28b+ Mdr49+ NaPi-T+)", celltype), levels = DISP_SN)]
dd15_l <- melt(dd15, id.vars = c("gene", "celltype"), measure.vars = c("before", "after"),
               variable.name = "stage", value.name = "value")
pe15 <- ggplot(dd15_l, aes(value, celltype, fill = stage)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  facet_wrap(~ gene, nrow = 1) +
  scale_fill_manual(values = c(before = "#b5b4ae", after = "#2a78d6"),
                    labels = c("Before decontX", "After decontX"), name = NULL) +
  xlab("Nuclei with >= 1 UMI (%)") + ylab(NULL) + theme_pub(8.5) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 7.5),
        strip.text = element_text(size = 9, face = "bold.italic"))
save_both(file.path(OUT_DIR, "Supplementary_Figure_9_ambient_RNA_decontX.pdf"), file.path(OUT_DIR, "Supplementary_Figure_9_ambient_RNA_decontX.png"), 7.5, 4.6, obj = pe15)
copy_renum(file.path(OUT_DIR, "Supplementary_Figure_9_ambient_RNA_decontX.pdf"), "SuppFig15_ambient_RNA.pdf")
rm(s); gc(FALSE)
cat("Done generating Supplementary Figs. 7a-11.\n")
