# Vector (PDF) copies of every figure used in the revised manuscript and the response letter.
# Single-panel figures are copied from the notebook outputs; composites are assembled from the panels.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(cowplot); library(magick) })
source("R/00_common.R")
OUTPDF <- file.path(dirname(ROOT), "paper_revised", "figures_pdf")
dir.create(OUTPDF, recursive = TRUE, showWarnings = FALSE)
F <- file.path(ROOT, "figures")

copy <- function(src, dest) file.copy(file.path(F, src), file.path(OUTPDF, dest), overwrite = TRUE)
pdf_of <- function(p, f, w, h) ggsave(file.path(OUTPDF, f), p, width = w, height = h, device = cairo_pdf)
panel <- function(f) ggdraw() + draw_image(file.path(F, f))

copy("Fig1a_species_tree.pdf", "Fig1a_species_tree.pdf")
copy("Fig1c.pdf", "Fig1c.pdf")
copy("Fig1e.pdf", "Fig1e.pdf")
copy("Fig1f_upset_up.pdf", "Fig1f_upregulated_DEGs.pdf")
copy("Fig2d_Gr28b_heatmap.pdf", "Fig2d_Gr28b_heatmap.pdf")
copy("Fig4.pdf", "Fig4_Gr28b_locus.pdf")
copy("SuppFig1_JSD.pdf", "SuppFig1_JSD.pdf")
copy("SuppFig2.pdf", "SuppFig2_PCA_per_tissue.pdf")
copy("SuppFig4_upset_down.pdf", "SuppFig4_downregulated_DEGs.pdf")
copy("SuppFig5_PCA_maturation_effects.pdf", "SuppFig5_PCA_maturation_effects.pdf")
copy("SuppFig6a_upset_up_all.pdf", "SuppFig6a_upregulated_DEGs.pdf")
copy("SuppFig6b_upset_down_all.pdf", "SuppFig6b_downregulated_DEGs.pdf")
copy("SuppFig7_enrichment.pdf", "SuppFig7_enrichment.pdf")
copy("SuppFig8_coupregulated_heatmap.pdf", "SuppFig8_coupregulated_genes.pdf")
copy("SuppFig13_Gr_tree.pdf", "SuppFig13_Gr_tree.pdf")
copy("SuppFig15_PCA_all.pdf", "SuppFig15a_PCA_all_genes.pdf")
copy("SuppFig15_PCA_smDEGs.pdf", "SuppFig15b_PCA_smDEGs.pdf")
copy("SuppFig15c_PCA_per_tissue.pdf", "SuppFig15c_PCA_per_tissue.pdf")
copy("SuppFig16a_clustering.pdf", "SuppFig16a_clustering.pdf")
copy("SuppFig16b_divergence.pdf", "SuppFig16b_divergence.pdf")
copy("SuppFig18_Gr_subfamily_branch_lengths.pdf", "SuppFig18_Gr_subfamily_branch_lengths.pdf")
copy("SuppFig19_Gr28b_isoform_bars.pdf", "SuppFig19_Gr28b_isoforms.pdf")

# Supp. Fig. 3: six per-tissue JSD heatmaps
pdf_of(wrap_plots(lapply(sprintf("SuppFig3_JSD_%s.pdf", c("ovipositor","foreleg","mouthpart","antenna","eye","brain")), panel), ncol = 3),
       "SuppFig3_JSD_per_tissue.pdf", 9.6, 6.4)

# Fig. 2: panels a-c are ggplot objects saved as PNG only; rebuild the same objects as PDF
suppressPackageStartupMessages({ library(Seurat); library(pheatmap); library(dplyr) })
s <- snrna_load(); ord <- SN_ORD
th <- theme_classic(base_size = 6) + theme(axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25))
num <- setNames(seq_along(ord), ord); n <- table(s$final_celltype)
s$lab <- factor(paste0(num[as.character(s$final_celltype)], "  ", as.character(s$final_celltype), " (", n[as.character(s$final_celltype)], ")"),
                levels = paste0(num, "  ", ord, " (", n[ord], ")"))
um <- as.data.frame(Embeddings(s, "umap")); um$ct <- s$final_celltype
cen <- um %>% group_by(ct) %>% summarise(x = median(umap_1), y = median(umap_2), n = dplyr::n()) %>%
  mutate(k = num[as.character(ct)], x = ifelse(n < 600, x + 1.3, x), y = ifelse(n < 600, y + 0.9, y))
pa <- DimPlot(s, group.by = "lab", pt.size = 0.05, raster = FALSE) + th + ggtitle(NULL) +
  geom_label(data = cen, aes(x, y, label = k), inherit.aes = FALSE, size = 1.9, label.size = 0.12,
             label.padding = unit(0.06, "lines"), fill = "white", alpha = 0.85) +
  guides(color = guide_legend(override.aes = list(size = 1.6), ncol = 1)) +
  theme(legend.text = element_text(size = 5.5), legend.key.size = unit(0.30, "cm"))
pdf_of(pa, "Fig2a_UMAP.pdf", 7.2, 3.8)
pdf_of(panel("Fig2b.png"), "Fig2b_focal_genes.pdf", 3.0, 1.9)
pdf_of(panel("Fig2c.png"), "Fig2c_Gr28b_NKCC.pdf", 4.2, 3.4)
pdf_of(panel("SuppFig9.png"), "SuppFig9_markers.pdf", 5.77, 7.4)
pdf_of(panel("SuppFig10.png"), "SuppFig10_focal_genes.pdf", 5.77, 4.6)
pdf_of(panel("SuppFig17_decontX_Mhc_elav.png"), "SuppFig17_ambient_RNA.pdf", 7.2, 3.8)
cat("wrote", length(list.files(OUTPDF, "\\.pdf$")), "PDFs to", OUTPDF, "\n")
