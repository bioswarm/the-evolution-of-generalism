#!/usr/bin/env Rscript
# Pure-source-code Supplementary Fig. 12 (enlarged fonts, NO abbreviated species names):
#   a: Row-aligned N-terminal sequence identity matrix (vs B. dorsalis X1-X5 and D. melanogaster Gr28b.A-E)
#      linked seamlessly to a shared-scale multi-tissue expression heatmap (log2(TPM + 1), nF vs mF)
#      across all six Bactrocera species and D. melanogaster (FlyBase r6.69 kallisto).
#   b: Mature vs newly emerged ovipositor Gr28b isoform TPM stacked bars across species.
#   c: Relative isoform composition (%) in mature ovipositors across species.
# 100% embedded TrueType Arial via cairo_pdf.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

REPO_ROOT <- Sys.getenv("REPO_ROOT", ".")
OUT_DIR   <- Sys.getenv("OUT_DIR", file.path(REPO_ROOT, "outputs", "figures"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
bact   <- fread(file.path(REPO_ROOT, "data", "supplementary_data", "SuppData4b_Gr28b_isoform_tpm_per_sample.tsv"))
id_mat <- fread(file.path(REPO_ROOT, "data", "Gr28b_Nterm_identity_matrix.tsv"))
fly    <- fread(file.path(REPO_ROOT, "data", "Dmel_Gr28b_tpm.csv"))

# Ordered row specification (33 isoforms: 28 Bactrocera + 5 D. melanogaster)
# Full standard italicized species names ONLY (never B. min. / B. cor. etc.)
row_spec <- rbindlist(list(
  data.table(sp = "Bcor", sp_lab = "B. correcta",     iso = c("Bc_003035-R1", "Bc_003035-R2"),
             disp = c("R1 (X1)", "R2 (X4-like)"),
             cls  = c("X1 / Gr28b.C", "X4 / X4-like")),
  data.table(sp = "Bdor", sp_lab = "B. dorsalis",     iso = c("Bd1G02958.1", "Bd1G02958.2", "Bd1G02958.3", "Bd1G02958.4", "Bd1G02958.5"),
             disp = c(".1 Gr28b-X1", ".2 X2", ".3 X3", ".4 X4", ".5 X5"),
             cls  = c("X1 / Gr28b.C", "X2 / X2-like", "X3 / X3-like", "X4 / X4-like", "X5 / Gr28b.D")),
  data.table(sp = "Zcuc", sp_lab = "B. cucurbitae",   iso = c("Zc1G00533.1", "Zc1G00533.2", "Zc1G00533.3", "Zc1G00533.4", "Zc1G00533.5"),
             disp = c(".1 (X1)", ".2 (Gr28b.B-like)", ".3 (X3-like)", ".4 (X3-like, ovi.)", ".5 (X5)"),
             cls  = c("X1 / Gr28b.C", "Gr28b.B-like", "X3 / X3-like", "X3 / X3-like", "X5 / Gr28b.D")),
  data.table(sp = "Ztau", sp_lab = "B. tau",          iso = c("Zt1G00596.1", "Zt1G00596.2", "Zt1G00596.3", "Zt1G00596.4", "Zt1G00596.5", "Zt1G00596.6", "Zt1G00596.7", "Zt1G00596.10", "Zt1G00596.8"),
             disp = c(".1 (X1)", ".2 (Gr28b.B-like)", ".3 (X2-like)", ".4 (X3-like)", ".5 (X3-like, ovi.)", ".6 (X3-like)", ".7 (X3-like)", ".10 (X4-like)", ".8 (X5)"),
             cls  = c("X1 / Gr28b.C", "Gr28b.B-like", "X2 / X2-like", "X3 / X3-like", "X3 / X3-like", "X3 / X3-like", "X3 / X3-like", "X4 / X4-like", "X5 / Gr28b.D")),
  data.table(sp = "Bmin", sp_lab = "B. minax",        iso = c("Bm1G00667.3", "Bm1G00667.2", "Bm1G00667.1", "Bm1G00667.4", "Bm1G00667.5"),
             disp = c(".3 (X1)", ".2 (Gr28b.B-like)", ".1 (X3-like)", ".4 (X4-like)", ".5 (X5)"),
             cls  = c("X1 / Gr28b.C", "Gr28b.B-like", "X3 / X3-like", "X4 / X4-like", "X5 / Gr28b.D")),
  data.table(sp = "Btsu", sp_lab = "B. tsuneonis",    iso = c("Bt_011990-R1", "Bt_011990-R2"),
             disp = c("R1 (X1)", "R2 (X4-like)"),
             cls  = c("X1 / Gr28b.C", "X4 / X4-like")),
  data.table(sp = "Dmel", sp_lab = "D. melanogaster", iso = c("Gr28b.C", "Gr28b.A", "Gr28b.B", "Gr28b.E", "Gr28b.D"),
             disp = c("Gr28b.C", "Gr28b.A", "Gr28b.B", "Gr28b.E", "Gr28b.D"),
             cls  = c("X1 / Gr28b.C", "Gr28b.A / E", "Gr28b.B-like", "Gr28b.A / E", "X5 / Gr28b.D"))
))

sp_order <- c("B. correcta", "B. dorsalis", "B. cucurbitae", "B. tau",
              "B. minax", "B. tsuneonis", "D. melanogaster")
row_spec[, sp_fac := factor(sp_lab, levels = sp_order)]
row_spec[, row_id := factor(iso, levels = rev(iso))]

cls_levels <- c("X1 / Gr28b.C", "X3 / X3-like", "X4 / X4-like", "X2 / X2-like",
                "X5 / Gr28b.D", "Gr28b.B-like", "Gr28b.A / E")
cls_cols <- c(
  "X1 / Gr28b.C" = "#e6ab02",
  "X3 / X3-like" = "#d95f02",
  "X4 / X4-like" = "#e78ac3",
  "X2 / X2-like" = "#8da0cb",
  "X5 / Gr28b.D" = "#2b83ba",
  "Gr28b.B-like" = "#7570b3",
  "Gr28b.A / E"  = "#bdbdbd"
)
row_spec[, cls_fac := factor(cls, levels = cls_levels)]

# ------------------------------------------------------------------------------
# Panel a-left: N-terminal sequence identity heatmap (Bdor X1-X5 | Dmel C, A, B, E, D)
# ------------------------------------------------------------------------------
seq_cols <- c("Bdor_X1", "Bdor_X2", "Bdor_X3", "Bdor_X4", "Bdor_X5",
              "Dmel_C",  "Dmel_A",  "Dmel_B",  "Dmel_E",  "Dmel_D")
seq_xlabs <- c("X1", "X2", "X3", "X4", "X5", "C", "A", "B", "E", "D")

id_long <- melt(id_mat[isoform %in% row_spec$iso],
                id.vars = c("species", "isoform"),
                measure.vars = seq_cols,
                variable.name = "ref", value.name = "pct")
id_long <- merge(id_long, row_spec[, .(iso, sp_fac, disp, row_id, cls_fac)],
                 by.x = "isoform", by.y = "iso", sort = FALSE)
id_long[, ref := factor(ref, levels = seq_cols)]
id_long[, xnum := as.numeric(ref) + ifelse(as.numeric(ref) > 5, 0.35, 0)]
id_long[, txt := ifelse(pct >= 35, sprintf("%.0f", pct), "")]

p_cls <- ggplot(row_spec, aes(x = 1, y = row_id, fill = cls_fac)) +
  geom_tile(width = 0.85, height = 0.88, colour = "white", linewidth = 0.28) +
  scale_fill_manual(values = cls_cols, name = "Isoform class", drop = FALSE) +
  scale_y_discrete(labels = setNames(row_spec$disp, row_spec$iso), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_grid(sp_fac ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL, title = "") +
  theme_minimal(base_size = 8.5, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(size = 7.6, colour = "#111111", hjust = 1, margin = margin(r = 3)),
    strip.placement = "outside",
    strip.text.y.left = element_text(size = 8.4, face = "bold.italic", angle = 0, hjust = 1, colour = "#111111",
                                     margin = margin(r = 5)),
    strip.background = element_blank(),
    panel.spacing.y = unit(4.5, "pt"),
    legend.position = "none",
    plot.margin = margin(16, 1, 4, 2)
  )

p_seq <- ggplot(id_long, aes(x = xnum, y = row_id, fill = pct)) +
  geom_tile(width = 0.94, height = 0.90, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = txt, colour = pct >= 60), size = 2.35, family = "Arial") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#1a1a1a"), guide = "none") +
  scale_fill_gradientn(
    colours = c("#f7fbff", "#deebf7", "#9ecae1", "#3182bd", "#084594"),
    limits = c(0, 100), breaks = c(0, 50, 100),
    name = "N-term identity (%)"
  ) +
  scale_x_continuous(
    breaks = c(1:5, (6:10) + 0.35),
    labels = seq_xlabs,
    expand = c(0.02, 0.02),
    position = "top"
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(sp_fac ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL,
       title = expression(bolditalic("B. dorsalis") * bold(" (X1–X5)   |   ") * bolditalic("D. melanogaster"))) +
  theme_minimal(base_size = 8.5, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 8.5, hjust = 0.5, margin = margin(b = 3)),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x.top = element_text(size = 7.6, face = "bold", colour = "#111111", margin = margin(b = 2)),
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.spacing.y = unit(4.5, "pt"),
    legend.position = "bottom",
    legend.key.height = unit(3.2, "mm"),
    legend.key.width  = unit(8.0, "mm"),
    legend.title = element_text(size = 7.8, face = "bold", margin = margin(r = 6)),
    legend.text  = element_text(size = 7.4),
    plot.margin = margin(2, 6, 4, 1)
  )

# ------------------------------------------------------------------------------
# Panel a-right: Shared-scale multi-tissue expression heatmap (nF | mF across 6 tissues)
# ------------------------------------------------------------------------------
tis_ord <- c("Ovipositor", "Foreleg", "Mouthpart", "Antenna", "Eye", "Brain")

b_sub <- bact[isoform %in% row_spec$iso & stage %in% c("newly emerged female", "mature female", "adult female (public data)")]
b_sub[, stg := ifelse(stage == "newly emerged female", "nF", "mF")]
b_agg <- b_sub[, .(mean_tpm = mean(tpm), log2_tpm = mean(log2(tpm + 1))), by = .(sp = species, iso = isoform, tissue, stg)]

fly_f <- fly[sex == "F"]
fly_f[, tissue := c(OV = "Ovipositor", LE = "Foreleg", PR = "Mouthpart", AN = "Antenna")[tissue]]
fly_f[, iso := sub("-R", ".", transcript)]
fly_agg <- fly_f[, .(sp = "Dmel", mean_tpm = mean(tpm), log2_tpm = mean(log2(tpm + 1)), stg = "mF"), by = .(iso, tissue)]

expr_obs <- rbind(b_agg, fly_agg[, .(sp, iso, tissue, stg, mean_tpm, log2_tpm)])

col_grid <- rbindlist(lapply(tis_ord, function(ti) data.table(tissue = ti, stg = c("nF", "mF"))))
col_grid[, col_idx := seq_len(.N)]
col_grid[, xnum := col_idx + (match(tissue, tis_ord) - 1) * 0.35]

expr_full <- CJ(iso = row_spec$iso, col_idx = col_grid$col_idx, sorted = FALSE)
expr_full <- merge(expr_full, col_grid, by = "col_idx", sort = FALSE)
expr_full <- merge(expr_full, row_spec[, .(iso, sp, sp_fac, row_id)], by = "iso", sort = FALSE)
expr_full <- merge(expr_full, expr_obs[, .(iso, tissue, stg, mean_tpm, log2_tpm)],
                   by = c("iso", "tissue", "stg"), all.x = TRUE, sort = FALSE)
expr_full[, val_cap := pmin(log2_tpm, 4.5)]
expr_full[, tpm_lbl := ifelse(!is.na(mean_tpm) & mean_tpm >= 0.8, sprintf("%.1f", mean_tpm), "")]

tis_hdr_pos <- col_grid[, .(x = mean(xnum)), by = tissue]

p_expr <- ggplot(expr_full, aes(x = xnum, y = row_id)) +
  geom_tile(aes(fill = val_cap), width = 0.94, height = 0.90, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = tpm_lbl, colour = val_cap >= 2.6), size = 2.25, family = "Arial") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#1a1a1a"), guide = "none") +
  scale_fill_gradientn(
    colours = c("#fffcf7", "#fee8c8", "#fdbb84", "#e34a33", "#7f0000"),
    na.value = "#eaeaea",
    limits = c(0, 4.5), breaks = c(0, 1.5, 3, 4.5), labels = c("0", "1.5", "3.0", "≥4.5"),
    name = expression(bold(log[2]*"(TPM + 1)"))
  ) +
  scale_x_continuous(
    breaks = col_grid$xnum,
    labels = col_grid$stg,
    expand = c(0.02, 0.02),
    position = "top",
    sec.axis = dup_axis(breaks = tis_hdr_pos$x, labels = tis_hdr_pos$tissue)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(sp_fac ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL, title = "Multi-tissue female expression (TPM ≥ 0.8 shown; grey = unsampled)") +
  theme_minimal(base_size = 8.5, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5, margin = margin(b = 3)),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x.top = element_text(size = 7.2, colour = "#222222", margin = margin(b = 1)),
    axis.text.x.bottom = element_text(size = 8.0, face = "bold", colour = "#111111", margin = margin(t = 2)),
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.spacing.y = unit(4.5, "pt"),
    legend.position = "bottom",
    legend.key.height = unit(3.2, "mm"),
    legend.key.width  = unit(9.0, "mm"),
    legend.title = element_text(size = 7.8, margin = margin(r = 6)),
    legend.text  = element_text(size = 7.4),
    plot.margin = margin(2, 4, 4, 2)
  )

# ------------------------------------------------------------------------------
# Panel b (left): Ovipositor Gr28b isoform TPM stacked bars (nF vs mF)
# Panel c (right): Mature ovipositor relative isoform composition (%)
# Full italicized species names on x-axis (ZERO abbreviations like B. min. or B. cor.!)
# ------------------------------------------------------------------------------
ovi_dt <- merge(expr_obs[tissue == "Ovipositor"], row_spec[, .(iso, sp_lab, cls_fac)], by = "iso")
ovi_sum <- ovi_dt[, .(tpm = sum(mean_tpm)), by = .(sp_lab, stg, cls_fac)]
ovi_sum[, sp_fac := factor(sp_lab, levels = sp_order)]
ovi_sum[, stg := factor(stg, levels = c("nF", "mF"))]

# Create grouped x positions on Panel b so full italic species names sit cleanly at 32 deg below!
sp_x_dt <- data.table(sp_fac = factor(sp_order, levels = sp_order), sp_i = seq_along(sp_order))
ovi_sum <- merge(ovi_sum, sp_x_dt, by = "sp_fac")
# Offset nF slightly left (-0.20) and mF slightly right (+0.20); for species with only mF (Btsu, Dmel), center mF at sp_i
has_nf <- unique(ovi_sum[stg == "nF"]$sp_fac)
ovi_sum[, xpos := ifelse(sp_fac %in% has_nf, sp_i + ifelse(stg == "nF", -0.21, 0.21), sp_i)]

ovi_tot <- ovi_sum[, .(tot = sum(tpm)), by = .(sp_fac, sp_i, stg, xpos)]
ovi_dom <- ovi_sum[stg == "mF", .SD[which.max(tpm)], by = sp_fac]
ovi_dom <- merge(ovi_dom, ovi_tot[stg == "mF"], by = c("sp_fac", "sp_i", "stg", "xpos"))
ovi_dom[, pct := ifelse(tot >= 0.1, round(100 * tpm / tot), NA_real_)]
ovi_dom[, short_cls := ifelse(cls_fac == "X1 / Gr28b.C",
                              ifelse(sp_fac == "D. melanogaster", "Gr28b.C", "X1"),
                              ifelse(cls_fac == "X3 / X3-like", "X3-like", as.character(cls_fac)))]
ovi_dom[, lbl := ifelse(tot < 0.1, sprintf("%.2f", tot), sprintf("%.1f\n%d%% %s", tot, pct, short_cls))]

# Small stage labels (nF / mF) just below baseline (y = -0.6) or inside/above bars
stg_lbl_dt <- unique(ovi_sum[, .(sp_fac, stg, xpos)])

p_ovi_abs <- ggplot(ovi_sum, aes(x = xpos, y = tpm, fill = cls_fac)) +
  geom_col(width = 0.36, colour = "white", linewidth = 0.28) +
  geom_text(data = ovi_dom, aes(x = xpos, y = tot + 0.55, label = lbl),
            inherit.aes = FALSE, family = "Arial", size = 2.22, lineheight = 0.90, vjust = 0, colour = "#111111") +
  geom_text(data = stg_lbl_dt, aes(x = xpos, y = -0.75, label = stg),
            inherit.aes = FALSE, family = "Arial", size = 2.35, colour = "#444444", vjust = 1) +
  scale_fill_manual(values = cls_cols, name = "Isoform class", drop = FALSE) +
  scale_x_continuous(breaks = sp_x_dt$sp_i, labels = sp_order, expand = expansion(add = c(0.45, 0.55))) +
  scale_y_continuous(limits = c(-2.2, 24.5), breaks = c(0, 5, 10, 15, 20), expand = c(0, 0)) +
  geom_hline(yintercept = 0, linewidth = 0.45, colour = "#111111") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Ovipositor TPM",
       title = expression(bold("Ovipositor ") * bolditalic("Gr28b") * bold(" expression (nF vs mF)"))) +
  theme_classic(base_size = 8.5, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 9.0, hjust = 0.5, margin = margin(b = 5)),
    axis.text.x = element_text(size = 8.2, face = "bold.italic", angle = 30, hjust = 1, vjust = 1, colour = "#111111",
                               margin = margin(t = 2)),
    axis.text.y = element_text(size = 8.0, colour = "#111111"),
    axis.title.y = element_text(size = 8.5, face = "bold", colour = "#111111"),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.45),
    axis.ticks.y = element_line(linewidth = 0.4),
    legend.position = "bottom",
    legend.key.size = unit(3.5, "mm"),
    legend.title = element_text(size = 8.0, face = "bold"),
    legend.text = element_text(size = 7.6),
    plot.margin = margin(6, 12, 2, 6)
  ) +
  guides(fill = guide_legend(nrow = 2))

# Panel c: Relative composition (%) in mature ovipositors (full italic species names!)
ovi_rel <- copy(ovi_sum[stg == "mF"])
ovi_rel <- merge(ovi_rel, ovi_tot[stg == "mF", .(sp_fac, tot)], by = "sp_fac")
ovi_rel[, pct := ifelse(tot >= 0.1, 100 * tpm / tot, 0)]
undet_df <- ovi_tot[stg == "mF" & tot < 0.1]

p_ovi_rel <- ggplot(ovi_rel, aes(x = sp_fac, y = pct, fill = cls_fac)) +
  geom_col(width = 0.66, colour = "white", linewidth = 0.28) +
  geom_text(data = undet_df, aes(x = sp_fac, y = 50, label = "n.d.\n(<0.1)"),
            inherit.aes = FALSE, family = "Arial", size = 2.55, colour = "#555555", lineheight = 0.92) +
  scale_fill_manual(values = cls_cols, guide = "none", drop = FALSE) +
  scale_y_continuous(limits = c(0, 105), breaks = c(0, 25, 50, 75, 100), expand = c(0, 0)) +
  labs(x = NULL, y = expression(bold("Share of mature ovipositor ") * bolditalic("Gr28b") * bold(" (%)")),
       title = "Mature ovipositor isoform share") +
  theme_classic(base_size = 8.5, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 9.0, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    axis.text.x = element_text(size = 8.2, face = "bold.italic", angle = 30, hjust = 1, vjust = 1, colour = "#111111"),
    axis.text.y = element_text(size = 8.0, colour = "#111111"),
    axis.title.y = element_text(size = 8.2, colour = "#111111"),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.4),
    plot.margin = margin(6, 6, 2, 6)
  )

# Assemble complete single-page Supplementary Fig. 12 (generous 9.6 x 11.8 in canvas -> crisp 8-9 pt fonts!)
top_row <- (p_cls | p_seq | p_expr) + plot_layout(widths = c(0.075, 1.00, 1.48))
bot_row <- (p_ovi_abs | p_ovi_rel) + plot_layout(widths = c(1.38, 1.00))

out_pdf <- file.path(OUT_DIR, "Supplementary_Figure_12_Gr28b_isoform_RNAseq.pdf")
out_png <- file.path(OUT_DIR, "Supplementary_Figure_12_Gr28b_isoform_RNAseq.png")

render_sf12 <- function() {
  grid.newpage()
  pushViewport(viewport(x = 0.5, y = 0.5, width = 0.985, height = 0.985))
  pushViewport(viewport(x = 0.008, y = 0.995, width = 0.985, height = 0.705, just = c("left", "top")))
  print(top_row, newpage = FALSE)
  popViewport()
  grid.text("a", x = unit(0.004, "npc"), y = unit(0.996, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 13, fontface = "bold", fontfamily = "Arial"))

  pushViewport(viewport(x = 0.008, y = 0.285, width = 0.985, height = 0.282, just = c("left", "top")))
  print(bot_row, newpage = FALSE)
  popViewport()
  grid.text("b", x = unit(0.004, "npc"), y = unit(0.286, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 13, fontface = "bold", fontfamily = "Arial"))
  grid.text("c", x = unit(0.576, "npc"), y = unit(0.286, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 13, fontface = "bold", fontfamily = "Arial"))
  popViewport()
}

cairo_pdf(out_pdf, width = 9.8, height = 11.8, family = "Arial")
render_sf12()
dev.off()

png(out_png, width = 9.8, height = 11.8, units = "in", res = 300, type = "cairo", family = "Arial")
render_sf12()
dev.off()

cat("Wrote enlarged Supplementary Fig. 12 ->", out_pdf, "\n")
