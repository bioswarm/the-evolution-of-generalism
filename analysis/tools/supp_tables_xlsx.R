#!/usr/bin/env Rscript
# Revised copy of the supplementary tables workbook (Supplementary Tables 1 and 3; the original in paper/ is never modified).
# Changes (changed cells filled yellow; a note row at the end of Supplementary Table 1 explains them):
#  1. B. correcta male foreleg / mouthpart: the raw files Bc-MM-L1-3 and Bc-MM-M1-3 carry exchanged labels (expression-based
#     identity check, R/01; corrected at load time via config.yaml `relabel`). Rows keep their sample ID and tissue and now
#     carry the read counts of the files that contain that tissue; column I names those raw files.
#  2. 'Lavae' -> 'Larvae'.   3. Hi-C %Mapped 77.599999999999994 -> 77.60.
#  4. ATAC: the B. dorsalis newly emerged female ovipositor library (Fig. 4, Supplementary Data 7) is added, with a Stage column.
# Usage: Rscript tools/supp_tables_xlsx.R [output.xlsx]   (paths from config.yaml: paper_dir, revised_dir, atac_root, atac)
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "..", "R", "00_common.R"))
suppressPackageStartupMessages({ library(openxlsx); library(jsonlite) })
SRC <- file.path(CFG$paper_dir, "supplemental_tables_Gr28b_20260306.xlsx")
DST <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else file.path(CFG$revised_dir, "supplemental_tables_Gr28b_revised.xlsx")   # optional output path
SRC_MD5 <- "694cf1df06d55f1e95b77ff1511918b2"
stopifnot(unname(tools::md5sum(SRC)) == SRC_MD5)

wb <- loadWorkbook(SRC); sh <- 1
x <- read.xlsx(SRC, sheet = sh, colNames = FALSE, skipEmptyRows = FALSE, skipEmptyCols = FALSE)
cell <- function(r, c) { v <- x[r, c]; if (is.null(v) || is.na(v)) NA else trimws(as.character(v)) }
YEL <- createStyle(fgFill = "#FFFF00"); BOLD <- createStyle(textDecoration = "bold"); ITAL <- createStyle(textDecoration = "italic", fgFill = "#FFFF00")
changed <- 0L
put <- function(r, c, v, style = YEL) { writeData(wb, sh, v, startRow = r, startCol = c, colNames = FALSE)
  addStyle(wb, sh, style, rows = r, cols = c, stack = TRUE); changed <<- changed + 1L }

# ---- 1. Bc-MM L <-> M ----
id <- vapply(seq_len(nrow(x)), function(r) cell(r, 1), "")
row_of <- setNames(seq_along(id), id)
orig <- lapply(setNames(nm = grep("^Bc-MM-[LM][123]$", id, value = TRUE)), function(s) as.numeric(c(x[row_of[[s]], 7], x[row_of[[s]], 8])))
put(3, 9, "Raw files (if different from sample ID)", createStyle(textDecoration = "bold", fgFill = "#FFFF00"))
for (k in 1:3) for (pr in list(c(sprintf("Bc-MM-L%d", k), sprintf("Bc-MM-M%d", k)), c(sprintf("Bc-MM-M%d", k), sprintf("Bc-MM-L%d", k)))) {
  r <- row_of[[pr[1]]]; v <- orig[[pr[2]]]
  put(r, 7, v[1]); put(r, 8, v[2]); put(r, 9, sprintf("%s.R1/R2.fastq.gz (labels exchanged; see note)", pr[2])) }

# ---- 2. Lavae -> Larvae ----
for (r in seq_len(nrow(x))) for (c in seq_len(ncol(x))) if (identical(cell(r, c), "Lavae")) put(r, c, "Larvae")

# ---- 3. Hi-C ----
hic <- which(id == "Hi-C sequencing data"); stopifnot(length(hic) == 1)
r <- hic + 3; put(r, 3, round(as.numeric(x[r, 3]), 2))          # 77.599999999999994 -> 77.6 (shown 77.60)
addStyle(wb, sh, createStyle(numFmt = "0.00"), rows = (hic + 2):(hic + 3), cols = c(3, 8, 9), gridExpand = TRUE, stack = TRUE)

# ---- 4. ATAC ----
at <- grep("^ATAC sequencing data", id); stopifnot(length(at) == 1)
qc <- function(d) { f <- sort(Sys.glob(file.path(CFG$atac_root, d, "atac/*/call-qc_report/execution/qc.json")))[1]
  s <- fromJSON(f)$align$samstat$rep1; c(s$total_reads, round(s$pct_mapped_reads, 1)) }
put(at, 1, "ATAC sequencing data for ovipositors")
hdr <- c("Species", "Stage", "#Reads", "% Mapped reads"); for (i in seq_along(hdr)) put(at + 1, i, hdr[i], createStyle(textDecoration = "bold", fgFill = "#FFFF00"))
lib <- list(c("Bactrocera minax", "Sexually mature female", CFG$atac$Bmin_MO), c("Bactrocera dorsalis", "Sexually mature female", CFG$atac$Bdor_MO),
            c("Bactrocera dorsalis", "Newly emerged female", CFG$atac$Bdor_NO))
for (i in seq_along(lib)) { q <- qc(lib[[i]][3]); r <- at + 1 + i
  put(r, 1, lib[[i]][1]); put(r, 2, lib[[i]][2]); put(r, 3, q[1]); put(r, 4, q[2]) }

# ---- note ----
note <- at + 1 + length(lib) + 2
put(note, 1, paste("Yellow cells were changed in the revision. B. correcta male foreleg (L) and mouthpart (M) libraries: the raw files",
  "Bc-MM-L1-3 and Bc-MM-M1-3 carry exchanged labels, identified from expression of tissue markers (Methods); read counts are now",
  "listed under the correct tissue and column I gives the raw file names. The B. dorsalis newly emerged ovipositor ATAC library",
  "(Fig. 4, Supplementary Data 7) was added."), ITAL)
saveWorkbook(wb, DST, overwrite = TRUE)
# The original workbook references drawing / vmlDrawing parts (cell comments) that openxlsx does not write back, which leaves
# dangling relationships (Excel asks to repair the file; openpyxl cannot open it). Remove every relationship whose target file is
# missing and the matching [Content_Types] overrides; sheet data and styles are untouched.
local({
  DSTA <- normalizePath(DST); td <- tempfile("xlsxfix"); dir.create(td); on.exit(unlink(td, recursive = TRUE), add = TRUE)
  utils::unzip(DSTA, exdir = td); have <- list.files(td, recursive = TRUE, all.files = TRUE)
  for (rf in grep("\\.rels$", have, value = TRUE)) {
    f <- file.path(td, rf); x <- paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    srcdir <- dirname(dirname(rf))                              # folder of the part that owns the .rels file
    for (r in regmatches(x, gregexpr("<Relationship\\b[^>]*/>", x))[[1]]) {
      if (grepl("TargetMode=\"External\"", r, fixed = TRUE)) next
      tg <- sub('.*Target="([^"]+)".*', "\\1", r)
      path <- if (startsWith(tg, "/")) file.path(td, sub("^/", "", tg)) else file.path(td, srcdir, tg)
      if (!file.exists(path)) { x <- sub(r, "", x, fixed = TRUE); log_msg(sprintf("removed dangling relationship in %s -> %s", rf, tg)) }
    }
    writeLines(x, f, useBytes = TRUE)
  }
  ct <- file.path(td, "[Content_Types].xml"); x <- paste(readLines(ct, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  for (o in regmatches(x, gregexpr("<Override\\b[^>]*/>", x))[[1]]) {
    pn <- sub('.*PartName="/?([^"]+)".*', "\\1", o); if (!file.exists(file.path(td, pn))) x <- sub(o, "", x, fixed = TRUE)
  }
  writeLines(x, ct, useBytes = TRUE)
  owd <- setwd(td); on.exit(setwd(owd), add = TRUE); unlink(DSTA)
  ord <- c("[Content_Types].xml", setdiff(list.files(".", recursive = TRUE, all.files = TRUE), "[Content_Types].xml"))
  utils::zip(DSTA, ord, flags = "-X -q", zip = Sys.which("zip"))
})
stopifnot(unname(tools::md5sum(SRC)) == SRC_MD5)
log_msg(sprintf("wrote %s; %d cells written; original unchanged", DST, changed))
