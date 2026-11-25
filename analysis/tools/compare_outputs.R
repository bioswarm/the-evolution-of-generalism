#!/usr/bin/env Rscript
# Compare the results/ and figures/ trees of two output roots (e.g. a frozen reference and a fresh notebook run).
# Usage: Rscript tools/compare_outputs.R <ref_root> <new_root> [regex ...]   (only paths matching any regex; default all)
# Tables: identical bytes, else equal values (tolerance 1e-8); PNG/PDF: identical pixels (PDF rasterised at 100 dpi);
# .rds: bytes only, a difference is a WARN (Seurat objects carry timestamps); anything else: bytes.
# Exit status 1 if any file is DIFF or MISSING (in ref, not in new).
suppressPackageStartupMessages({ library(data.table); library(magick); library(tools) })
a <- commandArgs(TRUE); stopifnot(length(a) >= 2)
ref <- normalizePath(a[1]); new <- normalizePath(a[2], mustWork = FALSE); pat <- if (length(a) > 2) a[-(1:2)] else "."
rel <- function(root) { f <- unlist(lapply(c("results", "figures"), function(d) if (dir.exists(file.path(root, d)))
  file.path(d, list.files(file.path(root, d), recursive = TRUE)))); f[Reduce(`|`, lapply(pat, grepl, x = f))] }
files <- sort(rel(ref))
pixels_equal <- function(f1, f2, pdf) {
  r <- function(f) if (pdf) image_read(f, density = 100) else image_read(f)
  i1 <- r(f1); i2 <- r(f2)
  if (length(i1) != length(i2)) return(sprintf("%d vs %d pages", length(i1), length(i2)))
  for (k in seq_along(i1)) {
    if (!identical(image_info(i1[k])[c("width", "height")], image_info(i2[k])[c("width", "height")])) return("size differs")
    d <- image_compare_dist(i1[k], i2[k], metric = "AE")$distortion
    if (d > 0) return(sprintf("%g pixels differ (page %d)", d, k)) }
  TRUE }
check <- function(f) {
  f1 <- file.path(ref, f); f2 <- file.path(new, f)
  if (!file.exists(f2)) return(c("MISSING", ""))
  if (md5sum(f1) == md5sum(f2)) return(c("identical", ""))
  ext <- tolower(file_ext(f))
  if (ext %in% c("tsv", "txt", "csv", "bedgraph")) {
    x <- tryCatch(all.equal(fread(f1, sep = "\t", header = TRUE), fread(f2, sep = "\t", header = TRUE), check.attributes = FALSE,
       tolerance = 1e-8), error = function(e) conditionMessage(e))
    return(if (isTRUE(x)) c("equal-values", "") else c("DIFF", paste(head(x, 2), collapse = "; "))) }
  if (ext %in% c("png", "pdf")) { x <- pixels_equal(f1, f2, ext == "pdf"); return(if (isTRUE(x)) c("pixel-identical", "") else c("DIFF", x)) }
  if (ext == "rds") return(c("WARN", "rds bytes differ"))
  c("DIFF", "bytes differ") }
res <- rbindlist(lapply(files, function(f) { r <- check(f); data.table(status = r[1], file = f, note = r[2]) }))
extra <- setdiff(rel(new), files)
if (nrow(res)) print(res[, .N, by = status])
bad <- res[status %in% c("DIFF", "MISSING", "WARN")]
if (nrow(bad)) { options(width = 250); print(bad, nrows = 500) }
if (length(extra)) cat("only in new:", paste(extra, collapse = ", "), "\n")
cat(sprintf("%d files compared, %d DIFF, %d MISSING\n", nrow(res), sum(res$status == "DIFF"), sum(res$status == "MISSING")))
quit(status = as.integer(any(res$status %in% c("DIFF", "MISSING"))))
