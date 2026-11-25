#!/usr/bin/env Rscript
# Species-tree supermatrix (Methods; called by shell/06): concatenate the trimmed per-locus protein alignments that
# contain exactly <n_species> sequences, in file-name order, one record per species (sequence ids = species names).
# Usage: Rscript tools/concatenate_loci.R <dir with *.trim or *.muscle alignments> <n_species> <out.fasta>
a <- commandArgs(TRUE); stopifnot(length(a) == 3); n <- as.integer(a[2])
fs <- sort(list.files(a[1], pattern = "\\.(trim|muscle|aln)$", full.names = TRUE)); seqs <- list(); used <- 0L
for (f in fs) {
  l <- readLines(f); h <- grep("^>", l); if (length(h) != n) next
  id <- sub("^>(\\S+).*", "\\1", l[h]); s <- vapply(seq_along(h), function(i) gsub("\\s", "", paste(l[(h[i] + 1):(c(h[-1], length(l) + 1)[i] - 1)], collapse = "")), "")
  for (i in seq_along(id)) seqs[[id[i]]] <- c(seqs[[id[i]]], s[i]); used <- used + 1L
}
id <- sort(names(seqs)); writeLines(unlist(lapply(id, function(k) c(paste0(">", k), paste(seqs[[k]], collapse = "")))), a[3])
message(used, " loci concatenated, ", length(id), " species, ", nchar(paste(seqs[[id[1]]], collapse = "")), " columns")
