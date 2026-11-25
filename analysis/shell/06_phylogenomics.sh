#!/bin/bash
# Species tree and aBSREL (Methods). NOT re-run for the revision (results reused; see README "What was re-run").
#   per locus: MAFFT v7.526 (default) -> trimAl v1.4 (-automated1) -> concatenate loci present in all 13 species
#   species tree: IQ-TREE v1.6.12, -m JTT+F+R10 -b 100
#   aBSREL: PAL2NAL codon alignment -> HyPhy 2.5 absrel --branches test (Bactrocera_minax{test}) per orthogroup
# Usage: shell/06_phylogenomics.sh <dir with fasta/*.fa per orthogroup> <out_dir> [threads]
set -euo pipefail
IN=${1:?}; O=${2:?}; T=${3:-20}; mkdir -p $O/aln $O/trim
for f in $IN/fasta/*.fa; do id=$(basename $f .fa)
  mafft --quiet $f > $O/aln/$id.aln; trimal -in $O/aln/$id.aln -out $O/trim/$id.trim -automated1; done
Rscript "$(dirname "$0")/../tools/concatenate_loci.R" $O/trim 13 $O/concatenated.fasta
iqtree -s $O/concatenated.fasta -m JTT+F+R10 -b 100 -nt $T -pre $O/species
