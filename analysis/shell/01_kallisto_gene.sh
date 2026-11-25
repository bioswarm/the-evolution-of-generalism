#!/bin/bash
# Gene-level quantification (Methods): kallisto v0.51.1, one representative CDS per gene, 100 bootstraps.
# Output layout expected by R/00_common.R: <kallisto_gene>/<Species>/<sample>/abundance.{h5,tsv}
# Usage: shell/01_kallisto_gene.sh <out_dir> [threads]
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
CFG=$ROOT/config.yaml
get() { grep -E "^$1:" "$CFG" | head -1 | sed -E "s/^$1: *//; s/ *#.*//; s/\"//g"; }
LISTS=$(get reads_lists); CDS=$(get cds_fasta)
OUT=${1:?output directory}; T=${2:-48}
for sp in Bcor Bdor Bmin Zcuc Ztau; do
  mkdir -p "$OUT/$sp"
  [ -f "$OUT/$sp/$sp.idx" ] || kallisto index -i "$OUT/$sp/$sp.idx" "$CDS/$sp.cds.fa"
  paste "$LISTS/$sp.reads.R1.list" "$LISTS/$sp.reads.R2.list" | while read -r r1 r2; do
    s=$(basename "$r1" .R1.fastq.gz)
    [ -f "$OUT/$sp/$s/abundance.h5" ] && continue
    kallisto quant -i "$OUT/$sp/$sp.idx" -o "$OUT/$sp/$s" -b 100 -t "$T" "$r1" "$r2"
  done
done
