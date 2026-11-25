#!/bin/bash
# B. tsuneonis public RNA-seq (SRR33000917-SRR33000937; adult female/male antenna, brain, foreleg, ovipositor):
# kallisto v0.51.1, same settings as the Methods (-b 100). Used for baseline expression only (R/11).
# Usage: shell/02_kallisto_Btsu.sh <out_dir> [threads]
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
get() { grep -E "^$1:" "$ROOT/config.yaml" | head -1 | sed -E "s/^$1: *//; s/ *#.*//; s/\"//g"; }
R=$(get btsu_reads); CDS=$(get btsu_cds)
O=${1:?output directory}; T=${2:-48}; mkdir -p "$O"
[ -f "$O/Btsu.idx" ] || kallisto index -i "$O/Btsu.idx" "$CDS"
for r1 in "$R"/Bt-*.R1.fastq.gz; do
  s=$(basename "$r1" .R1.fastq.gz)
  [ -f "$O/$s/abundance.h5" ] && continue
  kallisto quant -i "$O/Btsu.idx" -o "$O/$s" -b 100 -t "$T" "$r1" "${r1/.R1./.R2.}"
done
