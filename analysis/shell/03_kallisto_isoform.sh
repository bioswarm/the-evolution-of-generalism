#!/bin/bash
# Gr28b isoform analysis (R/26): kallisto v0.51.1 with the Methods settings (-b 100) against an index of ALL annotated
# transcripts (built by shell/03a_isoform_reference.sh); all six species (B. tsuneonis: public RNA-seq, btsu_reads).
# Usage: shell/03_kallisto_isoform.sh <out_dir> [parallel_jobs] [threads_per_job] [species ...]   (default: all six)
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
get() { grep -E "^$1:" "$ROOT/config.yaml" | head -1 | sed -E "s/^$1: *//; s/ *#.*//; s/\"//g"; }
LISTS=$(get reads_lists)
W=${1:?output directory}; P=${2:-4}; T=${3:-12}; shift 3 2>/dev/null || true
SP=${*:-Bdor Bmin Zcuc Ztau Bcor Btsu}; BT=$(get btsu_reads)
O=$W/kallisto; mkdir -p "$O"
for sp in $SP; do
  [ -f "$O/$sp.idx" ] || kallisto index -i "$O/$sp.idx" -t 8 "$W/ref/$sp.alltx.cds.fa"
done
jobs=$O/jobs.txt; : > "$jobs"
for sp in $SP; do
  if [ "$sp" = Btsu ]; then for r1 in "$BT"/Bt-*.R1.fastq.gz; do echo "$sp $(basename "$r1" .R1.fastq.gz) $r1 ${r1/.R1./.R2.}" >> "$jobs"; done
  else paste "$LISTS/$sp.reads.R1.list" "$LISTS/$sp.reads.R2.list" | while read -r r1 r2; do
    echo "$sp $(basename "$r1" .R1.fastq.gz) $r1 $r2" >> "$jobs"; done; fi
done
run() { sp=$1; s=$2; [ -f "$O/$sp/$s/abundance.h5" ] && return 0; mkdir -p "$O/$sp"
  kallisto quant -i "$O/$sp.idx" -o "$O/$sp/$s" -b 100 -t "$T" "$3" "$4" > "$O/$sp/$s.log" 2>&1 && echo "DONE $sp $s" || echo "FAIL $sp $s"; }
export -f run; export O T
xargs -P "$P" -L 1 bash -c 'run "$@"' _ < "$jobs"
echo "FINISHED $(ls "$O"/*/*/abundance.h5 2>/dev/null | wc -l) of $(wc -l < "$jobs")"
