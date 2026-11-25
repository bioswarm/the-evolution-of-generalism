#!/bin/bash
# All-transcript CDS reference for the Gr28b isoform analysis: gffread -x (genome + annotation) per species.
# B. dorsalis, B. minax, Z. cucurbitae: original annotation (Gr28b models unchanged in the re-annotation of 2026-11-04).
# B. correcta, Z. tau, B. tsuneonis: re-annotated Gr28b (Iso-Seq / RNA-seq-supported alternative 5' exons):
#   Bcor complete.genomic.gtf (Bc_003035-R1, -R2), Ztau step1.gff3 (Zt1G00596.1-.10), Btsu Bt.gff + Gr28b_R2.gff3 (Bt_011990-R1, -R2).
# Usage: shell/03a_isoform_reference.sh <out_dir>   (writes <out_dir>/ref/<Species>.alltx.cds.fa)
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
get() { grep -E "^$1:" "$ROOT/config.yaml" | head -1 | sed -E "s/^$1: *//; s/ *#.*//; s/\"//g"; }
G=$(get genesets); N=$(get gr28b_reannotation)
O=${1:?output directory}/ref; mkdir -p "$O"
for p in "Bdor Bdor.gtf" "Bmin Bmin.gtf" "Zcuc Bcuc.gtf"; do set -- $p
  gffread -x "$O/$1.alltx.cds.fa" -g "$G/${2%.*}.fa" "$G/$2"; done
grep -v 'transcript_id ""' "$N/Bcor/OGS/complete.genomic.gtf" > "$O/Bcor.complete.tx.gtf"          # drop gene-level records
gffread -x "$O/Bcor.alltx.cds.fa" -g "$(readlink -f "$N/Bcor/OGS/Bcor_final.fa")" "$O/Bcor.complete.tx.gtf"
gffread -x "$O/Ztau.alltx.cds.fa" -g "$(readlink -f "$N/Btau/OGS/Gr28b/Ztau.fa")" "$N/Btau/OGS/final/step1.gff3"
cat "$N/Btsu/OGS/Bt.gff" <(grep -v '^#' "$N/Btsu/OGS/Gr28b_AS/Gr28b_R2.gff3") > "$O/Btsu.gr28bR2.gff"
gffread -x "$O/Btsu.alltx.cds.fa" -g "$N/Btsu/OGS/Btsu.fa" "$O/Btsu.gr28bR2.gff"
for s in Bdor Bmin Zcuc Bcor Ztau Btsu; do echo "$s $(grep -c '>' "$O/$s.alltx.cds.fa") transcripts"; done
