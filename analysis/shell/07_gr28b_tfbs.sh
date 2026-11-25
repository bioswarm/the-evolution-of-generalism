#!/bin/bash
# TF-binding-site prediction at the ATAC peak (element E) upstream of the Gr28b
# C-clade first exon. Usage: shell/07_gr28b_tfbs.sh   (run from reproduce/)
# Outputs -> results/locus/tfbs/
#   accessibility.tsv          mean fc bigWig over E and the 428-bp control window (Bdor MO vs NO)
#   upstream_Cexon.fa          1,500 bp upstream + first 100 bp of the C-clade first exon (transcription orientation)
#   upstream_Cexon.mafft.fa    MAFFT E-INS-i alignment of the six regions
#   bg_<sp>.txt                order-1 Markov background (Gr28b gene +/- 100 kb)
#   fimo/<sp>_{1e-4,1e-3}/     FIMO scans (JASPAR2026 CORE insects, both strands)
#   jaspar_insects.meme        concatenated motif file used for scanning
#   Gr28b_Cexon_upstream_TFBS.tsv   one row per Bdor hit in E + per-species status + expression
#   Gr28b_Cexon_upstream_TFBS_counts.tsv    candidate-rule counts in E vs control window
#   sites_on_aln.tsv           every same-motif hit mapped to alignment columns (for the figure)
#   snrna_detection.tsv        % of Gr28b+ ovipositor nuclei expressing each Bd TF gene
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=results/locus/tfbs; mkdir -p "$OUT" "$OUT/fimo"
CHIP=/path/to/conda/envs/chipseq/bin/python
MEME=/path/to/conda/envs/anvio-9/bin
MAFFT=/path/to/conda/envs/genera/bin/mafft
GS=/path/to/data/genesets
JZIP=/path/to/project/TF/JASPAR_data/JASPAR2026_CORE_insects_non-redundant_pfms_meme.zip

# ---- element + control window (Bdor, + strand) ----------------------------------------------
E=NC_064303.1:108181027-108181454            # relaxed MACS2 peak, 38 bp upstream of exon 1
CT=NC_064303.1:108180599-108181026           # 428-bp control window immediately upstream of E
# in the 1600-bp extracted region (1500 up + 100 exon), 1-based:
E_RG=1036-1463;  CT_RG=608-1035

# ---- species table: name genome_fa chrom strand exon_start(transcription) gene_start gene_end
# exon_start = first-exon start in transcription orientation (max coord of the 5' segment for -).
cat > "$OUT/species.tsv" <<EOF
Bdor	$GS/Bdor.fa	NC_064303.1	+	108181492	108169689	108194489
Bmin	$GS/Bmin.fa	BmChr1	+	9498362	9477063	9508221
Zcuc	$GS/Bcuc.fa	NC_071668.1	+	10588317	10578893	10600247
Bcor	$GS/Bcor/genome_aln/final/Bcor_final.fa	Bcor_chr1	-	8407930	8395207	8419583
Btsu	$GS/Btsu/genome/ncbi/Btsu.fa	chr1	-	53408606	53398861	53408606
Ztau	/path/to/data/genesets/Btau/RNAseq/RSEM/Btau.fa	chr1	-	68077885	68068077	68113536
EOF

# ---- 1. accessibility at E vs control (fc bigWigs, mature vs newly-emerged ovipositor) -------
MO=$(ls /public_data/ATAC-seq/deprecate/Bdor/Wuhan/Mo/atac/*/call-macs2_signal_track/shard-0/execution/Bd_MO.*fc.signal.bigwig | head -1)
NO=$(ls /public_data/ATAC-seq/deprecate/Bdor/Wuhan/No/atac/*/call-macs2_signal_track/shard-0/execution/Bd_NO.*fc.signal.bigwig | head -1)
{ echo -e "track\tregion\tmean_fc\tn_bases";
  $CHIP shell/tfbs_bwmean.py "$MO" "$E" "$CT" | awk -v t=Bdor_MO 'BEGIN{OFS="\t"}{print t,$1,$2,$3}';
  $CHIP shell/tfbs_bwmean.py "$NO" "$E" "$CT" | awk -v t=Bdor_NO 'BEGIN{OFS="\t"}{print t,$1,$2,$3}';
} > "$OUT/accessibility.tsv"

# ---- 2. extract upstream regions + align ----------------------------------------------------
> "$OUT/upstream_Cexon.fa"; > "$OUT/locus_100kb_regions.txt"
while IFS=$'\t' read -r sp fa chr strand start gs ge; do
  [ -f "$fa.fai" ] || samtools faidx "$fa"
  if [ "$strand" = "+" ]; then
    lo=$((start-1500)); hi=$((start+99))
    samtools faidx "$fa" "$chr:$lo-$hi" | sed "1s/>.*/>$sp/" >> "$OUT/upstream_Cexon.fa"
  else
    lo=$((start-99)); hi=$((start+1500))
    samtools faidx "$fa" "$chr:$lo-$hi" | sed "1s/>.*/>$sp/" | \
      python3 -c 'import sys; from Bio import SeqIO; r=next(SeqIO.parse(sys.stdin,"fasta")); r.seq=r.seq.reverse_complement(); SeqIO.write([r],sys.stdout,"fasta")' >> "$OUT/upstream_Cexon.fa"
  fi
  blo=$((gs-100000)); [ $blo -lt 1 ] && blo=1; bhi=$((ge+100000))
  echo -e "$sp\t$fa\t$chr:$blo-$bhi" >> "$OUT/locus_100kb_regions.txt"
done < "$OUT/species.tsv"
"$MAFFT" --genafpair --maxiterate 1000 --thread 4 "$OUT/upstream_Cexon.fa" > "$OUT/upstream_Cexon.mafft.fa" 2> "$OUT/mafft.log"

# ---- 3. JASPAR motif file + per-species Markov background + FIMO ---------------------------
mkdir -p "$OUT/jaspar"; unzip -o -q "$JZIP" -d "$OUT/jaspar"
python3 - "$OUT/jaspar" "$OUT/jaspar_insects.meme" <<'PY'
import sys, glob
out = open(sys.argv[2], 'w'); head = None
for f in sorted(glob.glob(sys.argv[1] + '/*.meme')):
    lines = open(f).read().splitlines()
    if head is None:
        head = lines[:lines.index(next(l for l in lines if l.startswith('MOTIF')))]
        out.write('\n'.join(head) + '\n\n')
    i = lines.index(next(l for l in lines if l.startswith('MOTIF')))
    out.write('\n'.join(lines[i:]) + '\n\n')
out.close()
print('motifs:', open(sys.argv[2]).read().count('MOTIF'))
PY
while IFS=$'\t' read -r sp fa reg; do
  samtools faidx "$fa" "$reg" > "$OUT/locus_100kb_$sp.fa"
  $MEME/fasta-get-markov -m 1 "$OUT/locus_100kb_$sp.fa" "$OUT/bg_$sp.txt"
done < "$OUT/locus_100kb_regions.txt"
while IFS=$'\t' read -r sp fa chr strand start gs ge; do
  grep -A10000 ">$sp" "$OUT/upstream_Cexon.fa" > /tmp/up_$sp.fa
  for T in 1e-4 1e-3; do
    $MEME/fimo --oc "$OUT/fimo/${sp}_$T" --thresh $T --bgfile "$OUT/bg_$sp.txt" \
      "$OUT/jaspar_insects.meme" /tmp/up_$sp.fa > /dev/null 2>&1
  done
done < "$OUT/species.tsv"

# ---- 4-7. column mapping, classification, expression join, tables ---------------------------
python3 shell/tfbs_classify.py "$OUT" "$E_RG" "$CT_RG" /tmp
# snRNA: % of Gr28b+ ovipositor nuclei expressing each TF, then merge into the table
/path/to/conda/envs/pyscRNA/bin/Rscript shell/tfbs_snrna.R "$OUT" \
  > "$OUT/snrna_detection.log" 2>&1 || echo "snRNA detection step failed (see snrna_detection.log)"
python3 shell/tfbs_classify.py "$OUT" "$E_RG" "$CT_RG" /tmp
