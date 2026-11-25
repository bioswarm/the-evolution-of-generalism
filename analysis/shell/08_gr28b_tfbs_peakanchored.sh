#!/bin/bash
# Peak-anchored re-analysis of TFBS conservation at the Gr28b locus.
# Replaces the exon-anchored alignment of shell/07_gr28b_tfbs.sh: homologs of the
# 10 distinct relaxed B. dorsalis MO peaks (Gr28b +/-5 kb) and the control window
# are found by blastn inside each species' Gr28b locus (50 kb beyond the outermost
# alternative first exon -> 30 kb beyond the last shared exon, transcription
# orientation; the 3' margin was widened from 10 kb to 30 kb after peak_108198280's
# true B. correcta homolog turned out to lie ~13 kb beyond the old window edge).
# Usage: shell/08_gr28b_tfbs_peakanchored.sh   (run from reproduce/)
# Outputs -> results/locus/tfbs/peak_anchored/
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=results/locus/tfbs/peak_anchored; mkdir -p "$OUT" "$OUT/fimo" "$OUT/blast"
CHIP=/path/to/conda/envs/chipseq/bin/python
MEME=/path/to/conda/envs/anvio-9/bin
MAFFT=/path/to/conda/envs/genera/bin/mafft
BLAST=/opt/ncbi-blast-2.10.0+/bin
GS=/path/to/data/genesets
HIC=/path/to/data/Hic
TDIR=results/locus/tfbs

# species: name fasta locus_chr locus_lo locus_hi gene_strand N1_exon(transcription-oriented pos)
# locus = 50 kb beyond outermost first exon -> 30 kb beyond last shared exon (gene end)
cat > "$OUT/species_locus.tsv" <<EOF2
Bdor	$GS/Bdor.fa	NC_064303.1	108119689	108224489	+	108181492
Bcor	$GS/Bcor/genome_aln/final/Bcor_final.fa	Bcor_chr1	8365207	8469583	-	8407930
Bmin	$GS/Bmin.fa	BmChr1	9427063	9538221	+	9498362
Btsu	$GS/Btsu/genome/ncbi/Btsu.fa	chr1	53368861	53479779	-	53408606
Zcuc	$GS/Bcuc.fa	NC_071668.1	10528893	10630247	+	10588317
Ztau	/path/to/data/genesets/Btau/RNAseq/RSEM/Btau.fa	chr1	68038077	68163536	-	68077885
EOF2

# ---- 1. query set: distinct MO peak intervals in Gr28b +/-5 kb, each +/-100 bp, + control ----
python3 - "$HIC/Bdor/results_hicexplorer/Bd_MO.pval0.01.300K.narrowPeak" "$OUT/queries.bed" <<'PY'
import sys
seen = {}
for l in open(sys.argv[1]):
    v = l.split('\t')
    if v[0] != 'NC_064303.1': continue
    s, e = int(v[1]), int(v[2])
    if s >= 108164689 and e <= 108199489:
        seen.setdefault((s, e), []).append((v[3], float(v[8])))
out = open(sys.argv[2], 'w')
for (s, e), ids in sorted(seen.items()):
    best = max(ids, key=lambda x: x[1])
    out.write(f'NC_064303.1\t{s}\t{e}\tpeak_{s}\t{best[0]}\t{best[1]:.3f}\n')
# control window (no peak), same +/-100 flank convention
out.write('NC_064303.1\t108180599\t108181026\tcontrol_108180599\t.\t.\n')
out.close()
print('queries:', len(seen) + 1)
PY
python3 - "$OUT/queries.bed" "$OUT/queries_flank100.bed" <<'PY'
i, o = open(__import__('sys').argv[1]), open(__import__('sys').argv[2], 'w')
for l in i:
    v = l.rstrip('\n').split('\t')
    o.write(f'{v[0]}\t{int(v[1])-100}\t{int(v[2])+100}\t{v[3]}\t{v[4]}\t{v[5]}\t{v[1]}\t{v[2]}\n')
PY
samtools faidx "$GS/Bdor.fa" NC_064303.1 > /dev/null 2>&1 || true
bedtools getfasta -fi "$GS/Bdor.fa" -bed "$OUT/queries_flank100.bed" -name+ > "$OUT/queries_flank100.fa" 2>/dev/null || \
  python3 - "$GS/Bdor.fa" "$OUT/queries_flank100.bed" "$OUT/queries_flank100.fa" <<'PY'
import sys
from Bio import SeqIO
g = SeqIO.to_dict(SeqIO.parse(sys.argv[1], 'fasta'))
o = open(sys.argv[3], 'w')
for l in open(sys.argv[2]):
    v = l.rstrip('\n').split('\t')
    s = g[v[0]].seq[int(v[1]):int(v[2])]
    o.write(f'>{v[3]}\n{s}\n')
PY

# ---- 2. locus target sequences + blastn ------------------------------------------------------
while IFS=$'\t' read -r sp fa chr lo hi strand n1; do
  [ -f "$fa.fai" ] || samtools faidx "$fa"
  samtools faidx "$fa" "$chr:$lo-$hi" > "$OUT/locus_$sp.fa"
done < "$OUT/species_locus.tsv"

while IFS=$'\t' read -r sp fa chr lo hi strand n1; do
  [ "$sp" = Bdor ] && continue
  $BLAST/blastn -query "$OUT/queries_flank100.fa" -subject "$OUT/locus_$sp.fa" \
    -task blastn -word_size 7 -evalue 1e-2 -dust no -soft_masking false \
    -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore' \
    > "$OUT/blast/locus_$sp.tsv" || true
done < "$OUT/species_locus.tsv"

# ---- 3. chained full-length identity (peak intervals only) + locus background ---
# (runs first: the TFBS status rules below need the chained-HSP query intervals)
python3 shell/tfbs_fulllength.py "$OUT"

# ---- 4. classify hits, extract E/control homologs, MAFFT, Bmin ATAC, FIMO ---
python3 shell/tfbs_peakanchor.py "$OUT"

python3 shell/tfbs_section_peakanchored.py "$OUT"

# ---- 5. full MSAs for SuppFig_E_MSA (E + partially conserved peaks) ----------
python3 shell/tfbs_msa.py "$OUT"
