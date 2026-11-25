#!/bin/bash
# Hi-C (Methods): Juicer -> juicer_tools pre (.hic) -> hicConvertFormat (.cool) -> HiCExplorer:
#   hicCorrectMatrix (KR/ICE), hicPCA (100 kb, PC1/PC2 bigWig), hicFindTADs (>= 5 kb; --correctForMultipleTesting fdr
#   --minDepth 3*bin --maxDepth 10*bin --step bin --thresholdComparisons 0.05 --delta 0.01), hicDetectLoops (<= 50 kb, maxLoopDistance 2 Mb).
# The R scripts start from results_hicexplorer/<res>/{hic_corrected.<res>.cool, TAD_<res>_domains.bed, loops.<res>.bedgraph}.
# Usage: shell/05_hic.sh <inter.hic> <out_dir> [threads]
set -euo pipefail
HIC=${1:?.hic file}; O=${2:?output dir}; T=${3:-20}; mkdir -p $O
for r in 5000 10000 25000 100000; do k=$((r/1000))k; mkdir -p $O/$k
  hicConvertFormat -m $HIC --inputFormat hic --outputFormat cool -o $O/$k/raw.$k.cool --resolutions $r
  hicCorrectMatrix correct -m $O/$k/raw.${k}_$r.cool --correctionMethod KR -o $O/$k/hic_corrected.$k.cool
  hicFindTADs -m $O/$k/hic_corrected.$k.cool --outPrefix $O/$k/TAD_$k --correctForMultipleTesting fdr \
    --minDepth $((3*r)) --maxDepth $((10*r)) --step $r --thresholdComparisons 0.05 --delta 0.01 -p $T
  [ $r -le 50000 ] && hicDetectLoops -m $O/$k/hic_corrected.$k.cool -o $O/$k/loops.$k.bedgraph --maxLoopDistance 2000000 -t $T
done
hicPCA -m $O/100k/hic_corrected.100k.cool --outputFileName $O/100k/pca1.100k.bw $O/100k/pca2.100k.bw --format bigwig
