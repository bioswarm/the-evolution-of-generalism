#!/bin/bash
# ATAC-seq (Methods): ENCODE ATAC-seq pipeline v2.2.3 run with caper (Singularity image atac-seq-pipeline_v2.2.x.sif).
# One library per species/condition; outputs used downstream (R/22, R/31):
#   call-macs2_signal_track/*/<sample>.fc.signal.bigwig
#   call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz
#   call-call_peak/*/<sample>.pval0.01.300K.narrowPeak.gz
# Usage: shell/04_atac_encode.sh <input.json>   (template: tools/atac_input_template.json)
set -euo pipefail
caper run atac.wdl -i "${1:?ENCODE input JSON}" --singularity
