#!/usr/bin/env bash
# Reproduce all R-generated publication figures into OUT_DIR (default: outputs/figures)
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export REPO_ROOT
export OUT_DIR="${OUT_DIR:-${REPO_ROOT}/outputs/figures}"
RSCRIPT="${RSCRIPT:-Rscript}"

mkdir -p "${OUT_DIR}"

SCRIPTS=(
  "scripts/figures/plot_sf13_gr28b_isoform_qpcr.R"
  "scripts/figures/plot_sf12_gr28b_isoform_rnaseq.R"
  "scripts/figures/plot_fig1f_and_sf5.R"
  "scripts/figures/plot_sf02_global_jsd.R"
  "scripts/figures/plot_sf03_per_tissue_pca.R"
  "scripts/figures/plot_sf04_per_tissue_jsd.R"
  "scripts/figures/plot_sf06_and_sf20.R"
  "scripts/figures/plot_sf07_to_sf11.R"
  "scripts/figures/plot_fig2.R"
  "scripts/figures/plot_fig5.R"
)

echo "=== Running ${#SCRIPTS[@]} figure-generation scripts -> ${OUT_DIR} ==="
for s in "${SCRIPTS[@]}"; do
  echo ""
  echo ">>> [${s}]"
  "${RSCRIPT}" "${REPO_ROOT}/${s}"
done

echo ""
echo "=== All ${#SCRIPTS[@]} figure scripts completed successfully ==="
ls -lh "${OUT_DIR}"/*.pdf
