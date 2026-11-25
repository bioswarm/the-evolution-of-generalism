# Avoidance-based host filtering facilitates the evolution of generalism

Upstream bioinformatics pipelines and publication figure-generation scripts accompanying the manuscript **"Avoidance-based host filtering facilitates the evolution of generalism"**.

---

## Repository Structure

| Path | Description |
|---|---|
| `scripts/figures/` | Standalone R scripts reproducing computational main and supplementary figures |
| `run_all_figures.sh` | Batch runner executing all 10 figure scripts into `OUT_DIR` (default: `outputs/figures/`) |
| `analysis/shell/` | Upstream shell & Python/R pipelines (`kallisto`, ENCODE ATAC-seq v2.2.3, Juicer/HiCExplorer, IQ-TREE/FastTree, peak-anchored BLASTN + MAFFT + FIMO JASPAR 2026 TFBS scanning) |
| `analysis/R/00_common.R` | Shared sample sheet, color palettes, `ggplot2` themes, and helper functions |
| `analysis/config.yaml` | Input path configuration template |
| `analysis/renv.lock` | Pinned R 4.5.3 package versions (`renv::restore()`) |
| `analysis/tools/` | Supporting R helpers and external tool version manifest |

---

## Figure-to-Script Mapping

| Manuscript Figure | Script |
|---|---|
| **Fig. 1f** | `scripts/figures/plot_fig1f_and_sf5.R` |
| **Fig. 2a–e** | `scripts/figures/plot_fig2.R` |
| **Fig. 5a–d** | `scripts/figures/plot_fig5.R` |
| **Supplementary Fig. 2** | `scripts/figures/plot_sf02_global_jsd.R` |
| **Supplementary Fig. 3** | `scripts/figures/plot_sf03_per_tissue_pca.R` |
| **Supplementary Fig. 4** | `scripts/figures/plot_sf04_per_tissue_jsd.R` |
| **Supplementary Fig. 5** | `scripts/figures/plot_fig1f_and_sf5.R` |
| **Supplementary Fig. 6** | `scripts/figures/plot_sf06_and_sf20.R` |
| **Supplementary Figs. 7a–b, 8–11** | `scripts/figures/plot_sf07_to_sf11.R` |
| **Supplementary Fig. 12** | `scripts/figures/plot_sf12_gr28b_isoform_rnaseq.R` |
| **Supplementary Fig. 13** | `scripts/figures/plot_sf13_gr28b_isoform_qpcr.R` |
| **Supplementary Fig. 20** | `scripts/figures/plot_sf06_and_sf20.R` |

*Note:* Behavioral assays, HPLC chromatograms, SEM/RNA-FISH micrographs, electrophysiological traces, and schematics (Fig. 1b, 1d; Fig. 3; Fig. 4; Supplementary Figs. 1, 14–19) have their underlying quantitative measurements provided in the Source Data and Supplementary Data files accompanying the manuscript.

---

## Usage

### 1. Prepare Input Data
Raw sequencing reads and processed tables are not stored in this repository. Download:
- Raw sequencing reads from NCBI SRA (see **Data Availability** below)
- Published tabular inputs (`Source_Data.xlsx` and `Supplementary Data 1–10`) from the journal article page

Configure local input paths in `analysis/config.yaml` and set environment variables as needed:

```bash
export REVIEW_ROOT=/path/to/your/workspace
export REPRO_ROOT=/path/to/your/workspace/reproduce
export OUT_DIR=outputs/figures

./run_all_figures.sh
```

### 2. Software Environment

- **R**: `4.5.3` (exact package versions locked in `analysis/renv.lock`)
- **Key R packages**: `ggplot2`, `patchwork`, `data.table`, `dplyr`, `tidyr`, `ComplexHeatmap`, `pheatmap`, `philentropy`, `Seurat`, `celda` (`decontX`), `DoubletFinder`, `sleuth`, `mashr`, `ashr`, `ape`, `ggtree`, `rtracklayer`, `GenomicRanges`, `rhdf5`
- **External tools**: `kallisto v0.51.1`, `ENCODE ATAC-seq pipeline v2.2.3`, `Juicer` / `HiCExplorer`, `MAFFT v7.526`, `trimAl v1.4.rev15`, `IQ-TREE v1.6.12`, `FastTree v2.1.11`, `HyPhy (aBSREL)`, `MEME Suite v5.5.9 (FIMO)` + `JASPAR 2026 CORE insects` (full version manifest in `analysis/tools/external_tool_versions.txt`).

---

## Data Availability

- **RNA-seq (285 libraries across five *Bactrocera* species)**: NCBI SRA BioProject [`PRJNA1565472`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1565472) (`SRP758203`)
- **Mature-ovipositor ATAC-seq (*B. dorsalis* and *B. minax*)**: NCBI SRA BioProject [`PRJNA1562201`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1562201) (`SRP756271`; `SRR41771283–SRR41771284`)
- **Public *B. tsuneonis* RNA-seq re-analyzed**: NCBI SRA BioProject [`PRJNA1241247`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1241247) (`SRR33000917–SRR33000937`)
- **Reference genomes**: `GCF_023373825.1` (*B. dorsalis*), `GCF_028554725.1` (*B. cucurbitae*), `GCA_031772095.1` (*B. tau*), `GCA_027475135.1` (*B. correcta*), `GCA_021498325.1` (*B. minax*)

---

## License

Released under the [MIT License](LICENSE).
