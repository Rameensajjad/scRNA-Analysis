# Single-Cell RNA-seq Analysis Pipeline

> A complete end-to-end scRNA-seq workflow combining **STARsolo on Galaxy** for raw data
> preprocessing and **Scanpy (Python)** for downstream clustering and annotation of human PBMC data.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Part 1 — Raw Preprocessing (Galaxy)](#part-1--raw-preprocessing-galaxy)
- [Part 2 — Downstream Analysis (Python)](#part-2--downstream-analysis-python)
- [Key Findings](#key-findings)
- [Limitations](#limitations)
- [References](#references)

---

## Overview

Single-cell RNA sequencing (scRNA-seq) captures gene expression at the resolution of individual
cells, revealing rare populations and cellular heterogeneity invisible to bulk sequencing.

This repository implements a **two-part pipeline**:

| Part | Platform | Input | Output |
|---|---|---|---|
| **1 — Preprocessing** | STARsolo on Galaxy | Raw FASTQ files | Filtered gene expression matrix |
| **2 — Downstream Analysis** | Python / Scanpy | Filtered `.h5` matrix | Annotated clusters + figures |

**Sample:** Human PBMCs (peripheral blood mononuclear cells), 10x Genomics Chromium V3
chemistry — a standard reference dataset widely used in scRNA-seq benchmarking.

---

## Repository Structure

```
Preprocessing-of-scRNA/
│
├── README.md                    ← This file
├── SCRNA.ipynb                  ← Full Python analysis notebook
│
├── figures/                     ← 17 saved plots from Python analysis
├── results/
│   ├── scrna_results.h5ad       ← Final annotated AnnData object
│   └── cluster_annotations.csv ← Cluster to cell type mapping
│
└── galaxy_figures/              ← QC and alignment plots from Galaxy
```

---

## Part 1 — Raw Preprocessing (Galaxy)

### What This Step Does

Raw sequencing data from a 10x Chromium experiment consists of two paired FASTQ files per lane:
- **R1** — contains the cell barcode (16 bp) and UMI (10 bp)
- **R2** — contains the actual cDNA read (gene expression)

These cannot be used directly for analysis. They must be aligned to a reference genome, barcodes
must be validated, duplicates collapsed, and low-quality cells removed. This is done using
**STARsolo** on the Galaxy platform — a faster, open-source replacement for 10x Cell Ranger.

---

### Dataset

| Attribute | Value |
|---|---|
| Sample | PBMCs — 1k Healthy Donor (10x Genomics, V3) |
| Input | 2 lanes × R1 + R2 FASTQ files |
| Reference genome | hg19 / GRCh37 |
| Gene annotation | `Homo_sapiens.GRCh37.75.gtf` |
| Barcode whitelist | `3M-february-2018.txt.gz` (~3.7M valid barcodes) |
| Platform | usegalaxy.org |

---

### Preprocessing Steps

**1. Quality Assessment — FastQC / MultiQC**

Raw reads are assessed for per-base quality scores, GC content, adapter contamination, and
duplication levels. MultiQC aggregates logs from all lanes into a single HTML report.
For scRNA-seq, >70% uniquely mapped reads is the acceptable threshold.

**2. Alignment and UMI Counting — STARsolo**

STARsolo processes both FASTQ files simultaneously:
- Extracts barcodes and UMIs from R1
- Validates barcodes against the 10x whitelist (1 mismatch allowed)
- Aligns cDNA reads (R2) to hg19
- Assigns reads to genes and collapses PCR duplicates using UMI counting

Output is a **MEX-format matrix** (`matrix.mtx`, `barcodes.tsv`, `features.tsv`) compatible
with both Seurat and Scanpy.

**3. Cell Calling and QC — DropletUtils**

Not all detected barcodes represent real cells. DropletUtils uses a **knee plot** — barcodes
ranked by UMI count (log scale) — to separate real cells (above the inflection point) from
empty droplets (below). Additional filters remove damaged cells using three metrics:

| Metric | Interpretation |
|---|---|
| Low genes per cell | Likely empty droplet |
| High genes per cell | Likely doublet (two cells captured) |
| High % mitochondrial reads | Likely damaged or dying cell |

---

### Galaxy Results Summary

| Metric | Value |
|---|---|
| Total barcodes detected | 5,200 |
| Cells passing QC (knee plot) | 272 |
| Uniquely mapped reads | >70% |
| Output format | MEX (matrix.mtx + barcodes.tsv + features.tsv) |

The filtered matrix from this step feeds directly into Part 2.

---

## Part 2 — Downstream Analysis (Python)

### What This Step Does

The clean gene expression matrix is loaded into Python and taken through normalization,
dimensionality reduction, clustering, and cell type annotation using the **Scanpy** ecosystem.
Two samples (`s1d1` and `s1d3`) from the NeurIPS 2021 dataset are analyzed together.

---

### Environment

```python
!pip install anndata pooch scanpy scrublet leidenalg igraph celltypist decoupler
```

---

### Step 1 — Quality Control

Three biologically meaningful gene categories are flagged and used to assess cell quality:

```python
adata.var["mt"]   = adata.var_names.str.startswith("MT-")           # mitochondrial
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))  # ribosomal
adata.var["hb"]   = adata.var_names.str.contains(r"^HB[^EP]", regex=True)  # hemoglobin

sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)
```

Cells with too few genes (empty droplets) or too high mitochondrial content (damaged cells)
are removed. Violin and scatter plots visualize distributions before filtering.

---

### Step 2 — Doublet Detection

Doublets are droplets that accidentally captured two cells and appear as artificially
high-count cells. Scrublet scores each cell per batch:

```python
sc.pp.scrublet(adata, batch_key="sample")
adata = adata[~adata.obs["predicted_doublet"].to_numpy().astype(bool)].copy()
```

Predicted doublets are confirmed on UMAP before removal.

---

### Step 3 — Normalization

Raw counts are preserved, then the data is normalized and log-transformed:

```python
adata.layers["counts"] = adata.X.copy()  # save raw counts
sc.pp.normalize_total(adata)              # normalize to median library size
sc.pp.log1p(adata)                        # log1p transform
```

This ensures comparability across cells with different sequencing depths.

---

### Step 4 — Feature Selection and Dimensionality Reduction

2,000 highly variable genes (HVGs) are selected in a batch-aware manner, then PCA
reduces the data to 50 components before UMAP embedding:

```python
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
```

The PCA variance ratio plot guides the choice of how many components to use.

---

### Step 5 — Leiden Clustering

Graph-based Leiden clustering is applied at three resolutions to explore granularity:

```python
# flavor="igraph" is required in scanpy >= 1.10
sc.tl.leiden(adata, key_added="leiden_res0_02", resolution=0.02, flavor="igraph", n_iterations=2)
sc.tl.leiden(adata, key_added="leiden_res0_5",  resolution=0.5,  flavor="igraph", n_iterations=2)
sc.tl.leiden(adata, key_added="leiden_res2",    resolution=2,    flavor="igraph", n_iterations=2)
```

Low resolution (0.02) gives broad populations; high resolution (2.0) reveals fine substructure.

---

### Step 6 — Cell Type Annotation

Two complementary strategies are used and compared side-by-side:

**A) CellTypist** — model-based automated annotation trained on immune atlases:

```python
predictions = ct.annotate(adata, model="Immune_All_Low.pkl",
                          majority_voting=True,
                          over_clustering="leiden_res0_5")
adata = predictions.to_adata()
```

**B) Marker gene scoring** — manual scoring using known immune cell markers:

```python
for cell_type, genes in marker_dict.items():
    valid_genes = [g for g in genes if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=valid_genes, score_name=f"score_{cell_type}")
```

> **Note:** `decoupler 2.1.6` removed `run_mlm` / `run_ulm`. `sc.tl.score_genes` was used
> as a stable, equivalent alternative for cell type enrichment scoring.

Dotplots at multiple resolutions display marker gene expression per cluster,
and the two annotation results are overlaid on UMAP for comparison.

---

### Step 7 — Differential Expression

Wilcoxon rank-sum test identifies cluster-specific marker genes:

```python
sc.tl.rank_genes_groups(adata, groupby="leiden_res0_5", method="wilcoxon")
sc.tl.filter_rank_genes_groups(adata, min_fold_change=1.5)
```

Top 5 DEGs per cluster are visualized in a dotplot. Specific monocyte marker genes
(`LYZ`, `ACTB`, `S100A6`, `S100A4`, `CST3`) are explored on UMAP and violin plots.

---

## Key Findings

### Cell Populations Identified

| Cell Type | Key Markers | Notes |
|---|---|---|
| CD14+ Monocytes | LYZ, FCN1, CD14, S100A8 | Dominant myeloid population |
| CD16+ Monocytes | FCGR3A, TCF7L2, LYN | Non-classical monocytes |
| CD4+ T cells | CD4, IL7R, CCR7, LEF1 | Helper T cells |
| CD8+ T cells | CD8A, CD8B, GZMK, GZMA | Cytotoxic T cells |
| NK cells | GNLY, NKG7, KLRB1 | Natural killer cells |
| B cells | MS4A1, CD79A, PAX5 | Naive and memory B cells |
| Plasma cells | MZB1, JCHAIN, IGKC | Antibody-secreting cells |
| Dendritic cells | FCER1A, CLEC10A, CD1C | Antigen-presenting cells |
| Erythroblasts | HBA1, HBB, MKI67 | Erythroid precursors |

### Cluster 3 — Monocyte Signature

Cluster 3 showed the strongest monocyte identity with high expression of:

| Gene | Function |
|---|---|
| `LYZ` | Lysozyme — canonical monocyte marker |
| `S100A6` | Calcium-binding protein, inflammation |
| `S100A4` | Immune activation and migration |
| `CST3` | Cystatin C — myeloid lineage marker |
| `ACTB` | Beta-actin — cytoskeletal integrity |

### Annotation Comparison

Both CellTypist (automated) and marker gene scoring (manual) were applied independently.
Results showed strong agreement on major populations (T cells, B cells, Monocytes, NK cells),
with minor disagreements on rare or transitional populations — confirming the robustness
of the clustering at resolution 0.5.

---

## Limitations

| Issue | Impact |
|---|---|
| ~300 cells in Galaxy (subsampled) | Reduced statistical power vs full 1k dataset |
| chrX only on usegalaxy.org | Training mode — not genome-wide |
| No ambient RNA correction | SoupX / CellBender not applied |
| No batch correction | Harmony / scVI not used — samples integrated as-is |
| hg19 reference | hg38 provides better gene annotation coverage |
| decoupler 2.1.6 API breaking change | `run_mlm` / `run_ulm` removed — switched to `sc.tl.score_genes` |
| OmniPath / PanglaoDB server down | Used curated built-in marker dictionary as fallback |


---

## References

1. Bacon W. et al. *Pre-processing of 10X Single-Cell RNA Datasets.* Galaxy Training Network.
   https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing-tenx/tutorial.html

2. Wolf F.A., Angerer P., Theis F.J. (2018). SCANPY: large-scale single-cell gene expression data analysis.
   *Genome Biology*, 19, 15.

3. Dobin A. et al. (2013). STAR: ultrafast universal RNA-seq aligner.
   *Bioinformatics*, 29(1), 15–21.

4. Lun A. et al. (2019). EmptyDrops: distinguishing cells from empty droplets in droplet-based scRNA-seq.
   *Genome Biology*, 20(1), 63.

5. Dominguez Conde C. et al. (2022). Cross-tissue immune cell analysis reveals tissue-specific features in humans.
   *Science*, 376(6594).

6. Ewels P. et al. (2016). MultiQC: summarize analysis results for multiple tools in a single report.
   *Bioinformatics*, 32(19), 3047–3048.

7. Andrews S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data.
   https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
