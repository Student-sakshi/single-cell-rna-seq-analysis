# Human COVID-19 Immune scRNA-seq Analysis

Single-cell transcriptomic analysis of immune-cell populations from COVID-19 and healthy samples.

## Overview

This repository contains a Seurat-based analysis of human single-cell RNA sequencing (scRNA-seq) data from immune-cell samples from COVID-19 and healthy individuals.

The analysis covers:

- Quality control and preprocessing
- Doublet detection
- Normalization and highly variable gene selection
- Sample integration and batch correction
- PCA and UMAP dimensionality reduction
- Cell-type annotation
- Differential expression analysis
- GO Biological Process enrichment

The complete analysis and results are documented in the project report included in this repository.

## Dataset and Sample Metadata

The dataset contains four immune-cell samples:

| Sample | Condition |
|---|---|
| `covid_556` | COVID-19 |
| `covid_557` | COVID-19 |
| `covid_558` | COVID-19 |
| `HIP043` | Healthy |

The initial expression matrices contained:

| Sample | Cells | Genes |
|---|---:|---:|
| `covid_556` | 6,956 | 27,260 |
| `covid_557` | 15,102 | 29,509 |
| `covid_558` | 18,182 | 29,472 |
| `HIP043` | 11,943 | 29,205 |

The original expression matrices are not included in this repository.

The Seurat metadata included sample identity, total RNA counts, detected genes, donor, replicate, and sex.

## Analysis Workflow

The analysis followed this workflow:

**Expression Matrices → Quality Control → Doublet Detection → Normalization → Highly Variable Gene Selection → Batch Correction → PCA → UMAP → Clustering → Cell-Type Annotation → Cell-Type Proportions → Differential Expression → GO Pathway Enrichment**

## Quality Control and Preprocessing

Quality control was performed separately for each sample before merging.

Three QC metrics were used:

- `nFeature_RNA` — number of detected genes per cell
- `nCount_RNA` — total RNA counts per cell
- `percent.mt` — percentage of mitochondrial RNA

Cells with fewer than 200 detected genes were removed as low-quality cells or empty droplets.

Cells with more than 6,000 detected genes were removed as potential doublets or multiplets.

Cells with more than 10% mitochondrial RNA were removed because high mitochondrial content can indicate stressed or damaged cells.

### Doublet Detection

Doublet removal was performed using **DoubletFinder**. Doublets can occur when two cells are captured in the same droplet, producing an artificial expression profile containing signals from multiple cell types.

The optimal pK values identified for the four samples were:

| Sample | Optimal pK |
|---|---:|
| `covid_556` | 0.08 |
| `covid_557` | 0.09 |
| `covid_558` | 0.30 |
| `HIP043` | 0.30 |

## Normalization and Feature Selection

Normalization was performed using Seurat's **LogNormalize** method.

For each cell, gene counts were normalized by the total number of counts in that cell, multiplied by a scale factor, and log-transformed.

Highly variable genes (HVGs) were selected to focus downstream analyses on genes showing meaningful variation across cells rather than predominantly technical variation.

HVG selection was performed using Seurat's `FindVariableFeatures()` function with the **vst** method.

## Batch Correction and Sample Integration

The samples were compared before and after Seurat-based integration.

Before batch correction, the UMAP showed strong sample-associated separation, with cells from the four samples occupying distinct regions.

After integration, cells from different samples showed substantially improved mixing.

This indicates that sample-associated variation contributed strongly to the original embedding and that integration reduced this variation before downstream analysis.

## Dimensionality Reduction

Principal Component Analysis (PCA) was used to reduce the dimensionality of the gene-expression data before downstream visualization and clustering.

The PCA representation was used as the basis for UMAP visualization and clustering.

UMAP was subsequently used to visualize the low-dimensional representation of the integrated dataset.

## Clustering

Both **Louvain** and **Leiden** graph-based clustering were performed.

Leiden produced more refined clusters, separating some populations that appeared more broadly grouped with Louvain clustering.

The resulting clusters were used for downstream cell-type annotation.

## Cell-Type Annotation

Cell types were assigned using two complementary approaches.

### Automatic Annotation

**SingleR** was used for automatic cell-type annotation with the **Human Primary Cell Atlas** reference.

### Manual Annotation

Manual annotation was performed using cell-type marker gene expression.

The identified cell populations included:

- B cells
- T cells
- NK cells
- Plasma cells
- CD16 monocytes
- Plasmacytoid dendritic cells (pDC)
- Common lymphoid progenitors (CLP)
- Granulocyte-monocyte progenitors (GMP)
- Erythroblasts

## Cell-Type Proportions

Cell-type proportions were compared across the four samples.

HIP043 was dominated by NK cells, representing approximately 48% of the annotated cells.

The COVID-19 samples showed different cellular compositions. `covid_556` and `covid_557` had relatively large CLP populations, whereas `covid_558` was dominated by plasma cells, representing approximately 56% of the annotated population.

Smaller populations of T cells, B cells, erythroblasts, and pDCs were also observed.

## Differential Expression Analysis

Differential expression analysis was performed between selected immune-cell populations.

The main cell-type comparisons included:

- B cells vs T cells
- T cells vs CD16 monocytes

Differential expression was also performed between COVID-19 and healthy samples.

The analysis includes volcano-plot and top-DEG visualizations.

## Top Differentially Expressed Genes

The five most significant differentially expressed genes identified in the T-cell versus B-cell comparison were:

| Gene | log2 Fold Change | Adjusted p-value |
|---|---:|---:|
| `RP5-1028K7.2` | -7.04 | 8.83 × 10⁻⁷³ |
| `RASD2` | -11.98 | 1.55 × 10⁻⁷¹ |
| `CDK6` | -3.58 | 2.25 × 10⁻⁷¹ |
| `LAPTM4B` | -7.48 | 7.75 × 10⁻⁶⁶ |
| `SOX4` | -4.54 | 1.87 × 10⁻⁶⁵ |

## COVID-19 vs Healthy Differential Expression

Differential expression was performed between COVID-19 and healthy samples.

Genes showing increased expression in the COVID-19 group included:

- `IGHG4`
- `IGHG1`
- `IFI27`
- `S100A8`
- `S100A9`
- `LYZ`
- `VCAN`
- `IGHG3`
- `IGHG2`
- `IGHM`

Genes showing decreased expression included:

- `SH2D1B`
- `PRR20C`
- `SPON2`
- `KLRB1`
- `GNLY`
- `SLC38A7`
- `TRDC`
- `S1PR5`
- `MATK`
- `FCGR3A`

The observed differences included interferon-responsive, inflammatory myeloid, immunoglobulin-associated, and cytotoxic/NK/T-cell-associated genes.

## GO Biological Process Enrichment

Pathway analysis was performed for the T-cell versus B-cell comparison.

The pathway with the lowest p-value was **Cytoplasmic Translation (GO:0002181)**.

Cytoplasmic translation represents protein synthesis by cytoplasmic ribosomes and involves ribosomal components and translation-associated factors.

Other enriched biological processes included:

- Peptide Biosynthetic Process
- Macromolecule Biosynthetic Process
- Translation
- Gene Expression
- Ribonucleoprotein Complex Biogenesis
- Ribosomal Small Subunit Biogenesis
- Ribosome Biogenesis
- Ribosome Assembly
- rRNA Metabolic Process

The enrichment of translation- and ribosome-associated processes is consistent with increased protein-synthesis activity in the analyzed B-cell population.

## Key Findings

- Quality control and doublet detection were used to obtain high-quality single-cell profiles.
- Seurat integration reduced sample-associated variation and improved mixing between samples.
- PCA and UMAP provided low-dimensional representations of the integrated single-cell dataset.
- Distinct immune-cell populations were identified through marker-based cell-type annotation.
- The analysis identified B cells, T cells, NK cells, plasma cells, CD16 monocytes, pDCs, CLPs, GMPs, and erythroblasts.
- Cell-type composition differed substantially between samples, particularly for NK and plasma cells.
- Differential expression analysis identified distinct transcriptional programs between immune-cell populations.
- The COVID-19 versus healthy comparison showed differences in interferon-responsive, inflammatory, immunoglobulin-associated, and cytotoxic/NK/T-cell-associated genes.
- GO enrichment identified strong representation of translation- and ribosome-associated biological processes.

## Tools and Methods

- **Programming:** R
- **Single-cell analysis:** Seurat
- **Doublet detection:** DoubletFinder
- **Normalization:** LogNormalize
- **Feature selection:** Highly variable genes (`vst`)
- **Batch correction:** Seurat integration
- **Dimensionality reduction:** PCA, UMAP
- **Clustering:** Leiden, Louvain
- **Cell-type annotation:** SingleR, marker-based annotation
- **Reference:** Human Primary Cell Atlas
- **Differential expression:** Seurat
- **Pathway enrichment:** Enrichr / GO Biological Process

## Repository Contents

The repository contains:

- Analysis script(s)
- Project report
- Selected analysis figures
- Analysis results

The complete methodology and figures are described in the project report.

## Report

**Human COVID-19 Immune scRNA-seq Analysis — Project Report**

## Author

**Sakshi Parate**  
M.Sc. Bioinformatics, Saarland University
