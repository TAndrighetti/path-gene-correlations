# Path–Gene Correlation Analysis

This repository implements a **pathway–gene / pathway–pathway correlation framework**
using **microarray expression data** downloaded from GEO, combining **GSVA-based
enrichment scoring** with correlation analysis and visualization.

The project can be used in **two complementary ways**:
- an **interactive Jupyter notebook** for exploration and development
- a **fully reproducible command-line pipeline** (`.py` + `.sh`) for automated execution


## Scientific motivation

In many biological contexts, the key question is not whether a single gene is
differentially expressed, but whether:

- the expression of a **specific gene** is associated with the activity of a
  **biological pathway**, or
- two **biological pathways/programs** show coordinated behavior across samples.

This repository addresses this by:
- computing **gene set enrichment scores** using **GSVA**
- using **raw expression values** for single-gene analyses
- performing **Spearman correlation analysis**, optionally stratified by metadata
- generating **publication-ready scatter plots** with group-level statistics


## Core features

- Automatic download and parsing of **GEO microarray datasets**
- Probe-to-gene annotation via **Bioconductor (rpy2)**
- Flexible **metadata-based sample filtering**
- GSVA scoring for multi-gene pathways
- Direct expression usage for single-gene inputs
- Group-aware correlation analysis
- High-quality visualization with regression trend and confidence interval


## Repository structure

.
├── results/
│   └── correlation_IL10_neutrophils.png
│
├── codes/
|   ├── path_gene_correlation.ipynb
│   ├── run_path_gene_correlation.sh
│   └── path_gene_correlation.py
│
├── env/
│   ├── path_gene_correlation_env.yaml
│
└── README.md


## Usage options

### 1) Interactive notebook

Use `path_gene_correlation.ipynb` if you want to:
- inspect metadata and available biological groups
- iteratively adjust gene sets and filters
- explore correlations interactively

Typical workflow:
1. Download GEO dataset
2. Inspect metadata
3. Annotate probes
4. Define gene/pathway lists
5. Run correlation analysis and visualize results


### 2) Command-line pipeline

The script `path_gene_correlation.py` provides a **fully reproducible CLI interface**.

It supports:
- gene ↔ gene
- gene ↔ pathway
- pathway ↔ pathway correlations

Logic:
- If a gene set contains **more than one gene**, GSVA is applied
- If a gene set contains **a single gene**, raw expression values are used

The companion script `run_path_gene_correlation.sh` defines:
- GEO dataset accession
- gene/pathway definitions
- metadata filters (JSON)
- grouping variable
- color palette (JSON)
- output figure name

This design enables **portable and reproducible analyses** suitable for
batch execution and version control.


## Metadata filtering

Samples can be filtered using arbitrary metadata fields provided as a JSON string,
such as:
- disease status
- disease activity
- experimental condition

This allows biologically meaningful stratification and reproducible subgroup analyses.


## Visualization output

Generated scatter plots include:
- points colored by biological group
- a global regression line (OLS, visual aid)
- a 95% confidence interval
- a legend reporting **Spearman correlation (ρ), p-value, and sample size** per group


## Dependencies

Main dependencies include:
- pandas, numpy, scipy
- matplotlib, seaborn, statsmodels
- GEOparse
- decoupler
- rpy2
- Bioconductor annotation packages (platform-dependent)


It is intended as a **research-grade analytical scaffold**, not a general-purpose
package, and can be adapted to diverse datasets and biological questions.
