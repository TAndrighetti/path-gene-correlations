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

---
---
## Gene Set Variation Analysis (GSVA)

**GSVA (Gene Set Variation Analysis)** is a non-supervised method that transforms a gene expression matrix (genes × samples) into a matrix of **pathway or gene set activity scores** (gene sets × samples). Instead of testing predefined group contrasts, GSVA estimates the **relative activity of a biological process in each individual sample**.

In contrast to classical enrichment approaches (e.g. GSEA), GSVA:
- operates at the **single-sample level**,
- is robust to inter-individual heterogeneity,
- enables continuous association analyses (e.g. gene ↔ pathway, pathway ↔ pathway),
- is well suited for exploratory transcriptomic analyses.

In this project, GSVA is used to summarize the coordinated expression of biologically related genes into a **single continuous score**, enabling direct correlation with individual gene expression or other biological programs.


## Defining gene sets to represent biological pathways

The interpretability of GSVA-based analyses critically depends on the **biological definition of the gene set**. A gene set should represent a **coherent functional program**, rather than an arbitrary collection of genes.

Recommended principles for gene set definition:

- **Clear biological concept**  
  The gene set should correspond to a well-defined process (e.g. neutrophil activation, inflammatory response, metabolic pathway).

- **Literature-driven curation**  
  Genes should be selected based on:
  - review articles,
  - functional studies,
  - previously published signatures,
  - curated databases (e.g. KEGG, Reactome, MSigDB, VFDB).

- **Functional coherence over size**  
  Smaller, well-curated gene sets often yield more interpretable GSVA scores than large, heterogeneous lists.

- **Avoid circularity**  
  Gene sets should not be defined solely based on differential expression within the same dataset used for downstream analysis, to prevent bias.

- **Context awareness**  
  The same biological process may require different gene sets depending on tissue, organism, or experimental condition.

In this repository, gene sets are treated as **operational representations of biological programs**, allowing the evaluation of whether the **global activity of a pathway** is associated with the expression of specific genes or other functional modules.



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
