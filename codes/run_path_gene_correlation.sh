#!/usr/bin/env bash
set -euo pipefail

# ================================
# Path–gene / path–path correlation
# ================================

# --- Inputs ---
GEO_DATASET="GSE59071"
OUT_FIG="correlation_IL10_neutrophils.png"

# --- Gene / pathway definitions ---
SET1="IL10"
SET2="FCGR3B,CEACAM8,CSF3R,CXCR1,CXCR2,ITGAM,MPO,ELANE,PRTN3,CTSG,LTF,MMP8,OLFM4,CD66B"


METADATA_FILTERS='{"disease": ["ulcerative colitis", "control"],"disease_activity": ["active", "inactive", "normal"]}'
GROUPBY="disease_activity"
PALETTE='{"active": "#1f77b4", "inactive": "#ff7f0e", "normal": "#2ca02c","all": "black"}'

NAME1="IL10 expression"
NAME2="Neutrophil Enrichment Score (GSVA)"

# --- Plot title ---
TITLE="IL10 expression vs Neutrophil enrichment (GSVA)"



# ================================
# Run analysis
# ================================

python path_gene_correlation.py \
  --geo "${GEO_DATASET}" \
  --filters "${METADATA_FILTERS}" \
  --groupby "${GROUPBY}" \
  --palette "${PALETTE}" \
  --set1 "${SET1}" \
  --set2 "${SET2}" \
  --name1 "${NAME1}" \
  --name2 "${NAME2}" \
  --title "${TITLE}" \
  --out "${OUT_FIG}"

echo "[INFO] Correlation analysis completed."
echo "[INFO] Output figure: ${OUT_FIG}"

