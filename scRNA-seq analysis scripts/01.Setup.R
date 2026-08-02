# =============================================================================
# Regulation of the degranulation programme in healthy cytotoxic effectors
#
# Dataset: Kotliarov et al. CITE-seq PBMC, 20 healthy donors, 2 batches.
# Aim: map how the primary HLH-associated genes are regulated in NK cells and
#      CD8 TEMRA, and build signalling networks that explain the transcription
#      factor activity of those effector states.
#
# Pipeline: QC -> annotation (RNA + protein) -> pseudobulk -> differential
#           expression -> TF activity -> prior knowledge network -> CARNIVAL
# =============================================================================


# ============================================================
# 0. Load packages
# ============================================================

library(scRNAseq)
library(ggplotify)
library(Seurat)
library(clustree)
library(SingleCellExperiment)
library(celldex)
library(SingleR)
library(dplyr)
library(pheatmap)
library(DESeq2)
library(decoupleR)
library(tidyr)
library(RColorBrewer)
library(tibble)
library(ggplot2)
library(OmnipathR)
library(gt)
library(leidenbase)
library(mclust)
library(patchwork)
library(ggrepel)
library(patchwork)
library(ashr)
library(ggdist)
library(msigdbr)
library(fgsea)
library(tidytext)
library(stringr)

# Fixes the random number generator so that anything stochastic (UMAP,
# clustering, downsampling) gives the same answer every time the script runs.
set.seed(42)

# ============================================================
# 1. Load data
# ============================================================

# A healthy pbmc scRNA-seq dataset of 20 samples from 2 experimental batches (10 samples each batch)
kotliarov <- KotliarovPBMCData(
  mode = c("rna", "adt"), # load both RNA and ADT modalities
  ensembl = FALSE, # keep gene symbols rather than ENSEMBL IDs
  location = TRUE, # Attach chromosome location metadata to genes
  legacy = FALSE) # Use current version of the object

# Convert to a seurat object
pbmc <- CreateSeuratObject(
  counts = counts(kotliarov), # pulls the raw count matrix
  meta.data = as.data.frame(colData(kotliarov)), # pulls the per-cell metadata
  project = "pbmc")

# Inspect metadata columns to find donor and batch columns needed for pseudobulking later
colnames(pbmc@meta.data)