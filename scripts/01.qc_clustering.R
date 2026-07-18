# Import the 10K PBMC data
counts <- Read10X_h5("~/Documents/HLH Disease Modelling Project 2026/SC3_v3_NextGem_SI_PBMC_10K_filtered_feature_bc_matrix.h5")

# Create Seurat object
seurat <- CreateSeuratObject(counts = counts, # GEX matrix of unnormalised data with cells as columns and features as rows
                             project = "pbmc10k",
                             min.cells = 3, # Keep features expressed in >=3 cells
                             min.features = 200) # Keep cells expressing >= 200 cells

# Quality control
# Adds a per-cell QC metric to the Seurat object's metadata as a new column "percent.mt".
# Finds all features with the "MT-" prefix and sums their counts per cell, 
# divides by that cell's total counts, and multiplies by 100.
# Results: for each cell, the percentage of its transcripts coming from the mitochondrial genome.
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

# Plot the distribution of three QC metrics.
# nFeature_RNA: number of features (genes) expressed per cell
# nCount_RNA: number of individual transcripts (UMIs) per cell
# percent.mt: percentage of a cell's transcripts mapping to mitochondrial genes
# ncol = 3 arranges the three panels side by side in one row.
# Cells with very low nFeature_RNA (empty/poor droplets), very high counts (potential doublets),
# or high percent.mt (dying/lysed cells).
qc_vln <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("qc_violin.png", plot = qc_vln, width = 12, height = 5, dpi = 300)

# Plot the relationships between different QC metrics.
# plot1: nCount_RNA (x) vs percent.mt (y): checks whether mito % relates to sequencing depth.
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")

# plot2: nCount_RNA (x) vs nFeature_RNA (y): counts vs genes detected per cell.
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
feature_scatter <- plot1 + plot2
ggsave("feature_scatter.png", plot = feature_scatter, width = 12, height = 5, dpi = 300)

# to-do: try for doubletfinder methods
seurat <- subset(seurat, subset = nFeature_RNA > 200 & 
                   nFeature_RNA < 6000 & 
                   percent.mt < 10)

# Log-normalisation: divide each cell's gene counts by that cell's total counts,
# multiply by a scale factor (default 10,000), then natural-log transform.
# Corrects for differences in sequencing depth between cells and puts values on a log scale.
seurat <- NormalizeData(seurat)

# Identify the 2000 most variable genes (high cell-to-cell variance after mean-variance
# adjustment). 
# Only use these for downstream analysis to reduce noise and computation.
seurat <- FindVariableFeatures(seurat, nfeatures = 2000)

# Center and scale the variable features to mean = 0, variance = 1 (z-score per gene).
# Stops highly expressed genes from dominating PCA purely because of larger values.
seurat <- ScaleData(seurat)

# Principal component analysis on the scaled variable features.
# Finds orthogonal axes (PCs) that are linear combinations of genes, ordered by how much
# cell-to-cell variance they capture.
# Each cell gets a score on each PC; each gene a loading.
# The PCs are eigenvectors of the gene covariance matrix (computed via
# truncated SVD for the top 50 PCs). This compresses ~2000 genes into a few informative
# PCs, reducing noise, and speeding up downstream computation.
# Caveat: it's linear (misses nonlinear structure) and variance-driven, so uncorrected
# technical variance (batch, depth) can dominate a PC.
seurat <- RunPCA(seurat)

# Plots the variance captured by each PC in descending order.
# Used to decide the optimal number of PCs for downstream steps.
# Find the elbow: the point where the curve flattens and difference in variance captured plateaus.
# Keep the PCs before the elbow (often ~10-20).
elbow <- ElbowPlot(seurat)                     
ggsave("elbow_plot.png", plot = elbow, width = 6, height = 4, dpi = 300)

# Build a shared nearest-neighbor (SNN) graph from the first 15 PCs.
# For each cell, finds its k nearest neighbors in PC space, then weights edges by how many
# neighbors each pair shares. 
# dims = 1:x should match the number of PCs chosen from ElbowPlot.
seurat <- FindNeighbors(seurat, dims = 1:15)

# Cluster cells by community detection such as louvain on that SNN graph.
# Resolution controls granularity: higher = more, granular clusters. 
# The last resolution run becomes the active identity (Idents). 
# Lets you compare clustering granularity without recomputing.
seurat <- FindClusters(seurat, resolution = c(0.3, 0.5, 0.8, 1.0))

# Compute a 2D UMAP embedding from the same 15 PCs for visualization.
# Nonlinear projection that places similar cells near each other in a 2D space, capturing local structure.
# Use the same dims as FindNeighbors so the graph and the plot are consistent.
seurat <- RunUMAP(seurat, dims = 1:15)

# to-do: save plot
clustree <- clustree(seurat) # Justify which resolution to use
ggsave("clustree.png", plot = clustree, width = 8, height = 10, dpi = 300)

Idents(seurat) <- "RNA_snn_res.0.3"  # Set resolution
# check idents of seurat samples and correct for batch effects



