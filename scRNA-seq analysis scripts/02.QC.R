# ============================================================
# 2. Quality control
# ============================================================

# Mitochondrial percentage per cell.
# PercentageFeatureSet finds all genes whose name starts with "MT-" (the
# mitochondrial genome), sums their counts in each cell, divides by that cell's
# total counts, and multiplies by 100.
# A high value means most of the cell's RNA is mitochondrial, which happens when
# a cell is dying: the membrane ruptures, cytoplasmic RNA leaks out, and the
# mitochondrial RNA is left behind.
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Distributions of the three standard QC metrics.
#   nFeature_RNA = how many distinct genes were detected in that cell
#   nCount_RNA = how many transcripts (UMIs) total in that cell
#   percent.mt = the mitochondrial fraction computed above
# Interpretation: very low nFeature = empty or failed droplet; very high
# nCount = possible two cells in one droplet (doublet); high percent.mt = dying.
p_qc_vln <- VlnPlot(pbmc, c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    ncol = 3, pt.size = 0.05)

# Plot the relationships between different QC metrics.
# plot1: nCount_RNA (x) vs percent.mt (y): checks whether mito % relates to sequencing depth.
p1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")

# plot2: nCount_RNA (x) vs nFeature_RNA (y): counts vs genes detected per cell.
p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

metadata <- pbmc@meta.data

round(sapply(metadata[, c("nFeature_RNA", "nCount_RNA", "percent.mt")],
             quantile, probs = c(0.25, 0.5, 0.75)), 1)

# subset pbmcs based on the violin plots
pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 & # drop empty or failed droplets
                 nFeature_RNA < 4000 & # drop probable doublets
                 percent.mt < 10) # drop dying cells

# Run scDblFinder as sanity check
sce <- as.SingleCellExperiment(pbmc)
sce$capture <- paste(pbmc$batch, pbmc$tenx_lane, sep = "_")

set.seed(1)
sce <- scDblFinder(sce, samples = "capture",
                   BPPARAM = BiocParallel::MulticoreParam(4))

pbmc$scDblFinder.class <- sce$scDblFinder.class 
table(pbmc$scDblFinder.class)

rm(sce)

# ============================================================
# 3. Attach ADT
# ============================================================

# Inspect the feature names
adt_counts <- counts(altExp(kotliarov, "ADT"))
rownames(adt_counts)

# keep only cells that survived RNA QC
cells <- colnames(pbmc)
adt_counts <- adt_counts[, cells]
pbmc[["ADT"]] <- CreateAssayObject(counts = adt_counts)

# Normalise the protein data with CLR (centred log-ratio)

# Why not the RNA method: protein counts here are compositional. The total
# antibody reads per cell is set by how deeply that cell was sequenced, which is
# arbitrary, so only the *ratios* between antibodies carry real information.
#
# What CLR does: divides each value by the geometric mean of its group, then
# takes the log. Because log(a/b) = log(a) - log(b), ratios become differences,
# and the arbitrary total cancels out.
pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR", margin = 2)

# ============================================================
# 4. RNA preprocessing
# ============================================================
DefaultAssay(pbmc) <- "RNA"

# Log-normalisation: divide each cell's gene counts by that cell's total counts,
# multiply by a scale factor (default 10,000), then natural-log transform.
# Corrects for differences in sequencing depth between cells and puts values on a log scale.
pbmc <- NormalizeData(pbmc)

# Identify the 2000 most variable genes (high cell-to-cell variance after mean-variance adjustment). 
# Only use these for downstream analysis to reduce noise and computation.
pbmc <- FindVariableFeatures(pbmc, nfeatures = 2000)

# Center and scale the variable features to mean = 0, variance = 1 (z-score per gene).
# Stops highly expressed genes from dominating PCA purely because of larger values.
pbmc <- ScaleData(pbmc)

# Principal component analysis on the scaled variable features.
# Finds orthogonal axes (PCs) that are linear combinations of genes, ordered by how much
# cell-to-cell variance they capture.
# Each cell gets a score on each PC; each gene a loading.
# The PCs are eigenvectors of the gene covariance matrix (computed via
# truncated SVD for the top 50 PCs). This compresses ~2000 genes into a few informative
# PCs, reducing noise, and speeding up downstream computation.
# Caveat: it's linear (misses nonlinear structure) and variance-driven, so uncorrected
# technical variance (batch, depth) can dominate a PC.
pbmc <- RunPCA(pbmc, npcs = 50, seed.use = 42)
p3 <- DimPlot(pbmc, reduction = "pca") + 
  NoLegend() 

# Plots the variance captured by each PC in descending order.
# Used to decide the optimal number of PCs for downstream steps.
# Find the elbow: the point where the curve flattens and difference in variance captured plateaus.
# Keep the PCs before the elbow (often ~10-20).
p4 <- ElbowPlot(pbmc, ndims = 50)

fig1 <- wrap_elements(p_qc_vln) / (p1 + p2) / (p3 + p4) +
  plot_annotation(tag_levels = "A")

ggsave("Figure_01_qc.png", fig1, width = 10, height = 12, dpi = 300)