# ============================================================
# 0. Load packages
# ============================================================

library(scRNAseq)
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

set.seed(42)
# ============================================================
# 1. Load data
# ============================================================

# A healthy pbmc scRNA-seq dataset of 20 samples from 2 experimental batches (10 samples each batch)
kotliarov <- KotliarovPBMCData(
  mode = c("rna", "adt"),
  ensembl = FALSE,
  location = TRUE,
  legacy = FALSE)

# Create seurat object 
pbmc <- CreateSeuratObject(
  counts = counts(kotliarov),
  meta.data = as.data.frame(colData(kotliarov)))

# What is the donor/batch column for pseudobulking? 
colnames(pbmc@meta.data)

# ============================================================
# 2. Quality control
# ============================================================

# Adds a per-cell QC metric to the Seurat object's metadata as a new column "percent.mt".
# Finds all features with the "MT-" prefix and sums their counts per cell, 
# divides by that cell's total counts, and multiplies by 100.
# Results: for each cell, the percentage of its transcripts coming from the mitochondrial genome.
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Plot the distribution of three QC metrics.
# nFeature_RNA: number of features (genes) expressed per cell
# nCount_RNA: number of individual transcripts (UMIs) per cell
# percent.mt: percentage of a cell's transcripts mapping to mitochondrial genes
# ncol = 3 arranges the three panels side by side in one row.
# Cells with very low nFeature_RNA (empty/poor droplets), very high counts (potential doublets),
# or high percent.mt (dying/lysed cells).
p_qc_vln <- VlnPlot(pbmc, c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05)
ggsave("01_qc_violin.png", p_qc_vln, width = 10, height = 4, dpi = 300)

# Plot the relationships between different QC metrics.
# plot1: nCount_RNA (x) vs percent.mt (y): checks whether mito % relates to sequencing depth.
p1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")

# plot2: nCount_RNA (x) vs nFeature_RNA (y): counts vs genes detected per cell.
p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("02_qc_scatter.png", p1 + p2, width = 10, height = 4, dpi = 300)

# subset pbmcs based on the violin plots
pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 &
                nFeature_RNA < 4000 &
                percent.mt < 10)

# ============================================================
# 3. Attach ADT
# ============================================================

# Inspect the feature names
adt_counts <- counts(altExp(kotliarov, "ADT"))
rownames(adt_counts)

cells <- colnames(pbmc)
adt_counts <- adt_counts[, cells]
pbmc[["ADT"]] <- CreateAssayObject(counts = adt_counts)

# CLR (centered log-ratio) transforms each value by dividing by the geometric mean of its group, 
# then taking the log. Protein counts from CITE-seq are compositional. 
# The total counts per cell are set by sequencing depth, so only ratios between features carry information.
# CLR is the standard transform for compositional data because it maps ratios to differences in log space, 
# making the result independent of the arbitrary total.
# RNA-style log normalisation divides by total counts per cell and scales to a fixed factor.
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
ggsave("03_pca.png", DimPlot(pbmc, reduction = "pca") + NoLegend(), width = 6, height = 5, dpi = 300)

# Plots the variance captured by each PC in descending order.
# Used to decide the optimal number of PCs for downstream steps.
# Find the elbow: the point where the curve flattens and difference in variance captured plateaus.
# Keep the PCs before the elbow (often ~10-20).
ggsave("04_elbow.png", ElbowPlot(pbmc, ndims = 50), width = 6, height = 4, dpi = 300)

# ============================================================
# 5. Clustering
# ============================================================

# Build a shared nearest-neighbor (SNN) graph from the first 15 PCs.
# For each cell, finds its k nearest neighbors in PC space, then weights edges by how many
# neighbors each pair shares. 
# dims = 1:x should match the number of PCs chosen from ElbowPlot.
pbmc <- FindNeighbors(pbmc, dims = 1:15)
res_seq <- seq(0.3, 1.0, by = 0.1)

# Leiden clustering
pbmc <- FindClusters(pbmc, algorithm = 4, n.iter = 10, resolution = res_seq, random.seed = 42)

leiden_cols <- grep("^RNA_snn_res", colnames(pbmc@meta.data), value = TRUE)

for (col in leiden_cols) {
  pbmc@meta.data[[sub("^RNA_snn", "leiden", col)]] <- pbmc@meta.data[[col]]
}

# Louvain — overwrites RNA_snn_res.* columns, hence the copy above
pbmc <- FindClusters(pbmc, algorithm = 1, n.start = 10, n.iter = 10,
                     resolution = res_seq, random.seed = 42)
for (col in leiden_cols) {
  pbmc@meta.data[[sub("^RNA_snn", "louvain", col)]] <- pbmc@meta.data[[col]]
}

# Compute a 2D UMAP embedding from the same 15 PCs for visualization.
# Nonlinear projection that places similar cells near each other in a 2D space, capturing local structure.
# Use the same dims as FindNeighbors so the graph and the plot are consistent.
pbmc <- RunUMAP(pbmc, dims = 1:15, seed.use = 42)

comparison <- data.frame(
  resolution = res_seq,
  n_leiden = sapply(res_seq, function(r) length(unique(pbmc@meta.data[[paste0("leiden_res.", r)]]))),
  n_louvain = sapply(res_seq, function(r) length(unique(pbmc@meta.data[[paste0("louvain_res.", r)]]))),
  ARI = sapply(res_seq, function(r) adjustedRandIndex(
    pbmc@meta.data[[paste0("leiden_res.", r)]],
    pbmc@meta.data[[paste0("louvain_res.", r)]])))

print(comparison)

ggsave("05_leiden_resolutions.png",
       DimPlot(pbmc, group.by = paste0("leiden_res.", res_seq),
               label = TRUE, ncol = 4) & NoLegend(),
       width = 16, height = 8, dpi = 300)

ggsave("05b_louvain_resolutions.png",
       DimPlot(pbmc, group.by = paste0("louvain_res.", res_seq),
               label = TRUE, ncol = 4) & NoLegend(),
       width = 16, height = 8, dpi = 300)

# Find how cell clusters separate with each resolution
ggsave("06_clustree.png", clustree(pbmc, prefix = "leiden_res."), width = 10, height = 10, dpi = 300)

Idents(pbmc) <- "leiden_res.0.7"
ggsave("07_umap_res07.png",
       DimPlot(pbmc, reduction = "umap", group.by = "leiden_res.0.7", label = TRUE) + 
         NoLegend() +
         ggtitle("Leiden, resolution 0.7"),
       width = 7, height = 6, dpi = 300)

# ============================================================
# 6. Reference-based annotation
# ============================================================

expression <- GetAssayData(pbmc, layer = "data")

reference <- celldex::MonacoImmuneData()

# SingleR labels each cell by finding which reference cell type its gene expression ranks correlate best with (Spearman), 
# then iteratively narrows down among the top-scoring candidates using markers specific to just those labels.
pred.main <- SingleR(
  test = expression, 
  ref = reference, 
  labels = reference$label.main)

pred.fine <- SingleR(
  test = expression, 
  ref = reference, 
  labels = reference$label.fine)

table(pred.main$labels, pred.fine$labels)

pbmc$main <- pred.main$labels
pbmc$fine <- pred.fine$labels

p_main <- DimPlot(pbmc, group.by = "main", label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("label.main")

p_fine <- DimPlot(pbmc, group.by = "fine", label = TRUE, repel = TRUE, label.size = 2.5) +
  NoLegend() + ggtitle("label.fine")

ggsave("08_singler_main_vs_fine.png", p_main + p_fine, width = 14, height = 6, dpi = 300)

pbmc$singler.main <- pred.main$pruned.labels
pbmc$singler.fine <- pred.fine$pruned.labels

# Score heatmap: per-cell correlation against each reference label.
# One clearly dominant row per cell = confident call; several bright rows = ambiguous,
# usually meaning the reference can't resolve those subtypes at this granularity.
png("09_singler_scoreheatmap.png", width = 10, height = 8, units = "in", res = 300)
plotScoreHeatmap(pred.main)
dev.off()

# Delta = assigned label's score minus the median score across other labels.
# Large delta = confident assignment; small delta = marginal, pruned to NA by pruneScores().
# A label whose entire distribution sits low is suspect (cell type may be absent here).
ggsave("10_singler_delta.png", plotDeltaDistribution(pred.main, ncol = 4),
       width = 12, height = 8, dpi = 300)
summary(is.na(pred.main$pruned.labels))

# ============================================================
# 7. Cross-check labels against clusters
# ============================================================
tab <- table(Assigned = pred.main$pruned.labels, Cluster = pbmc$leiden_res.0.7)

png("11_singler_cluster_crosstab.png", width = 8, height = 6, units = "in", res = 300)
pheatmap(prop.table(tab, margin = 2))
dev.off()

# ============================================================
# 8. Protein-based annotation refinement
# ============================================================
set.seed(42)

# Monaco's "T cells" label is not a distinct population, but a compartment of
# cells that could not be resolved by RNA. Resolve with protein instead.
Idents(pbmc) <- "leiden_res.0.7"
DefaultAssay(pbmc) <- "ADT"

# Lineage panel
# Anything CD3-negative is not a T cell regardless of CD4 level,
# since CD4 protein is expressed on monocytes.
# CD14 + CD16 together split classical / intermediate / non-classical monocytes.
adt_lineage <- c("CD3-PROT", "CD4-PROT", "CD8-PROT", "CD19-PROT",
                 "CD56-PROT", "CD14-PROT", "CD16-PROT", "HLA-DR-PROT")

ggsave("12_adt_lineage_dotplot.png",
       DotPlot(pbmc, features = adt_lineage) + RotatedAxis(),
       width = 8, height = 7, dpi = 300)

ggsave("13_adt_lineage_violin.png",
       VlnPlot(pbmc, adt_lineage, stack = TRUE, flip = TRUE, pt.size = 0) + NoLegend(),
       width = 9, height = 8, dpi = 300)


# Subset markers, once lineage is settled
adt_subset <- c("CD45RA-PROT", "CD45RO-PROT", "CD62L-PROT", "CD197-PROT",
                "CD127-PROT", "CD25-PROT", "CD57-PROT", "CD161-PROT",
                "TCRgd-PROT", "CD279-PROT", "CD183-PROT",
                "CD11c-PROT", "CD123-PROT", "CD1c-PROT", "CD303-PROT",
                "CD20-PROT", "CD27-PROT", "IgD-PROT", "IgM-PROT", "CD34-PROT")

ggsave("14_adt_subset_dotplot.png",
       DotPlot(pbmc, features = adt_subset) + RotatedAxis(),
       width = 12, height = 7, dpi = 300)

# Background
# Signal sitting at isotype level is true negative; above it is real,
# or ambient contamination.
isotypes <- grep("sotype", rownames(pbmc[["ADT"]]), value = TRUE)

ggsave("15_isotypes.png",
       VlnPlot(pbmc, isotypes, stack = TRUE, flip = TRUE, pt.size = 0) + NoLegend(),
       width = 7, height = 5, dpi = 300)

ggsave("16_adt_cd4cd8_umap.png",
       FeaturePlot(pbmc, c("CD3-PROT", "CD4-PROT", "CD8-PROT", "CD56-PROT"),
                   ncol = 2, min.cutoff = "q05", max.cutoff = "q95"),
       width = 10, height = 9, dpi = 300)

# RNA markers as secondary evidence
DefaultAssay(pbmc) <- "RNA"

markers_rna <- c(
  "CD3D", "CD3E", # pan-T
  "CD4", "CD8A", "CD8B", # T lineage split
  "IL7R", "CCR7", "SELL", "TCF7", # naive / central memory
  "FOXP3", "IL2RA", # Treg
  "GZMK", "GZMB", "NKG7", "PRF1", # cytotoxic
  "KLRD1", "NCAM1", "KLRB1", # NK, MAIT
  "TRDC", "TRGC1", # gamma-delta
  "MS4A1", "CD79A", "TCL1A", # B
  "JCHAIN", "MZB1", # plasmablast
  "LYZ", "CD14", "FCGR3A", "MS4A7", # monocytes
  "FCER1A", "CLEC9A", "LILRA4", # cDC, pDC
  "CD34", "SPINK2", "PRSS57", # progenitors
  "PPBP", "PF4") # platelets

markers_rna <- intersect(markers_rna, rownames(pbmc))

ggsave("17_rna_markers_dotplot.png",
       DotPlot(pbmc, features = markers_rna) + RotatedAxis(),
       width = 14, height = 7, dpi = 300)

# Double and cell quality check
# Clusters showing exclusive lineage markers (CD3+ and CD14+, or
# CD4-high and CD8-high) are doublet candidates. 
# Elevated nCount/nFeature relative to neighbouring clusters supports this.
ggsave("18_counts_features.png",
       VlnPlot(pbmc, c("nCount_RNA", "nFeature_RNA"), pt.size = 0, ncol = 1) + NoLegend(),
       width = 10, height = 7, dpi = 300)

ggsave("19_percent_mt.png",
       VlnPlot(pbmc, "percent.mt", pt.size = 0) + NoLegend(),
       width = 8, height = 5, dpi = 300)

# Unbiased markers
all_markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  max.cells.per.ident = 500,
  random.seed = 42)

write.csv(all_markers, "all_markers_leiden_res07.csv", row.names = FALSE)

# Filter for genes with adjusted p-value < 0.05.
# Rank by fold-change but require the gene to be expressed in at least 50% of the cluster.
top_markers <- all_markers %>%
  filter(p_val_adj < 0.05, pct.1 > 0.5) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10) %>%
  ungroup()

print(top_markers, n = Inf)

ggsave("20_top_markers_heatmap.png",
       DoHeatmap(subset(pbmc, downsample = 100),
                 features = top_markers$gene, size = 3) + NoLegend(),
       width = 14, height = 16, dpi = 300)

# ============================================================
# 9. Final labels
# ============================================================
singler_by_cluster <- table(pbmc$leiden_res.0.7, pbmc$singler.main)
majority <- colnames(singler_by_cluster)[max.col(singler_by_cluster)]
names(majority) <- rownames(singler_by_cluster)

# Labels assigned from three pieces of evidence: SingleR/Monaco,
# ADT surface phenotype, and unbiased RNA marker. 
# Where they disagreed, protein and canonical markers took precedence over the reference call.
final_label <- c(
  # CD3+ CD4+ | CD62L-hi CD197(CCR7)-hi CD45RA+ CD127+ | RNA: CCR7, SELL, TCF7, LDHB, PIK3IP1, NOSIP
  "1"  = "CD4 naive T",
  
  # CD3+ CD4+ | CD45RO-hi CD45RA-lo CD25-hi CD127+ | RNA: IL7R, LTB, IL32, CD69, ITGB1; CCR7/SELL absent
  "2"  = "CD4 memory T",
  
  # CD14-hi CD16-neg HLA-DR+ CD11c+ | RNA: CD14, LYZ, S100A8/9/12, VCAN, LGALS2, MS4A6A
  "3"  = "classical monocytes",
  
  # CD3+ CD8+ | CD57-hi CD62L-neg CD45RO+ | RNA: GZMH, GNLY, FGFBP2, NKG7, CST7, GZMA/GZMM
  "4"  = "CD8 TEMRA",
  
  # CD3+ CD8+ | CD161-hi CD45RO-hi CD127+ CD279+ | RNA: GZMK, KLRB1, DUSP2, IL7R, CCL5
  "5"  = "MAIT / CD8 EM",
  
  # CD3-neg | CD56-hi CD16-hi | RNA: KLRF1, KLRD1, SPON2, CLIC3, PRF1, GZMB, GNLY
  "6"  = "NK",
  
  # CD19+ CD20-hi | IgD-hi IgM-hi CD27-neg | RNA: TCL1A, MS4A1, CD79A/B, HLA-DQ/DR
  "7"  = "Naive B",
  
  # CD3+ CD8-hi | CD62L+ CD197-hi CD45RA+ | RNA: CD8B, CCR7, TCF7, SELL, NOSIP, PIK3IP1, LDHB
  "8"  = "CD8 naive T",
  
  # No lineage-defining protein; percent.mt spans full 0-10% range
  # RNA: only MT-* genes and MALAT1, several with pct.2 > pct.1 — QC artifact, not a cell type
  "9"  = "Low quality",
  
  # CD19+ CD20-hi | CD27+ IgM+ IgD-lo CD11c+ | RNA: BANK1, MS4A1, CD79A/B, IGJ; TCL1A absent
  "10" = "Memory B",
  
  # CD16-hi CD14-lo CD11c-hi | RNA: FCGR3A, MS4A7, CDKN1C, CSF1R, LST1, LILRB2, TCF7L2
  "11" = "non-classical monocytes",
  
  # CD1c-hi CD11c+ HLA-DR-hi CD14-neg | RNA: FCER1A, CLEC10A, CD1C, CPVL, HLA-DPA1/DQA1
  "12" = "cDC2",
  
  # CD14-hi CD16-neg CD11c+ (same lineage protein as cl.3)
  # RNA: S100A8 + FCGR1A(CD64), FOLR3, GBP1, WARS, TNFSF10, TYMP — IFN-stimulated state
  "13" = "activated classical monocytes",
  
  # CD303(CLEC4C)-hi CD123-hi CD11c-neg HLA-DR+ | RNA: LILRA4, CLEC4C, SCT, SERPINF1, DNASE1L3, LRRC26
  "14" = "pDC",
  
  # CD3 + CD4 + CD8 + CD19 + CD16 all positive in one cluster (mutually exclusive lineages)
  # RNA: MKI67, TYMS, RRM2, TK1, PCNA, STMN1 — cycling signature, but protein says doublet
  "15" = "Doublets",
  
  # All lineage proteins flat
  # RNA: SDPR(CAVIN2), PPBP, PF4, HIST1H2AC, TSC22D1 — unambiguous platelet
  "16" = "Platelets")

lab <- unname(final_label[as.character(pbmc$leiden_res.0.7)])
names(lab) <- colnames(pbmc)
pbmc$celltype <- lab

celltype_levels <- c(
  "CD4 naive T", "CD4 memory T",
  "CD8 naive T", "CD8 TEMRA", "MAIT / CD8 EM",
  "NK",
  "Naive B", "Memory B",
  "classical monocytes", "activated classical monocytes", "non-classical monocytes",
  "cDC2", "pDC",
  "Platelets", "Doublets", "Low quality")

pbmc$celltype <- factor(pbmc$celltype, levels = celltype_levels)
Idents(pbmc) <- "celltype"

ggsave("21_umap_celltype.png",
       DimPlot(pbmc, label = TRUE, repel = TRUE, label.size = 3) + NoLegend(),
       width = 8, height = 7, dpi = 300)

# --- Remove non-cell and failed populations ---
pbmc.clean <- subset(pbmc, subset = celltype %in%
                       c("Low quality", "Doublets", "Platelets"), invert = TRUE)
saveRDS(pbmc.clean, "pbmc_annotated.rds")

# ============================================================
# 10. Pseudobulk
# ============================================================

# Donor 209 has cells in both batches. 
# Kept as two separate samples.
pbmc$pb_sample <- paste(pbmc$sampleid, pbmc$batch, sep = "_")

# Cells per celltype-sample pair
cells_per_group <- table(pbmc$celltype, pbmc$pb_sample)
rowSums(cells_per_group >= 10)

# AggregateExpression() sums all counts across all cells sharing celltype and sample.
# Returns a matrix of features x groups.
# Summing keeps the raw counts for DESeq2.
pseudo <- AggregateExpression(
  pbmc,
  assays = "RNA",
  slot = "counts",    
  group.by = c("celltype", "pb_sample"),
  return.seurat = FALSE)$RNA

# Sample sheet
meta_pb <- data.frame(column = colnames(pseudo), stringsAsFactors = FALSE) %>%
  mutate(
    celltype = sub("_[0-9]+-[0-9]+$", "", column),
    pb_sample = sub("^.*_([0-9]+-[0-9]+)$", "\\1", column),
    sampleid = sub("-.*$", "", pb_sample),
    batch = sub("^.*-", "", pb_sample))

rownames(meta_pb) <- meta_pb$column
meta_pb$celltype <- factor(meta_pb$celltype)
meta_pb$batch <- factor(meta_pb$batch)

# Filter and normalise 
# design = ~ 1 because no comparison is being run
dds <- DESeqDataSetFromMatrix(
  countData = pseudo,
  colData = meta_pb,
  design = ~ 1)

# Mild filter — HLH genes must survive to be perturbed downstream.
keep <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep, ]
nrow(dds)

dds <- estimateSizeFactors(dds)

vsd <- vst(dds, blind = TRUE)
mat <- assay(vsd)
dim(mat)

# QC the aggregates
# Samples should separate by cell type, not by batch. 
# Batch is confounded with donor in this design, so strong batch separation cannot be distinguished from
# donor variation.
ggsave("21_pca_celltype.png", plotPCA(vsd, intgroup = "celltype"), width = 7, height = 5, dpi = 300)
ggsave("22_pca_batch.png", plotPCA(vsd, intgroup = "batch"), width = 7, height = 5, dpi = 300)

# --- Expressed genes per cell type ---
# Used downstream to restrict the CARNIVAL PPI network to genes actually
# present in the cell type being modelled.
expressed <- lapply(split(seq_len(ncol(mat)), meta_pb$celltype), function(i) {
  rownames(mat)[rowMeans(mat[, i, drop = FALSE]) > 5]
})
sapply(expressed, length)

# Confirm HLH genes are detected in the cell types of interest
hlh <- c(
  "PRF1",   
  "UNC13D",  
  "STX11",  
  "STXBP2", 
  "RAB27A", 
  "LYST",    
  "SH2D1A",  
  "XIAP"
)

hlh[!hlh %in% rownames(mat)]

round(rowMeans(mat[intersect(hlh, rownames(mat)),
                   meta_pb$celltype %in% c("CD8 T effector", "NK")]), 2)

# ============================================================
# 11. TF activity (decoupleR + CollecTRI)
# ============================================================

# Get collecTRI regulon data
collectri <- get_collectri(organism = "human", split_complexes = FALSE)

# Cell-type contrast
# No condition variable exists in this dataset (all healthy donors, single
# timepoint), so TF activity is scored on each cell type's deviation from the
# mean across all cell types. This gives a signed contrast: what is up or down
# in NK relative to the PBMC compartment as a whole.
#
# Working from pseudobulk rather than per-cell also avoids the dropout noise
# that makes single-cell regulon scoring unreliable — most regulon targets are
# zero in any given cell, so per-cell scores are dominated by detection rate.
ct_means <- sapply(split(seq_len(ncol(mat)), meta_pb$celltype), function(i) {
  rowMeans(mat[, i, drop = FALSE])
})

# Center each gene across cell types.
ct_contrast <- ct_means - rowMeans(ct_means)

# TF activity calculation
# For each cell type and TF, ULM fits a linear regression where the cell type's
# centered expression values (response) are modelled against that TF's regulon,
# each target carrying a mode-of-regulation weight (+1 activation, -1 repression).
# The fitted t-value becomes the TF activity score.
#
# A high positive score means the TF's activating targets are elevated and its
# repressed targets reduced in that cell type relative to the PBMC average,
# implying the TF is active there.
#
# ULM is univariate: each TF is tested independently, which is fast and robust,
# though it does not account for TFs sharing target genes.
tf_acts <- run_ulm(
  mat = ct_contrast, # centered pseudobulk matrix (genes x cell types)
  net = collectri, # prior TF-target network (regulons)
  .source = "source", # column naming the TF
  .target = "target", # column naming each target gene
  .mor = "mor", # mode of regulation / weight
  minsize = 5) # keep TFs with at least 5 measured targets

# Compute BH-adjusted p-values on the ULM results
tf_acts <- tf_acts %>%
  dplyr::filter(statistic == 'ulm') %>% # keep ULM rows only
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) # FDR across all TF-cell tests

# TF x celltype matrix
tf_mat <- tf_acts %>%
  pivot_wider(id_cols = source, names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()

# Select the most differentially active TFs across cell types.
# Why SD? For each TF, mean = activity in each cell type. The SD of means
# across cell types measures how much a TF's activity varies between cell types.
# High SD = active in some types but not others.
# Low SD = roughly constant everywhere. 
# Ranking by SD shows the TFs that best separate cell types.
sd_ranked <- data.frame(
  source = rownames(tf_mat),
  std = apply(tf_mat, 1, sd)) %>%
  arrange(desc(std)) %>%
  mutate(rank = row_number())

top_sd <- ggplot(sd_ranked, aes(rank, std)) +
  geom_line() +
  geom_point(size = 0.6) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = "red") +
  labs(x = "TF rank by between-cell-type SD",
       y = "SD of mean activity across cell types",
       title = "Selection of differentially active TFs by standard deviation (SD)") +
  theme_minimal()

ggsave("23_tf_sd_elbow.png", top_sd, width = 7, height = 5, dpi = 300)

# SD of mean activity plateaus around 50, thus select the top 50 TFs for downstream analysis.

# Heatmap of differentially active TFs
top_tfs <- sd_ranked %>%
  slice_head(n = 50) %>%
  pull(source)

top_acts_mat <- t(tf_mat[top_tfs, ])

colors <- rev(brewer.pal(11, "RdBu"))
colors.use <- colorRampPalette(colors)(100)

lim <- quantile(abs(top_acts_mat), 0.95)
my_breaks <- c(seq(-lim, 0, length.out = 51),
               seq(lim / 50, lim, length.out = 50))
pheatmap(
  mat = top_acts_mat,
  color = colors.use,
  border_color = "white",
  breaks = my_breaks,
  cellwidth = 15,
  cellheight = 15,
  treeheight_row = 20,
  treeheight_col = 20,
  main = "Differential transcription factor activity across cell types",
  angle_col = 90,
  fontsize = 10,
  fontsize_row = 11,
  fontsize_col = 8,
  legend = TRUE,
  filename = "tf_heatmap_pseudo.pdf",
  width = 16, height = 5)

# Check against known lineage TFs 
# If these do not land where expected, the contrast construction is wrong
lineage_tfs <- c("TBX21", "EOMES", "RUNX3", "STAT4", # cytotoxic (NK, CD8 eff)
                 "PAX5", "POU2AF1", "EBF1", # B
                 "SPI1", "CEBPB", "CEBPA", # myeloid
                 "TCF7", "LEF1", "FOXP1", # naive T
                 "IRF7", "IRF8") # pDC

round(tf_mat[intersect(lineage_tfs, rownames(tf_mat)), ], 2)

# Top TFs per cell type
tf_acts %>%
  group_by(condition) %>%
  slice_max(abs(score), n = 10) %>%
  arrange(condition, desc(abs(score))) %>%
  print(n = Inf)

# How are primary HLH-associated genes expressed across cell types?
# Canonical fHL-associated genes and other primary HLH-associated genes
hlh <- c(
  "PRF1",   
  "UNC13D",  
  "STX11",  
  "STXBP2", 
  "RAB27A", 
  "LYST",    
  "SH2D1A",  
  "XIAP"
)

hlh <- factor(hlh, levels = hlh)

DefaultAssay(pbmc) <- "RNA"

dot_hlh <- DotPlot(pbmc, features = hlh) +
  RotatedAxis() +
  labs(title = "Expression of primary HLH-associated genes across cell types",
       x = "HLH-associated gene", y = "Cell type") +
  theme(plot.title  = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 9)) +
  geom_vline(xintercept = c(4.5, 6.5), linetype = "dashed", colour = "grey70")

ggsave("hlh_expression__psuedo_dotplot.png", dot_hlh, width = 8, height = 5)

# Same at pseudobulk level, for consistency with the TF scores.
pheatmap(t(ct_means[intersect(hlh, rownames(ct_means)), ]),
         scale = "column",
         cellwidth = 20,
         cellheight = 15,
         color = colors.use,
         border_color = "white",
         angle_col = 45,
         fontsize_col = 9,
         main = "HLH gene expression (pseudobulk, scaled per gene)",
         filename = "hlh_expression_pseudo.png",
         width = 8, height = 5)

# Plot activity of top 5 TFs by cell type
# First, select significant, most-characteristic TFs per cell type ----
top5_tfs <- tf_acts %>%
  dplyr::filter(p_adj < 0.05) %>%
  group_by(condition) %>%
  slice_max(abs(score), n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  pull(source) %>%
  unique()

top5_tfs <- intersect(top5_tfs, rownames(tf_mat))

# Fixed lineage ordering for rows for better readability
ct_order <- c("CD4 T naive", "CD4 T memory", "CD8 T naive", "CD8 T memory",
              "CD8 T effector", "NK", "B naive", "B memory",
              "pDC", "cDC2", "CD14 Mono", "CD14 Mono IFN", "CD16 Mono")
ct_order <- intersect(ct_order, colnames(tf_mat)) # guard against name mismatch

top5_mat <- t(tf_mat[top5_tfs, ct_order, drop = FALSE])

# Colour scale with less saturation 
limit <- quantile(abs(top5_mat), 0.99)
breaks <- c(seq(-lim, 0, length.out = 51),
               seq(lim / 50, lim, length.out = 50))

# Column annotation by functional module
module_map <- c(
  RFX5 = "MHC-II", RFXAP = "MHC-II", RFXANK = "MHC-II", CIITA = "MHC-II", RFX1 = "MHC-II",
  SPI1 = "Myeloid", CEBPA = "Myeloid", CEBPG = "Myeloid", CEBPB = "Myeloid",
  JUN = "Myeloid", JUND = "Myeloid", NFE2L2 = "Myeloid", SP1 = "Myeloid", PPARD = "Myeloid",
  EOMES = "Cytotoxic", STAT4 = "Cytotoxic", ZGLP1 = "Cytotoxic",
  TBX21 = "Cytotoxic", RUNX3 = "Cytotoxic",
  EBF1 = "B lineage", PAX5 = "B lineage", POU2AF1 = "B lineage",
  TRERF1 = "T lineage", RORC = "T lineage", NFKB2 = "T lineage",
  FOXP1 = "T lineage", SATB1 = "T lineage", IKZF1 = "T lineage",
  TCF7 = "T lineage", LEF1 = "T lineage", ETS1 = "T lineage",
  STAT1 = "Interferon", IRF1 = "Interferon", IRF5 = "Interferon",
  IRF7 = "Interferon", IRF8 = "Interferon", RELA = "NF-kB"
)

tf_annot <- data.frame(
  Module = unname(module_map[colnames(top5_mat)]),
  row.names = colnames(top5_mat)
)

module_cols <- list(Module = c(
  "MHC-II" = "#4C72B0",
  "Myeloid" = "#DD8452",
  "Cytotoxic" = "#C44E52",
  "B lineage" = "#55A868",
  "T lineage" = "#8172B3",
  "Interferon" = "#937860",
  "NF-kB" = "#8C8C8C"
))

pheatmap(
  mat = top5_mat,
  color = colors.use,
  breaks = my_breaks,
  border_color = "white",
  cluster_rows = FALSE,         
  cluster_cols = TRUE,          
  annotation_col = tf_annot,
  annotation_colors = module_cols,
  cellwidth = 20,
  cellheight = 20,
  treeheight_col = 20,
  main = "Cell-type-specific transcription factor activity (FDR < 0.05)",
  angle_col = 90,
  fontsize_col = 9,
  fontsize_row   = 10,
  filename = "24_tf_top5_heatmap.png",
  width = 9, height = 5
)

# ============================================================
# 12. Prepare CARNIVAL inputs
# ============================================================

Idents(pbmc) <- "celltype"

# Import prior knowledge network from OmniPath
ppi <- omnipath_interactions() # 85,217 interactions

# Build a signed, directed prior knowledge network (PKN) for CARNIVAL from OmniPath PPIs

# Signalling layer includes directed and unambiguously signed interactions.
# Curation_effort >= 2 was chosen.
# >=3 loses key interactions involved in HLH gene regulation.
# >=1 gives ~70k edges which is too big for solving a network.
sig <- ppi %>%
  dplyr::filter(consensus_direction == 1, # keep interactions with an agreed direction (source -> target)
                consensus_stimulation + consensus_inhibition == 1, # keep only signed edges
                curation_effort >= 2) %>% # keep only edges supported by >= 2 curation sources
  dplyr::mutate(interaction = ifelse(consensus_stimulation == 1, 1, -1)) %>% # encode sign as CARNIVAL expects: +1 activation, -1 inhibition
  dplyr::select(source = source_genesymbol, # regulator gene (edge start)
                target = target_genesymbol, # target gene (edge end)
                interaction) %>% # the signed edge weight
  dplyr::distinct() # remove duplicate edges

# nrow(sig) returns 14,947 interactions survived

# Can HLH genes be used as CARNIVAL inputs?
# Which are in the signalling layer at all?
hlh[hlh %in% c(sig$source, sig$target)]  # STX11, STXBP2, SH2D1A, XIAP

# Which are absent entirely?
setdiff(hlh, c(sig$source, sig$target)) # PRF1, UNC13D, RAB27A, LYST

# Which have outgoing edges? A source node needs these for CARNIVAL to propagate a perturbation
sig %>% 
  dplyr::filter(source %in% hlh) %>% 
  dplyr::count(source) 

# STXBP2 has 1 outgoing edges and XIAP has 5

# Which are sinks? Sinks have regulators (incoming edges) but regulate nothing (no outgoing edges)
intersect(hlh, setdiff(sig$target, sig$source)) # STX11 and SH2D1A

# Of the two with outgoing edges, do they reach any of the measured TFs?
# A source that reaches no target contributes nothing
g_sig <- igraph::graph_from_data_frame(sig %>% 
                                         dplyr::select(source, target),
                                       directed = TRUE)

