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

# ============================================================
# 1. Load data
# ============================================================

kotliarov <- KotliarovPBMCData(
  mode = c("rna", "adt"),
  ensembl = FALSE,
  location = TRUE,
  legacy = FALSE)

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

# subset pbmcs
pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 &
                nFeature_RNA < 4000 &
                percent.mt < 10)

# ============================================================
# 3. Attach ADT
# ============================================================

# Doing this ensures the protein assay is subset to the same cells
# Inspect the feature names
adt_counts <- counts(altExp(kotliarov, "ADT"))
rownames(adt_counts)

cells <- colnames(pbmc)
adt_counts <- adt_counts[, cells]
pbmc[["ADT"]] <- CreateAssayObject(counts = adt_counts)

# CLR across cells. Protein counts are compositional with high background,
# so RNA-style log-normalisation is not appropriate.
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
pbmc <- RunPCA(pbmc, npcs = 50)
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

# Cluster cells by community detection such as louvain on that SNN graph.
# Resolution controls granularity: higher = more, granular clusters. 
pbmc <- FindClusters(pbmc, resolution = seq(0.3, 1.0, by = 0.1))

# Compute a 2D UMAP embedding from the same 15 PCs for visualization.
# Nonlinear projection that places similar cells near each other in a 2D space, capturing local structure.
# Use the same dims as FindNeighbors so the graph and the plot are consistent.
pbmc <- RunUMAP(pbmc, dims = 1:15)

res_cols <- grep("^RNA_snn_res", colnames(pbmc@meta.data), value = TRUE)
ggsave("05_umap_resolutions.png",
       DimPlot(pbmc, group.by = res_cols, label = TRUE, ncol = 4) & NoLegend(),
       width = 16, height = 8, dpi = 300)

# Find how cell clusters separate with each resolution
ggsave("06_clustree.png", clustree(pbmc), width = 10, height = 10, dpi = 300)

Idents(pbmc) <- "RNA_snn_res.0.7"
ggsave("07_umap_res07.png", DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend(),
       width = 7, height = 6, dpi = 300)

# ============================================================
# 6. Reference-based annotation
# ============================================================

expression <- GetAssayData(pbmc, layer = "data")

reference <- celldex::MonacoImmuneData()

prediction <- SingleR(
  test = expression, 
  ref = reference, 
  labels = reference$label.main)

pbmc$singler <- prediction$pruned.labels

png("08_singler_scoreheatmap.png", width = 8, height = 8, units = "in", res = 300)
plotScoreHeatmap(prediction)
dev.off()

ggsave("09_singler_delta.png", plotDeltaDistribution(prediction, ncol = 3), width = 9, height = 8, dpi = 300)
summary(is.na(prediction$pruned.labels))

ggsave("10_umap_singler.png",
       DimPlot(pbmc, group.by = "singler", label = TRUE, repel = TRUE) + NoLegend(),
       width = 7, height = 6, dpi = 300)

# ============================================================
# 7. Cross-check labels against clusters
# ============================================================

tab <- table(Assigned = prediction$pruned.labels, Cluster = Idents(pbmc))

# Raw counts show which clusters are small enough that proportions mislead.
tab

# Column proportions: what fraction of each cluster carries each label.
# Preferable to log10(tab + 10), which compresses a 50-cell and a 5000-cell
# cluster into similar colours.
round(prop.table(tab, margin = 2), 2)

png("11_singler_cluster_crosstab.png", width = 8, height = 6, units = "in", res = 300)
pheatmap(prop.table(tab, margin = 2))
dev.off()

# ============================================================
# 8. Protein-based annotation refinement
# ============================================================

# Monaco's "T cells" label is not a distinct population, but a compartment of
# cells that could not be resolved by RNA. Resolve with protein instead.

DefaultAssay(pbmc) <- "ADT"

# --- Lineage gate ---
# CD3 first: anything CD3-negative is not a T cell regardless of CD4 level,
# since CD4 protein is genuinely expressed on monocytes.
adt_lineage <- c("CD3-PROT", "CD4-PROT", "CD8-PROT", "CD14-PROT",
                 "CD19-PROT", "CD56-PROT", "CD16-PROT", "HLA-DR-PROT")

ggsave("12_adt_lineage_dotplot.png", DotPlot(pbmc, features = adt_lineage) + RotatedAxis(),
       width = 8, height = 6, dpi = 300)

ggsave("13_adt_lineage_violin.png", VlnPlot(pbmc, adt_lineage, stack = TRUE, flip = TRUE) + NoLegend(),
       width = 8, height = 8, dpi = 300)

# --- Subset markers, once lineage is settled ---
adt_subset <- c("CD45RA-PROT", "CD45RO-PROT", "CD62L-PROT", "CD197-PROT",
                "CD127-PROT", "CD25-PROT", "CD57-PROT", "TCRgd-PROT",
                "CD11c-PROT", "CD123-PROT", "CD1c-PROT", "CD303-PROT",
                "CD20-PROT", "IgD-PROT", "CD34-PROT")

ggsave("14_adt_subset_dotplot.png", DotPlot(pbmc, features = adt_subset) + RotatedAxis(),
       width = 10, height = 6, dpi = 300)

# --- Background floor ---
# Signal sitting at isotype level is true negative; above it is real,
# or ambient contamination. Sets the threshold for the calls above.
isotypes <- grep("sotype", rownames(pbmc[["ADT"]]), value = TRUE)

ggsave("15_isotypes.png", VlnPlot(pbmc, isotypes, stack = TRUE, flip = TRUE) + NoLegend(),
       width = 7, height = 5, dpi = 300)

# --- RNA markers as secondary evidence ---
# CD4 mRNA drops out heavily in 10x and is near-zero even in true CD4 T cells —
# never call CD4 identity from RNA.

DefaultAssay(pbmc) <- "RNA"

markers_rna <- c(
  "CD3D", "CD3E",                  
  "CD8A", "CD8B",                  
  "IL7R", "CCR7", "SELL", "TCF7",  
  "FOXP3", "IL2RA",               
  "GZMK", "GZMB", "NKG7", "PRF1",  
  "KLRD1", "NCAM1",               
  "MS4A1", "CD79A",                
  "LYZ", "CD14", "FCGR3A",      
  "FCER1A", "CLEC9A",              
  "CD34", "SPINK2", "PRSS57")      

ggsave("16_rna_markers_dotplot.png", DotPlot(pbmc, features = markers_rna) + RotatedAxis(),
       width = 12, height = 6, dpi = 300)


# --- Doublet check ---
# Clusters showing mutually exclusive lineage markers (CD3+ and CD14+, or
# CD4-high and CD8-high) are doublet candidates. Elevated counts confirm.
ggsave("17_counts_features.png", VlnPlot(pbmc, c("nCount_RNA", "nFeature_RNA"), pt.size = 0) + NoLegend(),
       width = 8, height = 5, dpi = 300)

ggsave("18_percent_mt.png", VlnPlot(pbmc, "percent.mt", pt.size = 0) + NoLegend(),
       width = 6, height = 5, dpi = 300)

# --- Unbiased markers ---
# Downsampled: marker detection does not need every cell, and 55k is slow.
all_markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  max.cells.per.ident = 500)

top_markers <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10) %>%
  ungroup()

print(top_markers, n = Inf)

# ============================================================
# 9. Final labels
# ============================================================

# Labels assigned from three converging lines of evidence: SingleR/Monaco
# (lineage), ADT surface phenotype (CD3/CD4/CD8/CD14/CD19/CD56/CD16), and
# unbiased RNA markers (section 8). Where they disagreed, protein and canonical
# markers took precedence over the reference call.

cluster_labels <- c(
  "0" = "CD4 T naive", # FHIT, CCR7, LEF1, TCF7, MAL
  "1" = "CD4 T memory", # IL7R, LTB, IL32, GPR183, ITGB1
  "2" = "CD14 Mono", # S100A8/9/12, VCAN, LYZ, CD14
  "3" = "CD8 T effector", # GZMH, FGFBP2, KLRG1, NKG7, CCL5
  "4" = "NK", # KLRF1, SPON2, PRF1, GNLY, GZMB, CD56+CD16+
  "5" = "B naive", # TCL1A, FCER2, VPREB3, MS4A1
  "6" = "CD8 T memory", # GZMK, KLRB1, CD8B, DUSP2
  "7" = "CD8 T naive", # CD8B, NELL2, CCR7, LEF1
  "8" = "Low quality", # MALAT1 + MT- genes only; ~10% mito, lowest counts
  "9" = "B memory", # TNFRSF13B, POU2AF1, BANK1, BLK
  "10" = "CD16 Mono", # CDKN1C, MS4A7, C1QA, CSF1R, CD16-high CD14-low
  "11" = "cDC2", # FCER1A, CLEC10A, CD1C
  "12" = "CD14 Mono IFN", # FCGR1A/B, GBP1, WARS, APOBEC3A — activated state
  "13" = "pDC", # LILRA4, CLEC4C, SCT, LAMP5
  "14" = "Proliferating", # MKI67, TOP2A, CDK1, UBE2C — cell cycle, not doublets
  "15" = "Mast/basophil prog", # TPSAB1, CPA3, CD34
  "16" = "Platelet") # PF4, PPBP, GP9, ITGA2B, TUBB1

pbmc <- RenameIdents(pbmc, cluster_labels)
pbmc$celltype <- Idents(pbmc)

stopifnot(!any(is.na(pbmc$celltype)))
table(pbmc$celltype)

ggsave("19_umap_celltype.png",
       DimPlot(pbmc, group.by = "celltype", label = TRUE, repel = TRUE, label.size = 3) + NoLegend(),
       width = 8, height = 7, dpi = 300)


# --- Remove non-cell and failed populations ---
drop <- c("Low quality", "Platelet", "Mast/basophil prog", "Proliferating")

pbmc <- subset(pbmc, subset = celltype %in% drop, invert = TRUE)
pbmc$celltype <- droplevels(pbmc$celltype)

table(pbmc$celltype)


# --- Order levels sensibly for plotting ---
pbmc$celltype <- factor(pbmc$celltype, levels = c(
  "CD4 T naive", "CD4 T memory",
  "CD8 T naive", "CD8 T memory", "CD8 T effector",
  "NK",
  "B naive", "B memory",
  "CD14 Mono", "CD14 Mono IFN", "CD16 Mono",
  "cDC2", "pDC"))

Idents(pbmc) <- "celltype"

DimPlot(pbmc, group.by = "celltype", label = TRUE, repel = TRUE, label.size = 3) +
  NoLegend()

# --- Final sanity check ---
# Confirm labels against protein. Each population should be positive for the
# markers that define it and at isotype level for the others.
DefaultAssay(pbmc) <- "ADT"
ggsave("20_adt_final_dotplot.png",
       DotPlot(pbmc, features = c("CD3-PROT", "CD4-PROT", "CD8-PROT", "CD45RA-PROT",
                                  "CD45RO-PROT", "CD62L-PROT", "CD56-PROT", "CD16-PROT",
                                  "CD19-PROT", "IgD-PROT", "CD14-PROT", "HLA-DR-PROT",
                                  "CD123-PROT", "CD34-PROT")) + RotatedAxis(),
       width = 12, height = 6, dpi = 300)
DefaultAssay(pbmc) <- "RNA"

saveRDS(pbmc, "pbmc_pseudobulk_annotated.rds")

# ============================================================
# 10. Pseudobulk
# ============================================================

# Donor 209 has cells in both batches. 
# Kept as two separate samples.
pbmc$pb_sample <- paste(pbmc$sampleid, pbmc$batch, sep = "_")

# Cells per celltype-sample pair
cells_per_group <- table(pbmc$celltype, pbmc$pb_sample)
rowSums(cells_per_group >= 10)

# Aggregate raw counts 
pseudo <- AggregateExpression(
  pbmc,
  assays = "RNA",
  slot = "counts",    
  group.by = c("celltype", "pb_sample"),
  return.seurat = FALSE)$RNA

dim(pseudo)
head(colnames(pseudo))


# --- Sample sheet ---
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

