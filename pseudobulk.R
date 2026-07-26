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
    meta.data = as.data.frame(colData(kotliarov))) # pulls the per-cell metadata
  
  # Inspect metadata columns to find donor and batch columns needed for pseudobulking later
  colnames(pbmc@meta.data)
  
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
  #   nCount_RNA   = how many transcripts (UMIs) total in that cell
  #   percent.mt   = the mitochondrial fraction computed above
  # Interpretation: very low nFeature = empty or failed droplet; very high
  # nCount = possible two cells in one droplet (doublet); high percent.mt = dying.
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
                 subset = nFeature_RNA > 200 & # drop empty or failed droplets
                  nFeature_RNA < 4000 & # drop likely doublets
                  percent.mt < 10) # drop dying cells
  
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
  # For each cell, find its k nearest neighbors in PC space, then weight the edges beteen two cells by how many
  # neighbors they share. 
  # dims = 1:x should match the number of PCs chosen from ElbowPlot.
  pbmc <- FindNeighbors(pbmc, dims = 1:15)
  res_seq <- seq(0.3, 1.0, by = 0.1)
  
  # Leiden clustering finds groups of cells that are more densely
  # connected to each other than to the rest of the graph, and unlike Louvain it
  # guarantees every cluster is internally connected.
  # n.iter = 10 repeats the refinement to reach a stable partition.
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
  
  # Compare the two algorithms at each resolution.
  # ARI (adjusted Rand index) measures how much two clusterings agree.
  # 1 = identical grouping of cells.
  # 0 = no more agreement than random. 
  # High ARI means the cluster structure is a property of the data, not of the algorithm chosen.
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
  
  # extract the log-normalised counts for singleR annotation
  expression <- GetAssayData(pbmc, layer = "data")
  
  # Monaco: a bulk RNA-seq reference of sorted human immune populations.
  reference <- celldex::MonacoImmuneData()
  
  # How SingleR works, in two steps:
  #  1. For each cell, rank its genes by expression and compute Spearman
  #     correlation against every reference cell type. Ranks are used so the result
  #     does not depend on the absolute scale of either dataset.
  #  2. Take the top-scoring candidate labels and re-score using only genes that
  #     distinguish those specific candidates, repeating until one label wins.
  #     This fine-tuning is what separates similar subtypes.
  #
  # label.main = broad categories (T cells, monocytes)
  pred.main <- SingleR(
    test = expression, 
    ref = reference, 
    labels = reference$label.main)
  
  # label.fine = detailed subtypes (naive CD8 T, effector memory CD8 T)
  pred.fine <- SingleR(
    test = expression, 
    ref = reference, 
    labels = reference$label.fine)
  
  # Cross-tabulate to see how the fine labels nest inside the main ones.
  table(pred.main$labels, pred.fine$labels)
  
  pbmc$main <- pred.main$labels
  pbmc$fine <- pred.fine$labels
  
  p_main <- DimPlot(pbmc, group.by = "main", label = TRUE, repel = TRUE) +
    NoLegend() + ggtitle("label.main")
  
  p_fine <- DimPlot(pbmc, group.by = "fine", label = TRUE, repel = TRUE, label.size = 2.5) +
    NoLegend() + ggtitle("label.fine")
  
  ggsave("08_singler_main_vs_fine.png", p_main + p_fine, width = 14, height = 6, dpi = 300)
  
  # pruned.labels sets low-confidence calls to NA instead of forcing a label.
  # Preferred over $labels because a forced guess is worse than an admitted gap.
  pbmc$singler.main <- pred.main$pruned.labels
  pbmc$singler.fine <- pred.fine$pruned.labels
  
  # Score heatmap: each cell's correlation against every reference label.
  # One bright row per cell = confident call. Several bright rows = the reference
  # cannot tell those subtypes apart at this level of detail.
  png("09_singler_scoreheatmap.png", width = 10, height = 8, units = "in", res = 300)
  plotScoreHeatmap(pred.main)
  dev.off()
  
  # Delta = the assigned label's score minus the median score across all other
  # labels. It measures how much better the winner was than the field.
  # Large delta = confident. Small delta = the call was marginal and gets pruned.
  # If an entire label's delta distribution sits low, that cell type is probably
  # not actually present in the sample.
  ggsave("10_singler_delta.png", plotDeltaDistribution(pred.main, ncol = 4),
         width = 12, height = 8, dpi = 300)
  summary(is.na(pred.main$pruned.labels))
  
  # ============================================================
  # 7. Cross-check labels against clusters
  # ============================================================
  
  # SingleR labels each cell independently; clustering groups cells by similarity.
  # If the two agree, each cluster should be dominated by one label.
  tab <- table(Assigned = pred.main$pruned.labels, Cluster = pbmc$leiden_res.0.7)
  
  # margin = 2 converts to proportions within each column, i.e. what fraction of
  # each cluster carries each label. Without this, big clusters dominate the colour
  # scale regardless of how pure they are.
  png("11_singler_cluster_crosstab.png", width = 8, height = 6, units = "in", res = 300)
  pheatmap(prop.table(tab, margin = 2))
  dev.off()
  
  # ============================================================
  # 8. Protein-based annotation refinement
  # ============================================================
  
  set.seed(42)
  
  # Monaco's plain "T cells" label is not a real population. It is the set of cells
  # whose subtype could not be resolved from RNA. Surface protein resolves them,
  # because the markers immunologists use to define T cell subsets are proteins.
  Idents(pbmc) <- "leiden_res.0.7"
  DefaultAssay(pbmc) <- "ADT"
  
  # Lineage panel: the minimum set needed to assign a cell to a major lineage.
  # Note CD3 is decisive for T cells. A CD3-negative cell is not a T cell no matter
  # how much CD4 it has, because monocytes also carry CD4 protein.
  # CD14 and CD16 together split monocytes into classical, intermediate and
  # non-classical.
  adt_lineage <- c("CD3-PROT", "CD4-PROT", "CD8-PROT", "CD19-PROT",
                   "CD56-PROT", "CD14-PROT", "CD16-PROT", "HLA-DR-PROT")
  
  # DotPlot: dot colour = average expression in that cluster, dot size = fraction
  # of cells in the cluster expressing it. Both matter, because a high average
  # driven by a few cells means something different from uniform expression.
  ggsave("12_adt_lineage_dotplot.png",
         DotPlot(pbmc, features = adt_lineage) + RotatedAxis(),
         width = 8, height = 7, dpi = 300)
  
  # Stacked violins show the full distribution, which reveals bimodality (a
  # genuinely positive subset within a cluster) that a dot plot average hides.
  # pt.size = 0 hides individual points so the shapes stay readable.
  ggsave("13_adt_lineage_violin.png",
         VlnPlot(pbmc, adt_lineage, stack = TRUE, flip = TRUE, pt.size = 0) + NoLegend(),
         width = 9, height = 8, dpi = 300)
  
  
  # Subset markers: 
  # memory/naive state (CD45RA, CD45RO, CD62L, CD197)
  # regulatory (CD25, CD127), 
  # cytotoxic (CD57),
  # MAIT (CD161), 
  # gamma-delta (TCRgd), 
  # exhaustion (CD279), 
  # dendritic subsets
  # (CD11c, CD123, CD1c, CD303), 
  #B cell state (CD20, CD27, IgD, IgM).
  adt_subset <- c("CD45RA-PROT", "CD45RO-PROT", "CD62L-PROT", "CD197-PROT",
                  "CD127-PROT", "CD25-PROT", "CD57-PROT", "CD161-PROT",
                  "TCRgd-PROT", "CD279-PROT", "CD183-PROT",
                  "CD11c-PROT", "CD123-PROT", "CD1c-PROT", "CD303-PROT",
                  "CD20-PROT", "CD27-PROT", "IgD-PROT", "IgM-PROT", "CD34-PROT")
  
  ggsave("14_adt_subset_dotplot.png",
         DotPlot(pbmc, features = adt_subset) + RotatedAxis(),
         width = 12, height = 7, dpi = 300)
  
  # Isotype controls: antibodies that should bind nothing. They measure the
  # background level of the assay, so any real signal must sit clearly above them.
  # Signal at isotype level is a true negative, not a weak positive.
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
  
  # Doublet and quality checks.
  # A cluster showing markers of two mutually exclusive lineages (CD3 and CD14
  # together, or both CD4-high and CD8-high) is probably two cells captured in one
  # droplet. Elevated counts and genes relative to neighbouring clusters supports
  # that, since two cells contribute roughly twice the RNA.
  ggsave("18_counts_features.png",
         VlnPlot(pbmc, c("nCount_RNA", "nFeature_RNA"), pt.size = 0, ncol = 1) + NoLegend(),
         width = 10, height = 7, dpi = 300)
  
  ggsave("19_percent_mt.png",
         VlnPlot(pbmc, "percent.mt", pt.size = 0) + NoLegend(),
         width = 8, height = 5, dpi = 300)
  
  # Unbiased markers: for each cluster, test every gene for higher expression in
  # that cluster than in all other cells combined (one-vs-rest, Wilcoxon by default).
  all_markers <- FindAllMarkers(
    pbmc,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    max.cells.per.ident = 500,
    random.seed = 42)
  
  write.csv(all_markers, "all_markers_leiden_res07.csv", row.names = FALSE)
  
  # Top markers per cluster: significant, expressed in most of the cluster
  # (pct.1 > 0.5, so it describes the cluster as a whole rather than a fragment),
  # ranked by fold change.
  top_markers <- all_markers %>%
    filter(p_val_adj < 0.05, pct.1 > 0.5) %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 10) %>%
    ungroup()
  
  print(top_markers, n = Inf)
  
  # Heatmap of those markers. downsample = 100 takes 100 cells per cluster, so the
  # plot stays legible and small clusters are not visually swamped by large ones.
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
  
  # Map cluster IDs to labels. unname() strips the cluster-number names that
  # indexing carries over, then cell barcodes are attached instead.
  lab <- unname(final_label[as.character(pbmc$leiden_res.0.7)])
  names(lab) <- colnames(pbmc)
  pbmc$celltype <- lab
  
  # Fix the display order: related populations next to each other, artefacts last. 
  # Without an explicit factor order, plots and tables will sort alphabetically.
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
  
  # Remove non-cells and failed populations. Platelets are real but anucleate with
  # almost no transcriptome, so they are not informative here.
  # invert = TRUE keeps everything except the listed types.
  pbmc.clean <- subset(pbmc, subset = celltype %in%
                         c("Low quality", "Doublets", "Platelets"), invert = TRUE)
  
  saveRDS(pbmc.clean, "pbmc_annotated.rds")
  
  # ============================================================
  # 10. Pseudobulk
  #
  # Individual cells are not independent replicates. Cells from one donor share
  # that donor's genetics and handling. Treating them as independent inflates
  # significance. Summing counts per cell type per donor gives one value per
  # biological replicate, which is what a bulk method like DESeq2 expects.
  # ============================================================
  
  # Donor 209 has cells in both batches.
  # Kept as two separate samples.
  cells_donor <- table(pbmc$sampleid, pbmc$batch)
  names(which(rowSums(cells_donor > 0) > 1))  # donor 209
  cells_donor["209", ]                        # 1188 cells batch 1, 1977 batch 2
  
  # The pseudobulk grouping key: donor and batch together.
  pbmc.clean$pb_sample <- paste(pbmc.clean$sampleid, pbmc.clean$batch, sep = "_")
  
  # Merge the activated state back into classical monocytes.
  # Cluster 13 is an IFN-stimulated state, not a cell type: 531 of its 625
  # cells come from donor 209 alone, and only 3 donors contribute >=10 cells.
  # Merge the activated monocyte state back into classical monocytes.
  # Cluster 13 is an interferon-stimulated state, not a cell type: 531 of its 625
  # cells come from donor 209 alone, and only 3 donors contribute >=10 cells. Any
  # contrast on it would describe one individual rather than a population.
  pbmc.clean$celltype[pbmc.clean$celltype == "activated classical monocytes"] <-
    "classical monocytes"
  pbmc.clean$celltype <- droplevels(factor(pbmc.clean$celltype))
  
  # AggregateExpression sums raw counts across all cells sharing a cell type and
  # sample. Summing (not averaging) is required because DESeq2 models integer
  # counts and needs the total depth to estimate variance correctly.
  pseudo <- AggregateExpression(
    pbmc.clean,
    assays = "RNA",
    slot = "counts", # raw counts, not normalised values
    group.by = c("celltype", "pb_sample"), # one output column per combination
    return.seurat = FALSE)$RNA  # return a plain matrix
  
  # Rebuilding the metadata for those columns.
  # AggregateExpression names columns by replacing any "_" inside each group value
  # with "-", then joining the two values with "_". Parsing that back out is
  # error-prone, so the key is reconstructed from the metadata and matched instead.
  n_cells <- as.data.frame(
    table(celltype  = droplevels(factor(pbmc.clean$celltype)),
          pb_sample = factor(pbmc.clean$pb_sample)),
    stringsAsFactors = FALSE)
  
  n_cells$column <- paste(gsub("_", "-", n_cells$celltype),
                          gsub("_", "-", n_cells$pb_sample), sep = "_")
  
  # Fails loudly if any pseudobulk column could not have been produced by a
  # celltype-sample pair, which would mean the naming assumption above is wrong.
  stopifnot(all(colnames(pseudo) %in% n_cells$column))
  
  # Align metadata rows to matrix columns.
  meta_pb <- n_cells[match(colnames(pseudo), n_cells$column), ]
  rownames(meta_pb) <- meta_pb$column
  names(meta_pb)[names(meta_pb) == "Freq"] <- "n_cells"
  
  # Drop groups built from too few cells: a sum over 3 cells is not a stable
  # expression profile and adds noise without adding information.
  meta_pb <- meta_pb[meta_pb$n_cells >= 10, ]
  pseudo <- pseudo[, rownames(meta_pb)]
  
  meta_pb$celltype <- factor(meta_pb$celltype)
  meta_pb$pb_sample <- factor(meta_pb$pb_sample)
  meta_pb$sampleid <- sub("_.*$", "", meta_pb$pb_sample)
  meta_pb$batch <- factor(sub("^.*_", "", meta_pb$pb_sample))
  
  table(meta_pb$celltype)
  
  # ============================================================
  # 11. Differential expression
  #
  # Two contrasts, each comparing a cytotoxic effector to a naive baseline:
  #   NK vs pooled naive (CD4 + CD8 naive T)
  #   CD8 TEMRA vs CD8 naive T
  #
  # The naive populations are the resting counterparts, so the contrast isolates
  # the effector programme rather than lineage differences.
  # ============================================================
  
  DefaultAssay(pbmc.clean) <- "RNA"
  
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

dot_hlh <- DotPlot(pbmc.clean, features = levels(hlh)) +
  RotatedAxis() +
  labs(title = "Expression of primary HLH-associated genes across cell types",
       x = "HLH-associated gene",
       y = "Cell type") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10,), 
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 9)) + 
  geom_vline(xintercept = c(4.5, 6.5), linetype = "dashed", colour = "grey70")

ggsave("22_hlh_dotplot.png", dot_hlh, width = 8, height = 5)


# The four cell types used for the contrasts
effectors <- c("NK", "CD8 TEMRA")
naive <- c("CD8 naive T", "CD4 naive T")
used <- c(effectors, naive)

# Subset the pseudobulk matrix and metadata to those cell types, keeping the two
# objects aligned by using the metadata rownames to index the matrix columns.
keep_cols <- meta_pb$celltype %in% used
meta_sub <- droplevels(meta_pb[keep_cols, , drop = FALSE])
pseudo_sub <- pseudo[, rownames(meta_sub), drop = FALSE]

# Set CD8 naive T as the reference level, so model coefficients are expressed as
# "cell type X versus CD8 naive T" by default.
meta_sub$celltype <- relevel(factor(meta_sub$celltype), ref = "CD8 naive T")

# Keep only donors contributing at least two cell types.
# The design includes pb_sample as a term, which means each donor's baseline is
# estimated from that donor's own samples. A donor with only one column gives no
# within-donor comparison and makes the model unidentifiable.
cols_per_sample <- table(meta_sub$pb_sample)
usable_samples <- names(cols_per_sample)[cols_per_sample >= 2]
dropped <- setdiff(levels(factor(meta_sub$pb_sample)), usable_samples)

meta_sub <- droplevels(meta_sub[meta_sub$pb_sample %in% usable_samples, , drop = FALSE])
pseudo_sub <- pseudo_sub[, rownames(meta_sub), drop = FALSE]

# design = ~ pb_sample + celltype
#   pb_sample: a paired/blocking term. Each donor gets their own intercept, so
#   donor-to-donor differences are absorbed and do not inflate the
#   variance used to test cell type.
#   celltype: the term of interest, tested after that adjustment.
dds_ct <- DESeqDataSetFromMatrix(
  countData = pseudo_sub,
  colData = meta_sub,
  design = ~ pb_sample + celltype
)

# Gene filter: keep genes with >=5 counts in at least as many samples as the
# smallest cell type group. This ensures a gene can  be detected
# throughout one whole group, so a real group-specific gene is not removed,
# while genes seen in only a couple of samples are.
min_grp <- min(table(meta_sub$celltype))
keep <- rowSums(counts(dds_ct) >= 5) >= min_grp
dds_ct <- dds_ct[keep, ]

# Size factors correct for differences in total sequencing depth between samples.
# type = "poscounts" computes them using only non-zero counts per gene, which is
# needed here because pseudobulk matrices still contain many zeros and the
# default method requires a gene to be non-zero in every sample.
dds_ct <- estimateSizeFactors(dds_ct, type = "poscounts")

# Variance stabilising transformation.
# Raw counts have variance that grows with the mean, so high-expressed genes
# dominate distance calculations. VST puts everything on a roughly log scale with
# constant variance, which is what PCA and clustering need.
# blind = TRUE ignores the design, so the transform cannot be accused of having
# been shaped by the comparison it is used to inspect.
vsd_ct <- vst(dds_ct, blind = TRUE)
mat_ct <- assay(vsd_ct)

# PCA of the pseudobulk samples: checks that cell types separate and that batch
# is not the dominant source of variation.
pca_df <- plotPCA(vsd_ct, intgroup = c("celltype", "batch"), returnData = TRUE)
pv <- round(100 * attr(pca_df, "percentVar"))
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = celltype, shape = batch)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(x = sprintf("PC1: %d%%", pv[1]), y = sprintf("PC2: %d%%", pv[2]),
       title = "Lymphoid pseudobulk (VST)") +
  theme_bw()
ggsave("23_pca_lymphoid.png", p_pca, width = 8, height = 5, dpi = 300)

# Fit the model. DESeq() runs three steps: size factor estimation, dispersion
# estimation (how variable each gene is beyond what counting noise explains),
# and the negative binomial fit with Wald tests.
dds_ct <- DESeq(dds_ct, quiet = TRUE)

# Dispersion plot: Black points are per-gene estimates, the red
# line is the fitted trend, blue points are the final shrunken values. Genes
# should sit near the trend and shrink toward it; a cloud with no trend means
# the model has not fitted well.
png("24_dispersions_lymphoid.png", width = 7, height = 6, units = "in", res = 300)
plotDispEsts(dds_ct, main = "Dispersion plot: lymphoid fit")
dev.off()

# Sample-to-sample distances on the VST values. Samples of the same cell type
# should cluster together. Useful for spotting a mislabelled or outlier sample.
sampleDists <- dist(t(mat_ct))
ann <- as.data.frame(colData(dds_ct)[, c("celltype", "batch")])
png("25_sample_distance_lymphoid.png", width = 9, height = 8, units = "in", res = 300)
pheatmap(as.matrix(sampleDists),
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         annotation_col = ann, show_rownames = FALSE, show_colnames = FALSE,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample distances (VST)")
dev.off()

# Contrast 1: CD8 TEMRA vs CD8 naive T 
# A comparison of two levels of the same factor. Written as
# c(factor, numerator, denominator), so a positive fold change means higher in
# CD8 TEMRA.
# alpha = 0.05 sets the FDR level used for independent filtering and the summary.
res_temra <- results(dds_ct,
                     contrast = c("celltype", "CD8 TEMRA", "CD8 naive T"),
                     alpha = 0.05)

# Contrast 2: NK vs the average of two naive types
# With CD8 naive T as the reference, the group means are:
#   CD8 naive T = b0
#   CD4 naive T = b0 + b_CD4
#   NK = b0 + b_NK
# The pooled naive mean is the average of the two naive groups:
#   pooled = b0 + b_CD4 / 2
# So the difference of interest is:
#   NK - pooled = b_NK - 0.5 * b_CD4
# which is the weight vector: +1 on b_NK, -0.5 on b_CD4, 0 elsewhere.
rn <- resultsNames(dds_ct)

# DESeq2 mangles level names into coefficient names (spaces become dots), so the
# names are looked up rather than typed by hand. stopifnot fails loudly if the
# pattern matches zero or more than one coefficient.
coef_of <- function(lvl) {
  hit <- grep(paste0("^celltype_", make.names(lvl), "_vs_"), rn, value = TRUE)
  stopifnot(length(hit) == 1L)
  hit
}

con <- setNames(numeric(length(rn)), rn)
con[coef_of("NK")] <-  1
con[coef_of("CD4 naive T")] <- -0.5
# CD8 naive T is the reference, its coefficient is 0 by construction

# unname() because results() expects a plain numeric vector in coefficient order.
res_nk <- results(dds_ct, contrast = unname(con), alpha = 0.05)

# Volcano plots for the two contrasts
hlh_chr <- as.character(hlh)

# Well-known genes included purely for orientation, so a reader can confirm the
# plot is the right way round: cytotoxic genes on the effector side, naive
# markers (CCR7, SELL, TCF7) on the other.
anchors <- c("GNLY", "NKG7", "GZMB", "CCL5", "CD3D", "CCR7", "SELL", "TCF7")

volcano <- function(res, title, file, cap = 150) {
  
  df <- as.data.frame(res) |>
    rownames_to_column("gene") |>
    filter(!is.na(padj)) |> # DESeq2 sets padj to NA for filtered genes
    mutate(
      y = pmin(-log10(padj), cap), # Caps the axis to prevent small p-values dominating
      sig = case_when(padj < 0.05 & log2FoldChange > 1 ~ "Up in effector",
                      padj < 0.05 & log2FoldChange < -1 ~ "Up in naive",
                      TRUE ~ "n.s."))
  
  ggplot(df, aes(log2FoldChange, y)) +
    geom_point(aes(colour = sig), size = 0.7, alpha = 0.4) +
    scale_colour_manual(values = c("Up in effector" = "#C0504D",
                                   "Up in naive" = "#4F81BD",
                                   "n.s." = "grey80"), name = NULL) +
    geom_vline(xintercept = c(-1, 1), lty = 2, colour = "grey60") +
    geom_hline(yintercept = -log10(0.05), lty = 2, colour = "grey60") +
    
    # anchor markers, for orientation
    geom_text_repel(data = filter(df, gene %in% anchors),
                    aes(label = gene), size = 3, colour = "grey25",
                    fontface = "italic", max.overlaps = Inf, seed = 42,
                    ylim = c(NA, cap * 0.9), nudge_y = -8,
                    box.padding = 0.5,
                    segment.colour = "grey60", segment.size = 0.3) +
    
    # HLH genes, the subject
    geom_point(data = filter(df, gene %in% hlh_chr),
               colour = "#1B7837", size = 2.5, shape = 18) +
    geom_text_repel(data = filter(df, gene %in% hlh_chr),
                    aes(label = gene), size = 3.5, colour = "#1B7837",
                    fontface = "bold", max.overlaps = Inf, seed = 42,
                    ylim = c(NA, cap * 0.9), nudge_y = -8,
                    box.padding = 0.5,
                    segment.colour = "#1B7837", segment.size = 0.3) +
    
    labs(x = "log2 fold change", y = "-log10 (padj)", title = title,
         caption = sprintf("y-axis capped at %g", cap)) +
    theme_bw()
  
  ggsave(file, width = 9, height = 7, dpi = 300)
}

volcano(res_temra, "CD8 TEMRA vs CD8 naive T", "26_volcano_temra.png")
volcano(res_nk, "NK vs pooled naive (CD4 + CD8 naive T)", "27_volcano_nk.png")

# Extract the Wald statistics
# stat = log2FoldChange / standard error. It combines effect size and precision
# in one signed number: a large fold change measured noisily gets a small
# statistic, whereas the fold change alone would not distinguish the two.
# This is the right input for TF activity inference, which needs a ranking that
# already accounts for uncertainty.
stat_temra <- setNames(res_temra$stat, rownames(res_temra))
stat_nk <- setNames(res_nk$stat, rownames(res_nk))
stat_temra <- stat_temra[!is.na(stat_temra)]
stat_nk <- stat_nk[!is.na(stat_nk)]

saveRDS(stat_temra, "stat_temra.rds")
saveRDS(stat_nk, "stat_nk.rds")
write.csv(as.data.frame(res_temra), "de_temra.csv")
write.csv(as.data.frame(res_nk), "de_nk.csv")

# ============================================================
# 12. Transcription factor activity (decoupleR + CollecTRI)
#
# A TF's own mRNA level is a poor measure of its activity, because activity is
# controlled by phosphorylation, localisation and cofactors. Instead, activity is
# inferred from the behaviour of the genes it regulates.
# ============================================================

# CollecTRI is a curated set of TF-target relationships. Each edge carries a "mor"
# (mode of regulation): +1 if the TF activates that target, -1 if it represses.
# split_complexes = FALSE keeps complexes such as NF-kB as one entry rather than
# splitting them into subunits.
collectri <- get_collectri(organism = "human", split_complexes = FALSE)

# Activity per cell type
# Collapse the pseudobulk samples to one mean profile per cell type.
# tapply() groups the values of each gene by cell type and averages within group.
ct <- droplevels(factor(meta_sub[colnames(mat_ct), "celltype"]))
mat_ct_mean <- t(apply(mat_ct, 1, function(x) tapply(x, ct, mean)))

# Centre and scale each gene across the four cell types (z-score per row).
#
# Why this is necessary: the linear model below regresses expression on regulon
# membership. On raw VST values the fit is driven by absolute expression level,
# which is shared by all four cell types, so every TF would score high. After
# centring, each value says "how far above or below this gene's own average",
# which is the difference the analysis is about.
#
# What this means for interpretation: scores are relative to the mean of THESE
# FOUR cell types, not to a PBMC-wide average. A TF uniformly high in all
# lymphocytes but absent in monocytes now looks flat, not high.
mat_z <- t(scale(t(mat_ct_mean)))
mat_z <- mat_z[stats::complete.cases(mat_z), ]

# Univariate linear model (ULM).
# For each TF and each cell type, fit: expression ~ regulon membership, where
# membership is 0 for genes the TF does not regulate, +1 for activated targets
# and -1 for repressed ones. The t-value of the fitted slope is the activity score.
#
# A high positive score means activated targets are up and repressed targets are
# down, which is what an active TF produces.
#
# "Univariate" means each TF is tested on its own. That is fast and robust, but
# it does not account for TFs that share target genes, so overlapping regulons
# can give correlated scores.
#
# minsize = 5 excludes TFs with fewer than 5 measured targets, where the slope
# would be estimated from too little data to trust.
tf_acts <- run_ulm(
  mat = mat_z, 
  network = collectri,
  .source = "source", 
  .target = "target", 
  .mor = "mor",
  minsize = 5) %>%
  dplyr::filter(statistic == "ulm") %>%
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH"))

# Reshape to a TF x cell type matrix for plotting.
tf_mat <- tf_acts %>%
  pivot_wider(id_cols = source, names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()

# Ranked by largest absolute score, not by significance and not by variance.
#
# Not by FDR: with only four conditions, the z-scoring above compresses the
# t-values, so few tests clear BH correction. Those that do are the TFs with the
# largest regulons (MYC, SP1, JUN), because a bigger regulon gives a more precise
# slope, not a more meaningful one.
#
# Not by standard deviation: with n = 4, and two of those (CD4 and CD8 naive)
# being near-identical, an SD over four numbers is unstable and cannot tell
# "high in one cell type" apart from "splits 2 versus 2".
tf_ranked <- tf_acts %>%
  group_by(source) %>%
  summarise(max_abs = max(abs(score)), # strongest activity in any cell type
            best_padj = min(p_adj), .groups = "drop") %>%
  arrange(desc(max_abs))

top_tfs <- head(sig$source, 40)

# Transpose so cell types are rows and TFs are columns.
top_acts_mat <- t(tf_mat[top_tfs, ])

# Fix the row order naive -> effector, so the gradient reads down the plot
# instead of being rearranged by clustering.
row_order <- c("CD4 naive T", "CD8 naive T", "CD8 TEMRA", "NK")
top_acts_mat <- top_acts_mat[intersect(row_order, rownames(top_acts_mat)), ]

# module annotation
# Assigned from known biology
modules <- list(
  "Effector / IFN" = c("TBX21", "STAT1", "IRF1", "IRF3", "IRF5",
                           "NFKB", "RUNX1", "RUNX2", "PML"),
  "MHC class II" = c("CIITA", "RFX5", "RFXANK", "RFXAP"),
  "Growth / metabolic" = c("MYC", "MZF1", "SP1", "NFYA", "NFYB",
                           "SREBF1", "SREBF2", "NR1H3", "ATF6"),
  "Stress / AP-1" = c("JUN", "AP1", "ATF2", "DDIT3", "NFE2L2",
                           "TP53", "BACH1", "LITAF"))

# Everything starts as "Other" and is reassigned if it appears in a module.
module <- setNames(rep("Other", ncol(top_acts_mat)), colnames(top_acts_mat))
for (m in names(modules))
  module[names(module) %in% modules[[m]]] <- m

# pheatmap needs a data frame whose rownames match the column names of the matrix.
ann_col <- data.frame(Module = factor(
  module, levels = c(names(modules), "Other")),
  row.names = names(module))

ann_colors <- list(Module = c(
  "Effector / IFN" = "#C0504D",
  "MHC class II" = "#9970AB",
  "Growth / metabolic" = "#4F81BD",
  "Stress / AP-1" = "#E8A33D",
  "Other" = "grey85"))

table(ann_col$Module)   # how many fell through to "Other"

# significance stars 
# Recovers the FDR information dropped by ranking on effect size:
# the reader can see both magnitude (colour) and significance (star) together.
star_mat <- tf_acts %>%
  dplyr::filter(source %in% colnames(top_acts_mat)) %>%
  dplyr::mutate(lab = dplyr::case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01  ~ "**",
    p_adj < 0.05  ~ "*",
    TRUE ~ "")) %>%
  pivot_wider(id_cols = condition, names_from = source, values_from = lab) %>%
  column_to_rownames("condition") %>%
  as.matrix()

# Reorder to exactly match the heatmap, or the stars would land on wrong cells.
star_mat <- star_mat[rownames(top_acts_mat), colnames(top_acts_mat)]

# heatmap 
colors.use <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
lim <- quantile(abs(top_acts_mat), 0.95)
my_breaks <- c(seq(-lim, 0, length.out = 51),
               seq(lim / 50, lim, length.out = 50))

pheatmap(
  mat = top_acts_mat,
  color = colors.use,
  breaks = my_breaks,
  border_color = "white",
  cellwidth = 15, 
  cellheight = 15,
  cluster_rows = FALSE,         
  gaps_row = 2,                  
  cutree_cols = 4,               
  treeheight_col = 25,
  angle_col = 90,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  display_numbers = star_mat,
  number_color = "black",
  fontsize_number = 7,
  fontsize_row = 11,
  fontsize_col = 8,
  main = "TF activity: cytotoxic effectors vs naive T",
  filename = "28_tf_acts_lymphoid.pdf",
  width = ncol(top_acts_mat) * 15/72 + 5,
  height = nrow(top_acts_mat) * 15/72 + 3.5)

# TF activity per contrast
# The heatmap above describes cell types. This describes the two comparisons,
# which is what the network modelling needs: each column is already a difference
# from naive, so TFs shared by NK and TEMRA show up in both rather than
# cancelling against each other.

# The two contrasts were filtered independently, so their gene sets differ.
# Intersecting first prevents misalignment when they are bound into one matrix.
g <- intersect(names(stat_temra), names(stat_nk))
stat_mat <- cbind("CD8 TEMRA" = stat_temra[g], "NK" = stat_nk[g])

# No centring here. Unlike the expression matrix, these values are already
# differences relative to naive, so the baseline is built into the contrast.
# Centring across two columns would remove exactly the signal being measured.
tf_contrast <- run_ulm(
  mat = stat_mat, 
  network = collectri,
  .source = "source", 
  .target = "target", 
  .mor = "mor",
  minsize = 5) %>%
  dplyr::filter(statistic == "ulm") %>%
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH"))

# ============================================================
# 12. Prior knowledge network and CARNIVAL inputs
#
#     Signalling layer -> couple TF->HLH edges -> prune per cell type
#     The CollecTRI edges must be present before pruning so the
#     expression filter applies to them too, and baseline and
#     counterfactual runs must share identical topology so that a node
#     flip reflects the constraint, not a change in the network.
# 
# ============================================================

# TF activities as CARNIVAL measurements

# CARNIVAL is an ILP; runtime scales with the number of
# measurements, so take the most confident TFs rather than everything
# that clears FDR.
make_measobj <- function(contrast, n_top = 50) {
  tf_contrast %>%
    dplyr::filter(condition == contrast, p_adj < 0.05) %>%
    dplyr::slice_max(abs(score), n = n_top) %>%
    dplyr::select(source, score) %>%
    tidyr::pivot_wider(names_from = source, values_from = score) %>%
    as.data.frame()
}

meas_temra <- make_measobj("CD8 TEMRA")
meas_nk <- make_measobj("NK")

tf_all <- union(colnames(meas_temra), colnames(meas_nk))


# Signalling layer from OmniPath
ppi <- omnipath_interactions()      

sig_all <- ppi %>%
  dplyr::filter(consensus_direction == 1, # agreed direction source -> target
                consensus_stimulation + consensus_inhibition == 1) %>% # unambiguously signed
  dplyr::mutate(interaction = ifelse(consensus_stimulation == 1, 1L, -1L)) %>%
  dplyr::select(source = source_genesymbol, interaction,
                target = target_genesymbol, curation_effort) %>%
  dplyr::filter(source != "", # unmapped symbols would collapse into one phantom hub
                target != "",
                source != target) %>% # self-loops: one state per node, so vacuous or contradictory
  dplyr::group_by(source, interaction, target) %>%
  dplyr::summarise(curation_effort = max(curation_effort), .groups = "drop")
# nrow(sig_all) = 70,565 unique signed directed edges


# Choosing the curation and expression thresholds
# Two candidate gene universes:
#   strict = the DESeq2 filter (>=5 counts in >= min_grp samples). Designed
#   for stable dispersion estimates, so it is stricter than a
#   presence filter needs to be and drops low-abundance
#   signalling proteins.
#   lax = presence filter (>=1 count in >=5 lymphoid pseudobulk samples).
expressed_ct <- rownames(mat_ct)
expressed_lax <- rownames(pseudo_sub)[rowSums(pseudo_sub >= 1) >= 5]
c(strict = length(expressed_ct), lax = length(expressed_lax))

check <- function(ce, universe, label) {
  s <- sig_all %>%
    dplyr::filter(curation_effort >= ce,
                  source %in% universe, target %in% universe)
  nodes  <- unique(c(s$source, s$target))
  in_deg <- table(s$target)
  reach  <- intersect(tf_all, nodes)
  data.frame(universe = label, curation = ce,
             edges = nrow(s), nodes = length(nodes),
             tf_present = length(reach),
             # a TF with no incoming edge can never be explained by any solution
             tf_with_in = sum(reach %in% names(in_deg)))
}

grid <- do.call(rbind, c(
  lapply(1:4, check, universe = expressed_ct,  label = "strict"),
  lapply(1:4, check, universe = expressed_lax, label = "lax")))
grid

# Decision: lax universe, curation_effort >= 3.
#   The expression filter dominates: strict excludes 9-10 reachable TFs at
#   every curation level (32 vs 42 at curation effort >=1).
#   Within the lax universe, curation effort >=3 is the elbow: relaxing to >=2 buys one
#   extra TF for 2,331 extra edges; tightening to >=4 saves 1,283 edges but
#   costs two TFs.
sig <- sig_all %>%
  dplyr::filter(curation_effort >= 3,
                source %in% expressed_lax, target %in% expressed_lax) %>%
  dplyr::select(source, interaction, target)

nrow(sig)


# 12.4 Couple the transcriptional layer
# The HLH genes are terminal effectors. They have regulators but no downstream
# signalling. Adding TF -> HLH edges from CollecTRI lets them enter the
# model as measurable endpoints rather than as perturbation nodes.
# Require the TF to already sit in the signalling layer, otherwise it
# enters as a free upstream parameter with nothing constraining it.
trn <- collectri %>%
  dplyr::filter(target %in% hlh_chr,
                source %in% c(sig$source, sig$target)) %>%
  dplyr::mutate(interaction = as.integer(mor)) %>%
  dplyr::select(source, interaction, target) %>%
  dplyr::distinct()

pkn <- dplyr::bind_rows(sig, trn) %>% dplyr::distinct()

c(signalling = nrow(sig), transcriptional = nrow(trn), combined = nrow(pkn))

# restrict measurements to nodes the PKN contains
pkn_nodes <- unique(c(pkn$source, pkn$target))
meas_temra <- meas_temra[, colnames(meas_temra) %in% pkn_nodes, drop = FALSE]
meas_nk <- meas_nk[, colnames(meas_nk) %in% pkn_nodes, drop = FALSE]
dim(meas_temra); dim(meas_nk)

saveRDS(pkn, "pkn_full.rds")


# 12.5 Audit: how can each HLH gene enter the model?
# An input node needs outgoing edges that reach the
# measured TFs, so the perturbation can propagate.
# A measurement node needs only incoming edges, so a sink can be used.
g_pkn <- igraph::graph_from_data_frame(
  pkn %>% dplyr::select(source, target), directed = TRUE)

reach_to_tfs <- function(node) {
  if (!node %in% igraph::V(g_pkn)$name) return(0L)
  d <- igraph::distances(g_pkn, v = node,
                         to = intersect(tf_all, igraph::V(g_pkn)$name),
                         mode = "out")
  as.integer(sum(is.finite(d)))
}

hlh_audit <- tibble::tibble(gene = hlh_chr) %>%
  dplyr::mutate(
    in_signalling = gene %in% c(sig$source, sig$target),
    n_regulators = vapply(gene, function(g) sum(pkn$target == g), integer(1)),
    n_outgoing = vapply(gene, function(g) sum(pkn$source == g), integer(1)),
    tfs_reached = vapply(gene, reach_to_tfs, integer(1)),
    role = dplyr::case_when(
      n_outgoing > 0 & n_regulators > 0 ~ "both",
      n_outgoing > 0 ~ "source only",
      n_regulators > 0 ~ "sink (measurable only)",
      TRUE ~ "absent"),
    usable_as_input = ifelse(n_outgoing > 0 & tfs_reached > 0, "yes", "no")) %>%
  dplyr::arrange(dplyr::desc(n_outgoing), dplyr::desc(n_regulators))

hlh_audit


# Evidence that the exclusions are structural, not artefacts
# PRF1 has no usable outgoing edges anywhere in OmniPath, at any curation
# level: two of its four are unsigned, two have zero references.
all_int <- OmnipathR::import_all_interactions()

all_int %>%
  dplyr::filter(source_genesymbol == "PRF1") %>%
  dplyr::select(target_genesymbol, is_directed, is_stimulation,
                is_inhibition, curation_effort, n_references)

# usable outgoing edges per HLH gene across ALL OmniPath layers
all_int %>%
  dplyr::filter(source_genesymbol %in% hlh_chr,
                is_directed == 1,
                is_stimulation + is_inhibition == 1) %>%
  dplyr::count(source_genesymbol, name = "usable_outgoing") %>%
  dplyr::right_join(tibble::tibble(source_genesymbol = hlh_chr),
                    by = "source_genesymbol") %>%
  dplyr::mutate(usable_outgoing = tidyr::replace_na(usable_outgoing, 0L)) %>%
  dplyr::arrange(dplyr::desc(usable_outgoing))

# Genes with no regulators in CollecTRI either cannot enter at all
collectri %>% 
  dplyr::filter(target %in% hlh_chr) %>% 
  dplyr::count(target)

setdiff(hlh_chr, collectri$target)

# The two layers join: PRF1's regulators are embedded in the signalling layer
pkn %>%
  dplyr::filter(target == "PRF1") %>%
  dplyr::mutate(regulator_in_signalling = source %in% c(sig$source, sig$target))

# End-to-end reachability: cytokine receptor -> signalling -> TF -> PRF1
receptors <- intersect(c("IL12RB1","IL12RB2","IFNGR1","IFNGR2","IL2RB","IL18R1"),
                       igraph::V(g_pkn)$name)
sapply(receptors, function(r) {
  d <- igraph::distances(g_pkn, v = r, to = "PRF1", mode = "out")
  ifelse(is.finite(d), as.numeric(d), NA_real_)
})

# the canonical JAK-STAT route, recovered
igraph::V(g_pkn)$name[
  igraph::shortest_paths(g_pkn, from = "IL12RB1", to = "PRF1",
                         mode = "out")$vpath[[1]]]


# cell-type-specific regulation of PRF1
# Used ONLY to ask which edges differ between cell types, not to build
# the networks that get solved. min_pct = 0.1 is a marker-gene threshold
# (Seurat's FindMarkers default); with scRNA-seq dropout it is far too
# strict to use as a presence filter, and pruning the solved PKN with it
# removed two thirds of the measured TFs. The expression filter for the
# PKN was already applied at the pseudobulk level (>=1 count in >=5
# lymphoid samples).
expressed_in <- function(ct, min_pct = 0.05) {
  cells <- colnames(pbmc.clean)[pbmc.clean$celltype == ct]
  cnt <- GetAssayData(pbmc.clean, assay = "RNA", layer = "counts")[, cells]
  rownames(cnt)[Matrix::rowMeans(cnt > 0) >= min_pct]
}

genes_nk    <- expressed_in("NK")
genes_temra <- expressed_in("CD8 TEMRA")

prf1_regs <- pkn %>% dplyr::filter(target == "PRF1") %>%
  dplyr::mutate(in_nk    = source %in% genes_nk,
                in_temra = source %in% genes_temra)
prf1_regs

# regulators detected in NK but not CD8 TEMRA
prf1_regs %>% dplyr::filter(in_nk, !in_temra) %>% dplyr::pull(source)


# CARNIVAL inputs: one shared PKN, per-contrast measurements
# Both cell types are solved over the identical coupled PKN, so any
# difference between the solved networks reflects the measurements, not
# the topology. The same holds for baseline vs counterfactual later.
saveRDS(pkn, "pkn_carnival.rds")

stat_list <- list(nk = stat_nk, temra = stat_temra)
tf_list   <- list(nk = meas_nk, temra = meas_temra)
tags      <- c("nk", "temra")

for (tag in tags) {
  
  tf_v <- setNames(as.numeric(tf_list[[tag]][1, ]), colnames(tf_list[[tag]]))
  tf_v <- tf_v[names(tf_v) %in% pkn_nodes]
  
  h <- stat_list[[tag]][intersect(hlh_chr, pkn_nodes)]
  h <- h[!is.na(h)]
  
  # measObj magnitudes enter the ILP objective. If the HLH Wald statistics
  # span a wider range than the TF activity scores, a few genes dominate
  # the fit, so rescale them onto the TF range.
  cat(sprintf("%-6s TF [%.1f, %.1f] n=%d | HLH raw [%.1f, %.1f] n=%d\n",
              tag, min(tf_v), max(tf_v), length(tf_v),
              min(h), max(h), length(h)))
  
  if (max(abs(h)) > max(abs(tf_v)))
    h <- h * max(abs(tf_v)) / max(abs(h))
  
  meas <- as.data.frame(as.list(c(tf_v, h)))
  saveRDS(meas, sprintf("meas_%s_baseline.rds", tag))
  
  cat(sprintf("%-6s %d TFs + %d HLH genes = %d measurements\n",
              tag, length(tf_v), length(h), ncol(meas)))
}

# ---- validation ----
for (tag in tags) {
  p <- readRDS("pkn_carnival.rds")
  m <- readRDS(sprintf("meas_%s_baseline.rds", tag))
  stopifnot(identical(colnames(p), c("source", "interaction", "target")),
            all(p$interaction %in% c(-1L, 1L)),
            nrow(m) == 1L,
            all(colnames(m) %in% c(p$source, p$target)))
  cat(tag, "OK:", nrow(p), "edges,", ncol(m), "measurements\n")
}