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
  cells_donor <- table(pbmc$sampleid, pbmc$batch)
  names(which(rowSums(cells_donor > 0) > 1))  # donor 209
  cells_donor["209", ]                        # 1188 cells batch 1, 1977 batch 2
  
  pbmc.clean$pb_sample <- paste(pbmc.clean$sampleid, pbmc.clean$batch, sep = "_")
  
  # Merge the activated state back into classical monocytes.
  # Cluster 13 is an IFN-stimulated state, not a cell type: 531 of its 625
  # cells come from donor 209 alone, and only 3 donors contribute >=10 cells.
  # A one-vs-rest contrast on it would describe one individual.
  pbmc.clean$celltype[pbmc.clean$celltype == "activated classical monocytes"] <-
    "classical monocytes"
  pbmc.clean$celltype <- droplevels(factor(pbmc.clean$celltype))
  
  # AggregateExpression() sums all counts across all cells sharing celltype
  # and sample. Returns features x groups. Summing keeps raw counts for DESeq2.
  pseudo <- AggregateExpression(
    pbmc.clean,
    assays = "RNA",
    slot = "counts",
    group.by = c("celltype", "pb_sample"),
    return.seurat = FALSE)$RNA
  
  # Column names are built by rewriting any "_" inside each group value to "-",
  # then joining the two values with "_". Rebuild that key from the metadata
  # rather than parsing it back out of the column names.
  n_cells <- as.data.frame(
    table(celltype  = droplevels(factor(pbmc.clean$celltype)),
          pb_sample = factor(pbmc.clean$pb_sample)),
    stringsAsFactors = FALSE)
  
  n_cells$column <- paste(gsub("_", "-", n_cells$celltype),
                          gsub("_", "-", n_cells$pb_sample), sep = "_")
  
  # Fires only if a pseudobulk column exists that no celltype-donor pair
  # could have produced, i.e. the mangling assumption above is wrong.
  stopifnot(all(colnames(pseudo) %in% n_cells$column))
  
  meta_pb <- n_cells[match(colnames(pseudo), n_cells$column), ]
  rownames(meta_pb) <- meta_pb$column
  names(meta_pb)[names(meta_pb) == "Freq"] <- "n_cells"
  
  # Drop groups aggregated from too few cells to be a stable profile
  meta_pb <- meta_pb[meta_pb$n_cells >= 10, ]
  pseudo <- pseudo[, rownames(meta_pb)]
  
  meta_pb$celltype <- factor(meta_pb$celltype)
  meta_pb$pb_sample <- factor(meta_pb$pb_sample)
  meta_pb$sampleid <- sub("_.*$", "", meta_pb$pb_sample)
  meta_pb$batch <- factor(sub("^.*_", "", meta_pb$pb_sample))
  
  table(meta_pb$celltype)
  
  # ============================================================
  # 11. Differential expression for regulatory-network inference
  #
  #     Two healthy cytotoxic networks will be built downstream:
  #       - NK : contrast NK vs pooled naive (CD4 + CD8 naive T)
  #       - CD8 TEMRA : contrast  CD8 TEMRA vs CD8 naive T
  #
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

effectors <- c("NK", "CD8 TEMRA")
naive <- c("CD8 naive T", "CD4 naive T")
used <- c(effectors, naive)

# subset the pseudobulked counts + metadata to the cell types selected
keep_cols <- meta_pb$celltype %in% used
meta_sub <- droplevels(meta_pb[keep_cols, , drop = FALSE])
pseudo_sub <- pseudo[, rownames(meta_sub), drop = FALSE]

meta_sub$celltype <- relevel(factor(meta_sub$celltype), ref = "CD8 naive T")

cols_per_sample <- table(meta_sub$pb_sample)
usable_samples <- names(cols_per_sample)[cols_per_sample >= 2]
dropped <- setdiff(levels(factor(meta_sub$pb_sample)), usable_samples)

meta_sub <- droplevels(meta_sub[meta_sub$pb_sample %in% usable_samples, , drop = FALSE])
pseudo_sub <- pseudo_sub[, rownames(meta_sub), drop = FALSE]

# build the DESeq2 object
dds_ct <- DESeqDataSetFromMatrix(
  countData = pseudo_sub,
  colData = meta_sub,
  design = ~ pb_sample + celltype
)

min_grp <- min(table(meta_sub$celltype))
keep <- rowSums(counts(dds_ct) >= 5) >= min_grp
dds_ct <- dds_ct[keep, ]

dds_ct <- estimateSizeFactors(dds_ct, type = "poscounts")

# VST for visualisation
vsd_ct <- vst(dds_ct, blind = TRUE)
mat_ct <- assay(vsd_ct)

# PCA: cell types cleanly separate
pca_df <- plotPCA(vsd_ct, intgroup = c("celltype", "batch"), returnData = TRUE)
pv <- round(100 * attr(pca_df, "percentVar"))
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = celltype, shape = batch)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(x = sprintf("PC1: %d%%", pv[1]), y = sprintf("PC2: %d%%", pv[2]),
       title = "Lymphoid pseudobulk (VST)") +
  theme_bw()
ggsave("23_pca_lymphoid.png", p_pca, width = 8, height = 5, dpi = 300)

# run DESeq2
dds_ct <- DESeq(dds_ct, quiet = TRUE)

png("24_dispersions_lymphoid.png", width = 7, height = 6, units = "in", res = 300)
plotDispEsts(dds_ct, main = "Dispersion plot: lymphoid fit")
dev.off()

# Sample-distance heatmap (shared, VST-based)
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

# CD8 TEMRA vs CD8 naive T — direct level-vs-level, ref is already CD8 naive T
res_temra <- results(dds_ct,
                     contrast = c("celltype", "CD8 TEMRA", "CD8 naive T"),
                     alpha = 0.05)

# NK vs pooled naive — an average of two levels, so it needs a numeric contrast.
# Group means: CD8naive = b0, CD4naive = b0 + b_CD4, NK = b0 + b_NK
# Pooled naive  = b0 + b_CD4/2   ->   NK - pooled = b_NK - 0.5 * b_CD4
rn <- resultsNames(dds_ct)

coef_of <- function(lvl) {
  hit <- grep(paste0("^celltype_", make.names(lvl), "_vs_"), rn, value = TRUE)
  stopifnot(length(hit) == 1L)
  hit
}

con <- setNames(numeric(length(rn)), rn)
con[coef_of("NK")] <-  1
con[coef_of("CD4 naive T")] <- -0.5
# CD8 naive T is the reference, its coefficient is 0 by construction

res_nk <- results(dds_ct, contrast = unname(con), alpha = 0.05)

hlh_chr <- as.character(hlh)
anchors <- c("GNLY", "NKG7", "GZMB", "CCL5", "CD3D", "CCR7", "SELL", "TCF7")

volcano <- function(res, title, file, cap = 150) {
  
  df <- as.data.frame(res) |>
    rownames_to_column("gene") |>
    filter(!is.na(padj)) |>
    mutate(
      y = pmin(-log10(padj), cap),
      sig = case_when(padj < 0.05 & log2FoldChange >  1 ~ "Up in effector",
                      padj < 0.05 & log2FoldChange < -1 ~ "Up in naive",
                      TRUE ~ "n.s."))
  
  ggplot(df, aes(log2FoldChange, y)) +
    geom_point(aes(colour = sig), size = 0.7, alpha = 0.4) +
    scale_colour_manual(values = c("Up in effector" = "#C0504D",
                                   "Up in naive"    = "#4F81BD",
                                   "n.s."           = "grey80"), name = NULL) +
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

stat_temra <- setNames(res_temra$stat, rownames(res_temra))
stat_nk <- setNames(res_nk$stat, rownames(res_nk))
stat_temra <- stat_temra[!is.na(stat_temra)]
stat_nk <- stat_nk[!is.na(stat_nk)]

saveRDS(stat_temra, "stat_temra.rds")
saveRDS(stat_nk, "stat_nk.rds")
write.csv(as.data.frame(res_temra), "de_temra.csv")
write.csv(as.data.frame(res_nk), "de_nk.csv")

# ============================================================
# 11. TF activity (decoupleR + CollecTRI)
# ============================================================
collectri <- get_collectri(organism = "human", split_complexes = FALSE)

# collapse pseudobulk samples to cell-type means
ct <- droplevels(factor(meta_sub[colnames(mat_ct), "celltype"]))
mat_ct_mean <- t(apply(mat_ct, 1, function(x) tapply(x, ct, mean)))

# gene-centering and scaling
# Scores are relative to the mean of these four cell types, not to a
# PBMC-wide average. Without centering, ULM tracks absolute expression
# level (shared across all four) rather than difference between them.
mat_z <- t(scale(t(mat_ct_mean)))
mat_z <- mat_z[stats::complete.cases(mat_z), ]

tf_acts <- run_ulm(
  mat = mat_z, network = collectri,
  .source = "source", .target = "target", .mor = "mor",
  minsize = 5) %>%
  dplyr::filter(statistic == "ulm") %>%
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH"))

tf_mat <- tf_acts %>%
  pivot_wider(id_cols = source, names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()

# select TFs by effect size 
# Not by FDR: with four conditions, the z-scoring compresses t-values,
# so few tests clear BH and those that do are the largest regulons
# (MYC, SP1, JUN) rather than the most informative TFs. 
# Not by SD: n = 4, two of which (CD4/CD8 naive) are near-duplicates, so
# SD cannot distinguish "high in one type" from "splits 2 vs 2".
sig <- tf_acts %>%
  group_by(source) %>%
  summarise(max_abs   = max(abs(score)),
            best_padj = min(p_adj), .groups = "drop") %>%
  arrange(desc(max_abs))

top_tfs <- head(sig$source, 40)

top_acts_mat <- t(tf_mat[top_tfs, ])
row_order <- c("CD4 naive T", "CD8 naive T", "CD8 TEMRA", "NK")
top_acts_mat <- top_acts_mat[intersect(row_order, rownames(top_acts_mat)), ]

# module annotation
# Assigned from prior biology
modules <- list(
  "Effector / IFN" = c("TBX21", "STAT1", "IRF1", "IRF3", "IRF5",
                           "NFKB", "RUNX1", "RUNX2", "PML"),
  "MHC class II" = c("CIITA", "RFX5", "RFXANK", "RFXAP"),
  "Growth / metabolic" = c("MYC", "MZF1", "SP1", "NFYA", "NFYB",
                           "SREBF1", "SREBF2", "NR1H3", "ATF6"),
  "Stress / AP-1" = c("JUN", "AP1", "ATF2", "DDIT3", "NFE2L2",
                           "TP53", "BACH1", "LITAF"))

module <- setNames(rep("Other", ncol(top_acts_mat)), colnames(top_acts_mat))
for (m in names(modules))
  module[names(module) %in% modules[[m]]] <- m

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

# combine the Wald Statistics for both contrasts into one single matrix
g <- intersect(names(stat_temra), names(stat_nk))
stat_mat <- cbind("CD8 TEMRA" = stat_temra[g], "NK" = stat_nk[g])

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
# 12. Prepare CARNIVAL inputs
# ============================================================

# Extract TF activities
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

# export as rds for carnival
saveRDS(meas_temra, "meas_temra.rds")
saveRDS(meas_nk, "meas_nk.rds")

# Import prior knowledge network from OmniPath
ppi <- omnipath_interactions() # 85,217 interactions

# Build a signed, directed prior knowledge network (PKN) for CARNIVAL from OmniPath PPIs

# Signalling layer includes directed and unambiguously signed interactions.
# Curation_effort >= 2 was chosen.
# >=3 loses key interactions involved in HLH gene regulation.
# >=1 gives ~70k edges which is too big for solving a network.
sig_all <- ppi %>%
  dplyr::filter(consensus_direction == 1,
                consensus_stimulation + consensus_inhibition == 1) %>%
  dplyr::mutate(interaction = ifelse(consensus_stimulation == 1, 1L, -1L)) %>%
  dplyr::select(source = source_genesymbol, interaction,
                target = target_genesymbol, curation_effort) %>%
  dplyr::filter(source != "", target != "", source != target) %>%
  dplyr::group_by(source, interaction, target) %>%         
  dplyr::summarise(curation_effort = max(curation_effort), .groups = "drop")

# nrow(sig_all) returns 70,565 interactions

expressed_ct <- rownames(mat_ct)                                        
expressed_lax <- rownames(pseudo_sub)[rowSums(pseudo_sub >= 1) >= 5]     
c(strict = length(expressed_ct), lax = length(expressed_lax))

# ---- grid: curation threshold x universe ----
tf_all <- union(colnames(meas_temra), colnames(meas_nk))

check <- function(ce, universe, label) {
  s <- sig_all %>%
    dplyr::filter(curation_effort >= ce,
                  source %in% universe, target %in% universe)
  nodes <- unique(c(s$source, s$target))
  in_deg <- table(s$target)
  reach <- intersect(tf_all, nodes)
  data.frame(
    universe = label,
    curation = ce,
    edges = nrow(s),
    nodes = length(nodes),
    tf_present = length(reach),
    tf_with_in = sum(reach %in% names(in_deg)))   # TFs that can actually be explained
}

grid <- do.call(rbind, c(
  lapply(1:4, check, universe = expressed_ct,  label = "strict"),
  lapply(1:4, check, universe = expressed_lax, label = "lax")))
grid



