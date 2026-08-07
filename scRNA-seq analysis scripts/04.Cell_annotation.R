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
score_heatmap <- as.ggplot(plotScoreHeatmap(pred.main, silent = TRUE)$gtable)


# Delta = the assigned label's score minus the median score across all other
# labels. It measures how much better the winner was than the field.
# Large delta = confident. Small delta = the call was marginal and gets pruned.
# If an entire label's delta distribution sits low, that cell type is probably
# not actually present in the sample.
delta_distribution <- plotDeltaDistribution(pred.main, ncol = 4)

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
check_labels <- as.ggplot(pheatmap(prop.table(tab, margin = 2), silent = TRUE)$gtable)

label_check_supp <- score_heatmap / delta_distribution / check_labels +
  plot_layout(heights = c(1, 1, 0.9)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave("Supp_Figure_02_singler_qc.png", label_check_supp,
       width = 10, height = 16, dpi = 300, bg = "white")
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

# RNA markers as secondary evidence
DefaultAssay(pbmc) <- "RNA"

# Canonical PBMC markers
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

# Top markers per cluster: significant, expressed in most of the cluster
# (pct.1 > 0.5, so it describes the cluster as a whole rather than a fragment),
# ranked by fold change.
top_markers <- all_markers %>%
  filter(p_val_adj < 0.05, pct.1 > 0.5) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  ungroup() %>%
  distinct(gene, .keep_all = TRUE)

print(top_markers, n = Inf)

# Heatmap of those markers. downsample = 100 takes 100 cells per cluster, so the
# plot stays legible and small clusters are not visually swamped by large ones.
top_markers_heatmap <- DoHeatmap(subset(pbmc, downsample = 100),
                                 features = top_markers$gene,
                                 group.by = "seurat_clusters",
                                 size = 2.5, angle = 0, hjust = 0.5) +
  theme(axis.text.y = element_text(size = 7)) +
  NoLegend()

ggsave("supp_fig_04_marker_heatmap.png", top_markers_heatmap,
       width = 6.2, height = 0.12 * nrow(top_markers) + 1)

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
  
  # No lineage-defining protein
  # RNA: only MT-* genes and MALAT1
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

evidence <- tibble::tribble(
  ~cluster, ~protein, ~rna, ~label,
  "1", "CD3+ CD4+ CD62L-hi CD197-hi CD45RA+ CD127+", "CCR7, SELL, TCF7, LDHB, PIK3IP1, NOSIP", "CD4 naive T",
  "2", "CD3+ CD4+ CD45RO-hi CD45RA-lo CD25-hi CD127+", "IL7R, LTB, IL32, CD69, ITGB1; CCR7/SELL absent", "CD4 memory T",
  "3", "CD14-hi CD16-neg HLA-DR+ CD11c+", "CD14, LYZ, S100A8/9/12, VCAN, LGALS2, MS4A6A", "Classical monocytes",
  "4", "CD3+ CD8+ CD57-hi CD62L-neg CD45RO+", "GZMH, GNLY, FGFBP2, NKG7, CST7, GZMA/GZMM", "CD8 TEMRA",
  "5", "CD3+ CD8+ CD161-hi CD45RO-hi CD127+ CD279+", "GZMK, KLRB1, DUSP2, IL7R, CCL5", "MAIT / CD8 EM",
  "6", "CD3-neg CD56-hi CD16-hi", "KLRF1, KLRD1, SPON2, CLIC3, PRF1, GZMB, GNLY", "NK",
  "7", "CD19+ CD20-hi IgD-hi IgM-hi CD27-neg", "TCL1A, MS4A1, CD79A/B, HLA-DQ/DR", "Naive B",
  "8", "CD3+ CD8-hi CD62L+ CD197-hi CD45RA+", "CD8B, CCR7, TCF7, SELL, NOSIP, PIK3IP1, LDHB", "CD8 naive T",
  "9", "No lineage-defining protein", "MT-* genes and MALAT1 only; several with pct.2 > pct.1", "Low quality",
  "10", "CD19+ CD20-hi CD27+ IgM+ IgD-lo CD11c+", "BANK1, MS4A1, CD79A/B, IGJ; TCL1A absent", "Memory B",
  "11", "CD16-hi CD14-lo CD11c-hi", "FCGR3A, MS4A7, CDKN1C, CSF1R, LST1, LILRB2, TCF7L2", "Non-classical monocytes",
  "12", "CD1c-hi CD11c+ HLA-DR-hi CD14-neg", "FCER1A, CLEC10A, CD1C, CPVL, HLA-DPA1/DQA1", "cDC2",
  "13", "CD14-hi CD16-neg CD11c+", "S100A8, FCGR1A, FOLR3, GBP1, WARS, TNFSF10, TYMP", "Activated classical monocytes",
  "14", "CD303-hi CD123-hi CD11c-neg HLA-DR+", "LILRA4, CLEC4C, SCT, SERPINF1, DNASE1L3, LRRC26", "pDC",
  "15", "CD3+ CD4+ CD8+ CD19+ CD16+ (mutually exclusive)", "MKI67, TYMS, RRM2, TK1, PCNA, STMN1", "Doublets",
  "16", "All lineage proteins flat", "SDPR/CAVIN2, PPBP, PF4, HIST1H2AC, TSC22D1", "Platelets")

# add cell counts and SingleR majority call, so the table carries data not just prose
evidence <- evidence %>%
  dplyr::mutate(
    n_cells = as.integer(table(pbmc$leiden_res.0.7)[cluster]),
    singler = colnames(singler_by_cluster)[max.col(singler_by_cluster)][
      match(cluster, rownames(singler_by_cluster))])

evidence %>%
  dplyr::select(Cluster = cluster, `n cells` = n_cells,
                `SingleR (majority)` = singler,
                `Surface protein` = protein, `RNA markers` = rna,
                `Assigned label` = label) %>%
  gt::gt() %>%
  gt::tab_header(title = "Evidence supporting cell type assignment") %>%
  gt::tab_source_note("Clusters from Leiden clustering at resolution 0.7") %>%
  gt::gtsave("supp_table_annotation.png")