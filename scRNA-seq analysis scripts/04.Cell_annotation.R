# ============================================================
# Shared plot themes
# ============================================================
umap_theme <- theme(plot.title = element_text(size = 10, face = "bold"),
                    axis.title = element_text(size = 8),
                    axis.text  = element_text(size = 7))

qc_theme <- theme(plot.title   = element_text(size = 11, face = "bold"),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text    = element_text(size = 8))

no_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# ============================================================
# 6. Reference-based annotation
# ============================================================

DefaultAssay(pbmc) <- "RNA"
Idents(pbmc) <- "leiden_res.0.7"

# extract the log-normalised counts for SingleR. 
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

# pruned.labels sets low-confidence calls to NA instead of forcing a label.
# Preferred over $labels because a forced guess is worse than an admitted gap.
pbmc$singler.main <- pred.main$pruned.labels

# How many cells were assinged a label, how many were not?
summary(is.na(pred.main$pruned.labels))

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

# Panels for the main annotation figure
# A. UMAP of Leiden res_0.7: the substrate everything else is mapped onto
unsup_clus <- DimPlot(pbmc, group.by = "leiden_res.0.7",
                      label = TRUE, repel = TRUE, label.size = 3) +
  NoLegend() +
  umap_theme

# B. Reference-based call
ref_label <- DimPlot(pbmc, group.by = "singler.main",
                     label = TRUE, repel = TRUE, label.size = 2.5) +
  NoLegend() +
  umap_theme

# C. Where the RNA reference was not enough. Two failure modes: pruned to NA 
# (low confidence), or given a generic lineage label with no subtype. 
# These are what the protein panel addresses.
generic <- c("T cells", "CD4+ T cells", "CD8+ T cells")
pbmc$anno_gap <- dplyr::case_when(
  is.na(pbmc$singler.main) ~ "pruned (NA)",
  pbmc$singler.main %in% generic ~ "generic label",
  TRUE ~ "resolved")

gen_plot <- DimPlot(pbmc, group.by = "anno_gap",
                    cols = c("resolved" = "grey88",
                             "generic label" = "#C0504D",
                             "pruned (NA)" = "#4F81BD"),
                    order = c("pruned (NA)", "generic label")) +
  umap_theme +
  theme(legend.position = "bottom", legend.text = element_text(size = 7),
        legend.title = element_blank())

# ============================================================
# 7. Cross-check labels against clusters
# ============================================================

# SingleR labels each cell independently; clustering groups cells by similarity.
# If the two agree, each cluster should be dominated by one label.
# Rows = cluster, columns = assigned label.
singler_by_cluster <- table(Cluster  = pbmc$leiden_res.0.7,
                            Assigned = pbmc$singler.main)

# Transposed for plotting so labels are rows and clusters are columns, then
# margin = 2 converts to proportions within each cluster, i.e. what fraction of
# each cluster carries each label. Without this, big clusters dominate the colour
# scale regardless of how pure they are.
check_labels <- as.ggplot(
  pheatmap(prop.table(t(singler_by_cluster), margin = 2), silent = TRUE)$gtable)

# Majority reference call per cluster, for the evidence table in section 9.
# NA where a cluster has no labelled cells at all, so an all-pruned cluster is
# not handed an arbitrary label by max.col().
majority <- colnames(singler_by_cluster)[max.col(singler_by_cluster)]
majority[rowSums(singler_by_cluster) == 0] <- NA_character_
names(majority) <- rownames(singler_by_cluster)

# ============================================================
# 8. Protein-based annotation refinement
# ============================================================
set.seed(42)

# Monaco's plain "T cells" label is not a real population. It is the set of cells
# whose subtype could not be resolved from RNA. Surface protein resolves them,
# because the markers immunologists use to define T cell subsets are proteins.
DefaultAssay(pbmc) <- "ADT"

# D. Lineage panel: the minimum set needed to assign a cell to a major lineage.
# Note CD3 is decisive for T cells. A CD3-negative cell is not a T cell no matter
# how much CD4 it has, because monocytes also carry CD4 protein.
# CD14 and CD16 together split monocytes into classical, intermediate and
# non-classical.
lineage_features <- FeaturePlot(pbmc,
                                c("CD3-PROT", "CD8-PROT", "CD56-PROT", "CD45RA-PROT",
                                  "CD14-PROT", "CD16-PROT", "CD19-PROT", "HLA-DR-PROT"),
                                ncol = 4, min.cutoff = "q05", max.cutoff = "q95",
                                cols = c("grey90", "#08519C")) &
  umap_theme & NoAxes() & coord_fixed() &
  theme(plot.title = element_text(size = 9))

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

subset_adt <- DotPlot(pbmc, features = adt_subset, dot.scale = 4) +
  coord_flip() +
  labs(x = "Protein marker", y = "Cluster") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))

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

# E. Canonical PBMC marker dot plot
dot_canonical <- DotPlot(pbmc, features = markers_rna) + RotatedAxis()

# Doublet and quality checks.
# A cluster showing markers of two mutually exclusive lineages (CD3 and CD14
# together, or both CD4-high and CD8-high) is probably two cells captured in one
# droplet. Elevated counts and genes relative to neighbouring clusters supports
# that, since two cells contribute roughly twice the RNA.
v1 <- VlnPlot(pbmc, "nCount_RNA", group.by = "leiden_res.0.7", pt.size = 0) +
  NoLegend() + qc_theme + no_x

v2 <- VlnPlot(pbmc, "nFeature_RNA", group.by = "leiden_res.0.7", pt.size = 0) +
  NoLegend() + qc_theme + no_x

v3 <- VlnPlot(pbmc, "percent.mt", group.by = "leiden_res.0.7", pt.size = 0) +
  NoLegend() + qc_theme +
  xlab("Leiden cluster (res 0.7)") +
  theme(axis.title.x = element_text(size = 10))

# Unbiased marker check
# For each cluster, test every gene for higher expression in that cluster 
# than in all other cells combined (one-vs-rest, Wilcoxon by default).
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
                                 group.by = "leiden_res.0.7",
                                 size = 2.5, angle = 0, hjust = 0.5) +
  theme(axis.text.y = element_text(size = 7)) +
  NoLegend()



# ============================================================
# 9. Final labels
# ============================================================
# Labels assigned from three pieces of evidence: SingleR/Monaco, ADT surface
# phenotype, and unbiased RNA markers. Where they disagreed, protein and
# canonical markers took precedence over the reference call.
#
# This table is the single source of truth: the cluster -> label mapping applied
# below is derived from it, so the supplementary table cannot drift out of step
# with the labels used in the figures.
evidence <- tibble::tribble(
  ~cluster, ~protein, ~rna, ~label,
  "1",  "CD3+ CD4+ CD62L-hi CD197-hi CD45RA+ CD127+",
  "CCR7, SELL, TCF7, LDHB, PIK3IP1, NOSIP",
  "CD4 naive T",
  "2",  "CD3+ CD4+ CD45RO-hi CD45RA-lo CD25-hi CD127+",
  "IL7R, LTB, IL32, CD69, ITGB1; CCR7/SELL absent",
  "CD4 memory T",
  "3",  "CD14-hi CD16-neg HLA-DR+ CD11c+",
  "CD14, LYZ, S100A8/9/12, VCAN, LGALS2, MS4A6A",
  "Classical monocytes",
  "4",  "CD3+ CD8+ CD57-hi CD62L-neg CD45RO+",
  "GZMH, GNLY, FGFBP2, NKG7, CST7, GZMA/GZMM",
  "CD8 TEMRA",
  "5",  "CD3+ CD8+ CD161-hi CD45RO-hi CD127+ CD279+",
  "GZMK, KLRB1, DUSP2, IL7R, CCL5",
  "MAIT / CD8 EM",
  "6",  "CD3-neg CD56-hi CD16-hi",
  "KLRF1, KLRD1, SPON2, CLIC3, PRF1, GZMB, GNLY",
  "NK",
  "7",  "CD19+ CD20-hi IgD-hi IgM-hi CD27-neg",
  "TCL1A, MS4A1, CD79A/B, HLA-DQ/DR",
  "Naive B",
  "8",  "CD3+ CD8-hi CD62L+ CD197-hi CD45RA+",
  "CD8B, CCR7, TCF7, SELL, NOSIP, PIK3IP1, LDHB",
  "CD8 naive T",
  "9",  "No lineage-defining protein",
  "MT-* genes and MALAT1 only; several with pct.2 > pct.1",
  "Low quality",
  "10", "CD19+ CD20-hi CD27+ IgM+ IgD-lo CD11c+",
  "BANK1, MS4A1, CD79A/B, IGJ; TCL1A absent",
  "Memory B",
  "11", "CD16-hi CD14-lo CD11c-hi",
  "FCGR3A, MS4A7, CDKN1C, CSF1R, LST1, LILRB2, TCF7L2",
  "Non-classical monocytes",
  "12", "CD1c-hi CD11c+ HLA-DR-hi CD14-neg",
  "FCER1A, CLEC10A, CD1C, CPVL, HLA-DPA1/DQA1",
  "cDC2",
  "13", "CD14-hi CD16-neg CD11c+",
  "S100A8, FCGR1A, FOLR3, GBP1, WARS, TNFSF10, TYMP",
  "Activated classical monocytes",
  "14", "CD303-hi CD123-hi CD11c-neg HLA-DR+",
  "LILRA4, CLEC4C, SCT, SERPINF1, DNASE1L3, LRRC26",
  "pDC",
  "15", "CD3+ CD4+ CD8+ CD19+ CD16+ (mutually exclusive)",
  "MKI67, TYMS, RRM2, TK1, PCNA, STMN1",
  "Doublets",
  "16", "All lineage proteins flat",
  "SDPR/CAVIN2, PPBP, PF4, HIST1H2AC, TSC22D1",
  "Platelets")

# Display order: related populations next to each other, artefacts last.
# Without an explicit factor order, plots and tables sort alphabetically.
celltype_levels <- c(
  "CD4 naive T", "CD4 memory T",
  "CD8 naive T", "CD8 TEMRA", "MAIT / CD8 EM",
  "NK",
  "Naive B", "Memory B",
  "Classical monocytes", "Activated classical monocytes", "Non-classical monocytes",
  "cDC2", "pDC",
  "Platelets", "Doublets", "Low quality")

# Guards: every cluster in the data is described, and the display order names
# exactly the labels used. Either mismatch would silently produce NAs.
stopifnot(setequal(evidence$cluster, levels(factor(pbmc$leiden_res.0.7))),
          setequal(evidence$label, celltype_levels))

final_label <- setNames(evidence$label, evidence$cluster)
lab <- unname(final_label[as.character(pbmc$leiden_res.0.7)])
names(lab) <- colnames(pbmc)
pbmc$celltype <- factor(lab, levels = celltype_levels)
Idents(pbmc) <- "celltype"

# Final annotated UMAP (its own figure, not a panel of the composite)
final_annot <- DimPlot(pbmc, label = TRUE, repel = TRUE, label.size = 5.5) +
  NoLegend() +
  theme(axis.title = element_text(size = 13),
        axis.text  = element_text(size = 11))

# Cell counts and majority reference call, so the table carries data not prose
evidence <- evidence %>%
  dplyr::mutate(n_cells = as.integer(table(pbmc$leiden_res.0.7)[cluster]),
                singler = unname(majority[cluster]))

# ============================================================
# 10. Figures
# ============================================================

# Main figures
# Clusters -> reference call -> where it fell short -> protein and RNA evidence
# -> final annotated populations
fig_annotate_step <- (unsup_clus | ref_label | gen_plot) /
  wrap_elements(lineage_features) /
  wrap_elements(dot_canonical) /
  final_annot +
  plot_layout(heights = c(4.5, 6, 9, 12)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave("Annotation_Step.png", fig_annotate_step,
       width = 15, height = 23, dpi = 300, bg = "white")

# Supplementary figures 
# Reference-based annotation QC
supp_singler <- score_heatmap / delta_distribution / check_labels +
  plot_layout(heights = c(1, 1, 0.9)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave("Singler_qc.png", supp_singler,
       width = 10, height = 16, dpi = 300, bg = "white")

# Subset-defining surface proteins
ggsave("sub_adt_dotplot.png", subset_adt,
       width = 6.2, height = 7, dpi = 600)

# Unbiased top markers per cluster
ggsave("marker_heatmap.png", top_markers_heatmap,
       width = 6.2, height = 0.12 * nrow(top_markers) + 1)

# Per-cluster QC supporting the doublet and low-quality calls
supp_qc <- v1 / v2 / v3 + plot_annotation(tag_levels = "A")

ggsave("qc_per_cluster.png", supp_qc,
       width = 11, height = 9, dpi = 300)

# ---- Supplementary table ---------------------------------------------------
evidence %>%
  dplyr::select(Cluster = cluster, `n cells` = n_cells,
                `SingleR (majority)` = singler,
                `Surface protein` = protein, `RNA markers` = rna,
                `Assigned label` = label) %>%
  gt::gt() %>%
  gt::tab_header(title = "Evidence supporting cell type assignment") %>%
  gt::tab_source_note("Clusters from Leiden clustering at resolution 0.7") %>%
  gt::gtsave("Table_annotation_features.png")

# ============================================================
# 11. Export
# ============================================================
# Remove non-cells and failed populations. Platelets are real but anucleate with
# almost no transcriptome, so they are not informative here.
pbmc.clean <- subset(pbmc,
                     idents = c("Low quality", "Doublets", "Platelets"),
                     invert = TRUE)

saveRDS(pbmc.clean, "pbmc_annotated.rds")