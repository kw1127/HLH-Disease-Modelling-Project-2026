# Cell type annotation: SingleR + manual γδ T refinement
# SingleR annotation; compare two references
expr <- GetAssayData(seurat, layer = "data")

# to-do: save plots
ref_dice <- celldex::DatabaseImmuneCellExpressionData(cell.ont = "nonna")
ref_monaco <- celldex::MonacoImmuneData()

# cell type annotation with SingleR
# With DICE labels
pred_dice <- SingleR(
  test = expr, 
  ref = ref_dice, 
  labels = ref_dice$label.main)

# With monaco labels
pred_monaco <- SingleR(
  test = expr, 
  ref = ref_monaco, 
  labels = ref_monaco$label.main)

# Assign to seurat object
seurat$celltype_dice <- pred_dice$labels
seurat$celltype_monaco <- pred_monaco$labels

# Monaco resolves DCs that DICE absorbs into monocytes
tab <- table(seurat$celltype_dice, seurat$celltype_monaco)

# Converts the table into a formatted image.
as.data.frame.matrix(tab) |> # table -> data frame; keeps row/col names as a proper matrix
  tibble::rownames_to_column("DICE_celltype") |> # moves the DICE labels from row names into a real column
  gt() |> # builds a gt table object
  gtsave("annotation_table.png") # saves the table as an image

# Row-normalise so each DICE label sums to 1 (shows proportion, not raw counts),
# fairer when cell types have very different sizes)
tab_prop <- prop.table(tab, margin = 1)

# Heatmap of the row-normalised concordance matrix
annot_heatmap <- pheatmap(tab_prop,
              cluster_rows = FALSE, # keep rows in their given order
              cluster_cols = FALSE, # keep columns in their given order
              display_numbers = TRUE, # print the numeric value inside each cell
              number_format = "%.2f", # format those numbers to 2 decimal places
              main = "DICE vs Monaco annotation concordance", # add title
              filename = "annotation_concordance.png", # save
              width = 8, 
              height = 7)

table(seurat$RNA_snn_res.0.3, seurat$celltype_monaco)

dice <- DimPlot(seurat, group.by = "celltype_dice", label = TRUE,
        label.size = 3, repel = TRUE) + NoLegend()

monaco <- DimPlot(seurat, group.by = "celltype_monaco", label = TRUE, 
        label.size = 3, repel = TRUE) + NoLegend()

annotation_compare <- dice + monaco
ggsave("annotation_compare.png", plot = annotation_compare, width = 14, height = 6, dpi = 300)

# Monaco annotation leaves an unassigned T cell cluster
# Find positive marker genes for cluster 5 (the unassigned T cluster)
markers_c5 <- FindMarkers(seurat, ident.1 = 5, only.pos = TRUE) %>%
  filter(p_val_adj < 0.05) %>% # keep significant
  mutate(pct_diff = pct.1 - pct.2) %>% # specificity
  arrange(desc(pct_diff), desc(avg_log2FC)) %>% # specific first, then strong
  head(15) # keep the top 15

# Format the marker table and save it as an image.
markers_c5 %>%
  tibble::rownames_to_column("gene") %>% # move gene names from row names into a "gene" column
  dplyr::mutate(across(where(is.numeric), ~signif(.x, 3))) %>% # round all numeric columns to 3 significant figures
  gt() %>% # build a styled gt table
  gt::tab_header(title = "Top markers defining cluster 5") %>% # table header
  gt::cols_label(avg_log2FC = "log2FC", p_val_adj = "p.adj") %>% # column names
  gtsave("cluster5_markers.png") # save the table as an image
  gt() %>% # build a styled gt table
  gt::tab_header(title = "Top markers defining cluster 5") %>% # table header
  gt::cols_label(avg_log2FC = "log2FC", p_val_adj = "p.adj") %>% # column names
  gtsave("cluster5_markers.png") # save the table as an image

# Create a final cell-type label from the Monaco annotation
seurat$celltype_final <- as.character(seurat$celltype_monaco)

# cluster 5 reassigned to γδT cells.
# Top markers are the γδ TCR chain genes TRGV2 (log2FC 6.5), TRDV2 (6.2), and TRGV4 (5.9),
# expressed almost exclusively in this cluster (pct.2 ~ 0), the defining γδ signature, consistent with the circulating Vγ9Vδ2 subset. 
# ZNF683 (Hobit) supports an innate-like effector phenotype. 
# Monaco leaves these as generic "T cells" since γδ T are neither
# CD4 nor CD8 and don't match its reference profiles.
seurat$celltype_final[seurat$RNA_snn_res.0.3 == "5"] <- "γδT cells"

# Set the active identity of the seurat object to the final cell-type labels
Idents(seurat) <- "celltype_final"

# Keep cell populations with sufficient cells for TF inference
# Dropped: Basophils (n=1), Progenitors (n=43)
keep <- c("Monocytes", "CD4+ T cells", "CD8+ T cells", "γδT cells",
          "B cells", "NK cells", "Dendritic cells")

seurat_filt <- subset(seurat, idents = keep)

as.data.frame(table(Idents(seurat))) %>% # all cell populations with counts
  dplyr::rename(`Cell type` = Var1, `n cells` = Freq) %>%
  dplyr::mutate(Status = dplyr::case_when(
    `Cell type` %in% keep ~ "Retained",
    `Cell type` == "T cells" ~ "Dropped (unresolved identity)", 
    TRUE ~ "Dropped (n < 50)" # too few cells for stable inference
  )) %>%
  dplyr::arrange(desc(`n cells`)) %>% # largest populations first
  gt() %>%
  gt::tab_header(title = "Cell populations before TF inference") %>%
  gtsave("populations_filtering.png")

# UMAP of the cell type annotation with SingleR, monaco ref data, and manual annotation
umap_ct <- DimPlot(seurat_filt, label = TRUE, repel = TRUE) +
  NoLegend()

ggsave("umap_celltypes.pdf", umap_ct, width = 7, height = 6)

# Markers per cell type (validation for annotation).
# Find features differentially expressed in each cell type vs all other cells, 
# keeping only upregulated genes. 
# min.pct = 0.25 skips genes detected in <25% of cells in either group. 
# logfc.threshold = 0.25 skips weak-effect genes
all_markers <- FindAllMarkers(seurat_filt, 
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)

# what are the top 5 differentially-expressed genes per cell type?
top_markers <- all_markers %>%
  dplyr::filter(p_val_adj < 0.05, pct.1 > 0.25) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(n = 5, order_by = avg_log2FC, with_ties = FALSE) %>%
  dplyr::pull(gene) %>%
  unique()

# Dot plot of the top 5 markers genes per cell type
# Plots the z-scored expression of each gene in each cell type
DotPlot(seurat_filt, features = top_markers) +
  ggplot2::coord_flip() + # genes on y, cell types on x 
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

ggsave("celltype_dotplot.png", width = 10, height = 12, dpi = 300)

# Heatmap of top 10 differentially-expressed genes per cell type
top_markers_hm <- all_markers %>%
  dplyr::filter(p_val_adj < 0.05, pct.1 > 0.25) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(n = 10, order_by = avg_log2FC, with_ties = FALSE)

DoHeatmap(seurat_filt,
          features = top_markers_hm$gene,
          angle = 45,
          size = 4,
          hjust = 0) +
  ggplot2::ggtitle("Top 10 marker genes per cell type") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold")) +
  NoLegend()

ggsave("celltype_heatmap.png", width = 12, height = 14, dpi = 300)
