# Cell type annotation: SingleR + manual γδ T refinement
# SingleR annotation; compare two references
expr <- GetAssayData(seurat, layer = "data")

ref_dice <- celldex::DatabaseImmuneCellExpressionData(cell.ont = "nonna")
ref_monaco <- celldex::MonacoImmuneData()

pred_dice <- SingleR(test = expr, ref = ref_dice, labels = ref_dice$label.main)
pred_monaco <- SingleR(test = expr, ref = ref_monaco, labels = ref_monaco$label.main)

seurat$celltype_dice <- pred_dice$labels
seurat$celltype_monaco <- pred_monaco$labels

# Monaco resolves DCs that DICE absorbs into monocytes
table(seurat$celltype_dice, seurat$celltype_monaco)
table(seurat$RNA_snn_res.0.3, seurat$celltype_monaco)

DimPlot(seurat, group.by = "celltype_dice", label = TRUE,
        label.size = 3, repel = TRUE) + NoLegend()

DimPlot(seurat, group.by = "celltype_monaco", label = TRUE, 
        label.size = 3, repel = TRUE) + NoLegend()

# Manual check: cluster 5 = γδ T cells (TRGC2 log2FC 5.4; PRF1, GZMA, NKG7)
# Monaco leaves these as generic "T cells" -- γδ T are neither CD4 nor CD8
FindMarkers(seurat, ident.1 = 5, only.pos = TRUE) %>% head(15)

seurat$celltype_final <- as.character(seurat$celltype_monaco)
seurat$celltype_final[seurat$RNA_snn_res.0.3 == "5"] <- "gdT cells"

Idents(seurat) <- "celltype_final"

# Keep cell populations with sufficient cells for TF inference
# Dropped: Basophils (n=1), Progenitors (n=43), unassigned "T cells" (n=195)
keep <- c("Monocytes", "CD4+ T cells", "CD8+ T cells", "gdT cells",
          "B cells", "NK cells", "Dendritic cells")
seurat_filt <- subset(seurat, idents = keep)

table(Idents(seurat_filt))

umap_ct <- DimPlot(seurat_filt, label = TRUE, repel = TRUE) + NoLegend()
ggsave("umap_celltypes.pdf", umap_ct, width = 7, height = 6)

# Markers per cell type (evidence for annotation)
markers <- FindAllMarkers(seurat_filt, only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.25)

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)