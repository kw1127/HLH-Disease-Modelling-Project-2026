# ============================================================
# 5. Clustering
# ============================================================

# Build a shared nearest-neighbor (SNN) graph from the first 15 PCs.
# For each cell, find its k nearest neighbors in PC space, then weight the edges between two cells by how many
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

clustering_tbl <- comparison |>
  gt() |>
  cols_label(resolution = "Resolution",
             n_leiden = "Leiden clusters",
             n_louvain = "Louvain clusters",
             ARI = "ARI") |>
  fmt_number(columns = ARI, decimals = 3) |>
  fmt_number(columns = resolution, decimals = 1) |>
  cols_align(align = "center", columns = everything()) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = resolution == 1.0)) |>
  tab_footnote(
    footnote = "Adjusted Rand index (ARI) between Leiden and Louvain clustering at the same resolution.",
    locations = cells_column_labels(columns = ARI)) |>
  tab_options(table.font.size = 11,
              data_row.padding = px(4))

# Export as both a png and editable docx
gtsave(clustering_tbl, "table_ari_leiden_louvain.png", vwidth = 700, expand = 5)
gtsave(clustering_tbl, "table_ari_leiden_louvain.docx")
gtsave(clustering_tbl, "table1.html")

leiden_res <- DimPlot(pbmc, group.by = paste0("leiden_res.", res_seq),
                      label = TRUE, label.size = 3, ncol = 4) &
  NoLegend() & 
  NoAxes() &
  theme(plot.title = element_text(size = 10))

louvain_res <- DimPlot(pbmc, group.by = paste0("louvain_res.", res_seq),
                       label = TRUE, label.size = 3, ncol = 4) &
  NoLegend() & 
  NoAxes() &
  theme(plot.title = element_text(size = 10))

sup1 <- wrap_elements(leiden_res) / wrap_elements(louvain_res) +
  plot_annotation(tag_levels = "A")

ggsave("supp_fig_01_dimred.png", sup1, width = 14, height = 14, dpi = 300)

# Find how cell clusters separate with each resolution
ggsave("06_clustree.png", clustree(pbmc, prefix = "leiden_res."), width = 10, height = 10, dpi = 300)

Idents(pbmc) <- "leiden_res.0.7"
ggsave("07_umap_res07.png",
       DimPlot(pbmc, reduction = "umap", group.by = "leiden_res.0.7", label = TRUE) + 
         NoLegend() +
         ggtitle("Leiden, resolution 0.7"),
       width = 7, height = 6, dpi = 300)

# Unsupervised clusters
leiden_res07 <- DimPlot(pbmc, reduction = "umap", group.by = "leiden_res.0.7",
                        label = TRUE, repel = TRUE, label.size = 3) +
  guides(colour = "none") +
  coord_fixed() +
  NoAxes()