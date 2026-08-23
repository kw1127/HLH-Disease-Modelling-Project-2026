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
# connected to each other than to the rest of the graph, and unlike Louvain, it
# guarantees every cluster is internally connected through refinement.
# n.iter = 10 repeats the refinement to reach a stable partition.
pbmc <- FindClusters(pbmc, algorithm = 4, n.iter = 10, resolution = res_seq, random.seed = 42)

# Rename RNA_snn_res.* -> leiden_res.* for clarity
snn_cols <- grep("^RNA_snn_res\\.", colnames(pbmc@meta.data), value = TRUE)
leiden_cols <- sub("^RNA_snn", "leiden", snn_cols)
pbmc@meta.data[leiden_cols] <- pbmc@meta.data[snn_cols]
pbmc@meta.data[snn_cols] <- NULL

# Compute a 2D UMAP embedding from the same 15 PCs for visualization.
# Nonlinear projection that places similar cells near each other in a 2D space, capturing local structure.
# Use the same dims as FindNeighbors so the graph and the plot are consistent.
pbmc <- RunUMAP(pbmc, dims = 1:15, seed.use = 42)

# Find how cell clusters separate with each resolution
ggsave("06_clustree.png", clustree(pbmc, prefix = "leiden_res."), width = 10, height = 10, dpi = 300)