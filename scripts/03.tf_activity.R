# TF activity inference (decoupleR ULM, per cell)
# In:  02_annotated.rds
# Out: 03_with_tf.rds, 03_tf_by_celltype.rds, Fig 1, Fig 2

source("scripts/00_setup.R")

seurat_filt <- readRDS(file.path(data_proc, "02_annotated.rds"))

# Import CollecTRI regulon network
net <- decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE)

saveRDS(net, file.path(data_proc, "03_collectri.rds"))

# Calculate TF activities per cell (ULM)
# Note: per-cell, not per-contrast. FindAllMarkers p-values are all exactly 0,
# so sign(log2FC) * -log10(p) saturates and stops discriminating.
mat <- as.matrix(seurat_filt@assays$RNA$data) 

tf_acts <- decoupleR::run_ulm(mat = mat, 
                              net = net, 
                              .source = 'source', 
                              .target = 'target',
                              .mor = 'mor', 
                              minsize = 5)

rm(mat); gc()   # Free memory

# Store TF activities as a new assay
seurat_filt[['tfsulm']] <- tf_acts %>%
  filter(statistic == 'ulm') %>%
  tidyr::pivot_wider(id_cols = 'source',
                     names_from = 'condition',
                     values_from = 'score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

DefaultAssay(seurat_filt) <- "tfsulm"
seurat_filt <- ScaleData(seurat_filt)
seurat_filt@assays$tfsulm$data <- seurat_filt@assays$tfsulm$scale.data

# Mean TF activity per cell type
df <- t(as.matrix(seurat_filt@assays$tfsulm$data)) %>%
  as.data.frame() %>%
  mutate(celltype = Idents(seurat_filt)) %>%
  pivot_longer(cols = -celltype, names_to = "source", values_to = "score") %>%
  group_by(celltype, source) %>%
  summarise(mean = mean(score), .groups = "drop")

saveRDS(seurat_filt, file.path(data_proc, "03_with_tf.rds"))
saveRDS(df,          file.path(data_proc, "03_tf_by_celltype.rds"))

# --- Fig 1: differentially active TFs across cell types ---
tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(50) %>%
  pull(source)

top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = celltype, names_from = source, values_from = mean) %>%
  column_to_rownames('celltype') %>%
  as.matrix()

pheatmap(mat = top_acts_mat,
         color = colors.use,
         border_color = "white",
         breaks = my_breaks,
         cellwidth = 15, cellheight = 15,
         treeheight_row = 20, treeheight_col = 20,
         filename = file.path(fig_dir, "fig1_tf_activity_heatmap.pdf"),
         width = 14, height = 5)

# --- Fig 2: HLH gene expression per cell type ---
DefaultAssay(seurat_filt) <- "RNA"

dot_hlh <- DotPlot(seurat_filt, features = hlh) + RotatedAxis()

ggsave(file.path(fig_dir, "fig2_hlh_gene_expression_dotplot.pdf"), 
       dot_hlh, width = 8, height = 5)