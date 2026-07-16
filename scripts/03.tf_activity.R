# Import the CollecTRI regulon network
# split_complexes = FALSE keeps regulating complexes as whole
# Setting it to TRUE might help pinpoint specific TF-target interactions however
net <- decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE)

# Extract the scaled GEX counts as a matrix
mat <- as.matrix(seurat_filt@assays$RNA$data) 

# Per-cell transcription factor (TF) activity inference using a univariate linear model (ULM).
#
# Why ULM on the expression matrix rather than ranking by DE stats:
# FindAllMarkers reports p_val_adj of exactly 0 because with thousands of cells the tests are
# massively overpowered. That breaks any sign(log2FC) * -log10(p) gene-ranking statistic, 
# since -log10(0) is infinite and every gene saturates to the same value, losing all discrimination. 
# Running ULM on the expression matrix directly sidesteps this: it uses the continuous expression 
# values against the regulon rather than a saturated p-value transform, so it stays discriminative.
#
# For each cell and TF, ULM fits a simple linear regression where
# the cell's GEX values (response) are modelled against that TF's known regulon (target genes)
# each with a "mode of regulation" weight (+1 activation, -1 repression). 
# The fitted slope (t-value) becomes the TF activity score.
# A high positive score means the TF's activating targets are up and its repressed targets
# are down in that cell, implying the TF is active. 
# Because it's run per cell, you get a cells x TFs activity matrix that can then be summarised per cell type.
# ULM is univariate: each TF is tested independently (one regression per TF), which is fast
# and robust, though it doesn't account for TFs sharing target genes.
tf_acts <- decoupleR::run_ulm(mat = mat, # # the scaled expression matrix (genes x cells)
                              net = net, # prior TF-target network (regulons) from collecTRI
                              .source = 'source', # column in `net` naming the TF
                              .target = 'target', # column in `net` naming each target gene
                              .mor = 'mor', # column in `net` giving mode of regulation / weight
                              minsize = 5) # ignore TFs with fewer than 5 measured targets

# Compute BH-adjusted p-values on the ULM results
tf_acts <- tf_acts %>%
  dplyr::filter(statistic == 'ulm') %>% # keep ULM rows only
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) # FDR across all TF-cell tests

# Store the per-cell TF activity scores as a new assay
seurat_filt[['tfsulm']] <- tf_acts %>%
  tidyr::pivot_wider(id_cols = 'source', # reshape to wide: one row per TF
                     names_from = 'condition', # one column per cell
                     values_from = 'score') %>% # fill with the ULM activity score
  tibble::column_to_rownames('source') %>% # move TF names into row names (assay features = TFs)
  Seurat::CreateAssayObject(.) # build a Seurat assay (TFs x cells) from the wide matrix

DefaultAssay(seurat_filt) <- "tfsulm"
seurat_filt <- ScaleData(seurat_filt)
seurat_filt@assays$tfsulm$data <- seurat_filt@assays$tfsulm$scale.data

# Mean TF activity per immune cell type of interest
celltypes <- c("NK cells", 
               "gdT cells", 
               "CD8+ T cells")

# Collapse the per-cell TF activity scores from ULM to one mean activity value per cell type-TF pair
df <- t(as.matrix(seurat_filt@assays$tfsulm$data)) %>% # # TFs x cells -> cells x TFs (cells as rows)
  as.data.frame() %>% # convert to data frame
  mutate(celltype = Idents(seurat_filt)) %>% # tag each cell with its cell-type label
  pivot_longer(cols = -celltype, names_to = "source", values_to = "score") %>% # wide -> long: one row per cell-TF pair
  group_by(celltype, source) %>% # group by cell type and TF
  summarise(mean = mean(score), .groups = "drop") # averages the scaled activity across all cells of that type

# How many TFs to select for analysis?

# Select the most differentially active TFs across cell types.
# Why SD? For each TF, mean = activity in each cell type. The SD of means
# across cell types measures how much a TF's activity varies between cell types.
# High SD = active in some types but not others.
# Low SD = roughly constant everywhere. 
# Ranking by SD shows the TFs that best separate cell types.
sd_ranked <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(desc(std)) %>%
  mutate(rank = row_number())

library(ggplot2)
ggplot(sd_ranked, aes(x = rank, y = std)) +
  geom_line() +
  geom_point(size = 0.6) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = "red") +
  labs(x = "TF rank by between-cell-type SD",
       y = "SD of mean activity across cell types",
       title = "Selection of differentially active TFs based on their standard deviation (SD)") +
  theme_minimal()

ggsave("tf_elbow.png", width = 7, height = 5, dpi = 300)

# Differentially active TFs across cell types
top_tfs <- sd_ranked %>%
  slice_head(n = 50) %>% # top 50 by SD based off elbow plot
  pull(source) # extract just the TF names

top_acts_mat <- df %>%
  filter(source %in% top_tfs) %>%
  pivot_wider(id_cols = celltype, names_from = source, values_from = mean) %>%
  column_to_rownames('celltype') %>%
  as.matrix()

colors <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
colors.use <- grDevices::colorRampPalette(colors)(100)

my_breaks  <- c(seq(-2, 0, length.out = 51), seq(0.05, 2, length.out = 50))

pheatmap(mat = top_acts_mat,
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
         legend_breaks = c(-2, -1, 0, 1, 2),
         legend_labels = c("-2\n(lower\nactivity)", "-1", "0", "1", "2\n(higher\nactivity)"),  # annotate what colour means
         filename = "tf_activity_heatmap.pdf",
         width = 16, height = 5)           

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

DefaultAssay(seurat_filt) <- "RNA"

# Plot mean expression of HLH-associated genes across the cell types
# Separate the fHL-associated genes from the other primary HLH genes
dot_hlh <- DotPlot(seurat_filt, features = levels(hlh)) +
  RotatedAxis() +
  labs(title = "Expression of primary HLH-associated genes across cell types",
       x = "HLH-associated gene",
       y = "Cell type") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10,), 
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 9)) + 
  geom_vline(xintercept = c(4.5, 6.5), linetype = "dashed", colour = "grey70")

ggsave("hlh_gene_expression_dotplot.png", dot_hlh, width = 8, height = 5)

# Mean TF activity per cell type of top 50 TFs
DefaultAssay(seurat_filt) <- "tfsulm"

dot_tfs <- DotPlot(seurat_filt, features = top_tfs) +
  RotatedAxis() +
  labs(title = "Mean activity of the top 50 TFs across cell types",
       x = "Transcription factor",
       y = "Cell type") +
  theme(axis.text.x = element_text(size = 7)) +
  scale_colour_gradient2(
    low = "#2166AC",      # repressed (negative activity) 
    mid = "white",        # zero - no activity
    high = "#B2182B",     # active (positive activity) 
    midpoint = 0,
    name = "Mean activity"
  )       

ggsave("tf_activity_dotplot.png", dot_tfs, width = 16, height = 5)

# Top 5 TFs per cell type by mean activity
top_per_type <- df %>%
  group_by(celltype) %>%
  slice_max(mean, n = 5, with_ties = FALSE) %>% 
  pull(source) %>%
  unique()

DefaultAssay(seurat_filt) <- "tfsulm"

dot_tfs_focused <- DotPlot(seurat_filt, features = top_per_type) +
  RotatedAxis() +
  scale_colour_gradient2(
    low = "#2166AC", 
    mid = "white", 
    high = "#B2182B",
    midpoint = 0, 
    name = "Mean activity") +
  labs(title = "Mean activity of the top cell-type-defining TFs",
       x = "Transcription factor",
       y = "Cell type") +
  theme(plot.title  = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 10))

ggsave("tf_activity_dotplot_focused.png", dot_tfs_focused, width = 10, height = 5)
