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