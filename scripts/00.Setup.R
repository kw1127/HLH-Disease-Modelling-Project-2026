# Load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(hdf5r)
library(clustree)
library(celldex)
library(SingleR)
library(decoupleR)
library(pheatmap)
library(tibble)
library(tidyr)
library(OmnipathR)

data_raw  <- "data/raw"
data_proc <- "data/processed"
fig_dir   <- "outputs/figures"
carn_dir  <- "outputs/carnival"

colors     <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
colors.use <- grDevices::colorRampPalette(colors)(100)
my_breaks  <- c(seq(-2, 0, length.out = 51), seq(0.05, 2, length.out = 50))

hlh <- c("PRF1", "UNC13D", "STX11", "STXBP2",
         "LYST", "SH2D1A", "XIAP", "RAB27A")

celltypes <- c("NK cells", "gdT cells", "CD8+ T cells")

expressed <- function(obj, ct, min_pct = 0.1) {
  cells <- WhichCells(obj, idents = ct)
  cnt   <- GetAssayData(obj, assay = "RNA", layer = "counts")[, cells]
  rownames(cnt)[Matrix::rowMeans(cnt > 0) >= min_pct]
}