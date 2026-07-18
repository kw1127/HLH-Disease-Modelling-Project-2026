# Install packages ------------------------------------------------------------
# BiocManager coordinates both CRAN and Bioconductor installs
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_pkgs <- c(
  "dplyr", "Seurat", "patchwork", "hdf5r", "clustree",
  "pheatmap", "tibble", "tidyr", "gt"
)
bioc_pkgs <- c("celldex", "SingleR", "decoupleR", "OmnipathR")

# BiocManager::install() pulls CRAN packages too, so install everything together
to_install <- c(cran_pkgs, bioc_pkgs)
to_install <- to_install[!vapply(to_install, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) {
  BiocManager::install(to_install)
}

# Load packages ---------------------------------------------------------------
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
library(gt)
library(webshot2)
