source("scripts/00_setup.R")

# Import PBMC data
counts <- Read10X_h5("~/Documents/HLH Disease Modelling Project 2026/SC3_v3_NextGem_SI_PBMC_10K_filtered_feature_bc_matrix.h5")

# Create Seurat object
seurat <- CreateSeuratObject(counts = counts,
                             project = "pbmc10k",
                             min.cells = 3,
                             min.features = 200)

# Quality control
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

seurat <- subset(seurat, subset = nFeature_RNA > 200 & 
                   nFeature_RNA < 6000 & 
                   percent.mt < 10)

# Normalise, scale, and PCA
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, nfeatures = 2000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)

ElbowPlot(seurat)                     

# Cluster & UMAP 
seurat <- FindNeighbors(seurat, dims = 1:15)
seurat <- FindClusters(seurat, resolution = c(0.3, 0.5, 0.8, 1.0))
seurat <- RunUMAP(seurat, dims = 1:15)

clustree(seurat) # Justify which resolution to use
Idents(seurat) <- "RNA_snn_res.0.3"  # Set resolution