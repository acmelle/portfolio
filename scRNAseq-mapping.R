#/Users/amanda.mitchell/opt/miniconda3/bin/R

# install pacakges #
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c('GenomicRanges', 'rhdf5', 'BiocNeighbors', 'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment','SummarizedExperiment', 'batchelor', 'Matrix.utils'))

library(remotes)
install_github("mojaveazure/seurat-disk")
install_github('satijalab/seurat-data')
install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())

library(remotes)
BiocManager::install(c("AUCell", "RcisTarget"))
BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost

## Optional (but highly recommended):
BiocManager::install(c("zoo", "mixtools", "rbokeh"))
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
BiocManager::install(c("doMC", "doRNG"))

# SCENIC install
install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
install_github("aertslab/SCENIC") 
install_github("aertslab/SCENIC@v1.1.2")


# load pacakges #
library(Seurat)
library(harmony) 
library(zellkonverter)
library(SeuratDisk)
library(SeuratData)
library(loomR)
library(caret) #do not load with ArchR
#library(ArchR) #do not load with caret
library(data.table)
library(Biobase)
library(SingleCellExperiment)
library(GEOquery)
library(SCopeLoomR)
library(SCENIC)


# all references integrated #
myclusters = "/mnt/bh1/storage/external/external/sc-datasets/integration/seurat"
tomesave = "tome_int.rds"
markersave = "mmparkers.txt"
setwd(myclusters)
#query = readRDS("/mnt/bh1/data/user/amanda.mitchell/data/RXRX/tenx/runa/Seurat/dopamine_all/tome.rds")

reference_a = readRDS("/mnt/bh1/storage/external/external/sc-datasets/neurodev/paulsen-2022-ASD/Seurat/final/neurodvtome.rds") #CellType
reference_b = readRDS("/mnt/bh1/storage/external/external/sc-datasets/AD/grubman-2019-uglia/seurat/tome.rds") #cellType
reference_c = readRDS("/mnt/bh1/storage/external/external/sc-datasets/atlases/AIBS/AIBS-M1-hs/seurat/tome.rds") #subclass_label
reference_d = readRDS("/mnt/bh1/storage/external/external/sc-datasets/AD/saddick-2022-uglia/seurat/tome.rds") #celltype
reference_e = readRDS("/mnt/bh1/storage/external/external/sc-datasets/ALS/ho-2020-iMN/seurat/tome.rds") #celltype
reference_f = readRDS("/mnt/bh1/storage/external/external/sc-datasets/AD/polioudakis-2019/seurat/tome.rds") #Cluster
reference_g = readRDS("/mnt/bh1/storage/external/external/sc-datasets/PD/kamath-2021-hSN/seurat/tome.rds") #Cell_Type
reference_h = readRDS("/mnt/bh1/storage/external/external/sc-datasets/AD/schoernig-2020-NGN2/seurat/tome.rds") #celltype
#reference_g = reference_g[, sample(colnames(reference_g), size = 200000, replace=F)]

reference_a$celltype = reference_a$CellType
reference_b$celltype = reference_b$cellType
reference_c$celltype = reference_c$subclass_label
reference_f$celltype = reference_f$Cluster
reference_g$celltype = reference_g$Cell_Type

reference_a$dataset = "Paulsen" #iPSC cortical neurons, 300,000
reference_b$dataset = "Grubman" #cortical microglia 10,000
reference_c$dataset = "AIBS" #MTG 24,000
reference_d$dataset = "Sadick" #spinal cord microglia 77,000
reference_e$dataset = "Ho" #iPSC motor neurons 5,970
reference_f$dataset = "Polioudakis" #development 2,433
reference_g$dataset = "Kamath" #midbrain 387,000
reference_h$dataset = "Schoernig" #21,4500

reference_a$type = "human iPSC cortical organoids"
reference_b$type = "human cortical microglia"
reference_c$type = "human MTG"
reference_d$type = "human spinal cord microglia"
reference_e$type = "human iPSC motor neurons"
reference_f$type = "human developmental cortex"
reference_g$type = "human midbrain"
reference_h$type = "Schoernig"

reference_a = reference_a[, sample(colnames(reference_a), size = 100000, replace=F)]
reference_g = reference_g[, sample(colnames(reference_g), size = 100000, replace=F)]

mapping.list = list(reference_a, reference_b, reference_c, reference_d, reference_e, reference_f, reference_g, reference_h)
mapping.anchors <- FindIntegrationAnchors(object.list = mapping.list, dims = 1:30, reduction = "rpca")
mapping.integrated <- IntegrateData(anchorset = mapping.anchors, dims = 1:30)
DefaultAssay(mapping.integrated) <- "integrated"

mapping.integrated <- ScaleData(mapping.integrated, verbose = FALSE)
mapping.integrated <- RunPCA(mapping.integrated, npcs = 30, verbose = FALSE)
mapping.integrated <- RunUMAP(mapping.integrated, reduction = "pca", dims = 1:30)
umap_dataset = DimPlot(mapping.integrated, reduction = "umap", group.by = "dataset", label = TRUE, repel = TRUE)
umap_celltype = DimPlot(mapping.integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
umap_type = DimPlot(mapping.integrated, reduction = "umap", group.by = "type", label = TRUE, repel = TRUE)

#GRAPH UMAP
png(file="umap_dataset.png", width = 900, height = 900, units = "px", pointsize = 12, bg = "white",  res = NA)
umap_dataset
dev.off()

png(file="umap_celltype.png", width = 2000, height = 900, units = "px", pointsize = 12, bg = "white",  res = NA)
umap_celltype
dev.off()

png(file="umap_type.png", width = 2000, height = 900, units = "px", pointsize = 12, bg = "white",  res = NA)
umap_type
dev.off()

setwd(myclusters)
saveRDS(mapping.integrated, tomesave)


# cleanup #
reference = readRDS("/mnt/bh1/storage/external/external/sc-datasets/integration/seurat/tome_int.rds")
reference[['CellType', 'treat', 'class_label', 'con', 'Subcluster', 'Cells', 'ciortical_layer_label', 'Donor', 'org', 'condition', 'cluster_label',
        'age', 'genotype', 'Phase', 'cellType', 'cluster_color', 'Sex', 'name', 'subclusterID', 'Cell_Type', 'subclass_label', 'info',
        'Cluster', 'class', 'region_label', 'Layer', 'sample', 'geno', 'seurat_clusters', 'cell.names']] = NULL
saveRDS(reference, file="/mnt/bh1/storage/external/external/sc-datasets/integration/seurat/reference.rds")
reference = readRDS("/mnt/bh1/storage/external/external/sc-datasets/integration/seurat/reference.rds")
reference_small = reference[, sample(colnames(reference), size = 50000, replace=F)]
saveRDS(reference_small, "/mnt/bh1/storage/external/external/sc-datasets/integration/seurat/reference_small.rds")


# export to python
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(loomR)
reference = readRDS("/mnt/bh1/storage/external/external/sc-datasets/integration/seurat/reference.rds")
reference_small = readRDS("/mnt/bh1/storage/external/external/sc-datasets/integration/seurat/reference_small.rds")

#Seurat to h5ad
SaveH5Seurat(reference, filename = "/mnt/bh1/storage/external/external/sc-datasets/integration/files/reference.h5Seurat")
Convert("/mnt/bh1/storage/external/external/sc-datasets/integration/files/reference.h5Seurat", dest = "h5ad")
SaveH5Seurat(reference_small, filename = "/mnt/bh1/storage/external/external/sc-datasets/integration/files/reference_small.h5Seurat")
Convert("/mnt/bh1/storage/external/external/sc-datasets/integration/files/reference_small.h5Seurat", dest = "h5ad")

#h5ad to loom
Convert("//mnt/bh1/storage/external/external/sc-datasets/integration/files/reference.h5ad", dest = "h5seurat", overwrite = TRUE)
reference  <- LoadH5Seurat("/mnt/bh1/storage/external/external/sc-datasets/integration/files/reference.h5Seurat")
SaveLoom(reference, filename = "/mnt/bh1/storage/external/external/sc-datasets/integration/files/reference.loom")
saveRDS(reference, file="/mnt/bh1/storage/external/external/sc-datasets/integration/files/reference.rds")

Convert("//mnt/bh1/storage/external/external/sc-datasets/integration/files/reference_small.h5ad", dest = "h5seurat", overwrite = TRUE)
reference_small  <- LoadH5Seurat("/mnt/bh1/storage/external/external/sc-datasets/integration/files/reference_small.h5Seurat")
SaveLoom(reference_small, filename = "/mnt/bh1/storage/external/external/sc-datasets/integration/files/reference_small.loom")
saveRDS(reference_small, file="/mnt/bh1/storage/external/external/sc-datasets/integration/files/reference_small.rds")
