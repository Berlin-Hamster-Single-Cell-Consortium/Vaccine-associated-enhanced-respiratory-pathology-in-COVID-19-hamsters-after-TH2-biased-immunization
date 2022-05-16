library(ggplot2)
library(stringr)
library(Seurat)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
require(gridExtra)

data_dir <- "./Alum-S_1/Alum-S_1/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
AlumS_1 = CreateSeuratObject(counts = data, project="AlumS_1", min.cells=5, min.features=1000)

data_dir <- "./Alum-S_2/Alum-S_2/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
AlumS_2 = CreateSeuratObject(counts = data, project="AlumS_2", min.cells=5, min.features=1000)

data_dir <- "./Alum-S_3/Alum-S_3/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
AlumS_3 = CreateSeuratObject(counts = data, project="AlumS_3", min.cells=5, min.features=1000)

data_dir <- "./Alum-S_4/Alum-S_4/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
AlumS_4 = CreateSeuratObject(counts = data, project="AlumS_4", min.cells=5, min.features=1000)

data_dir <- "./MeV-S_1/MeV-S_1/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
MeVS_1 = CreateSeuratObject(counts = data, project="MeVS_1", min.cells=5, min.features=1000)

data_dir <- "./MeV-S_2/MeV-S_2/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
MeVS_2 = CreateSeuratObject(counts = data, project="MeVS_2", min.cells=5, min.features=1000)

data_dir <- "./MeV-S_3/MeV-S_3/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
MeVS_3 = CreateSeuratObject(counts = data, project="MeVS_3", min.cells=5, min.features=1000)

data_dir <- "./MeV-S_4/MeV-S_4/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
MeVS_4 = CreateSeuratObject(counts = data, project="MeVS_4", min.cells=5, min.features=1000)

data_dir <- "./Mock_1/Mock_1/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
Mock_1 = CreateSeuratObject(counts = data, project="Mock_1", min.cells=5, min.features=1000)

data_dir <- "./Mock_2/Mock_2/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
Mock_2 = CreateSeuratObject(counts = data, project="Mock_2", min.cells=5, min.features=1000)

data_dir <- "./Mock_3/Mock_3/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
Mock_3 = CreateSeuratObject(counts = data, project="Mock_3", min.cells=5, min.features=1000)

data_dir <- "./Mock_4/Mock_4/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
Mock_4 = CreateSeuratObject(counts = data, project="Mock_4", min.cells=5, min.features=1000)



hamster_all <- merge(AlumS_1, y = c(AlumS_2, AlumS_3, AlumS_4, MeVS_1, MeVS_2, MeVS_3, MeVS_4, Mock_1, Mock_2, Mock_3, Mock_4), add.cell.ids = c("AlumS_1", "AlumS_2", "AlumS_3", "AlumS_4", "MeVS_1", "MeVS_2", "MeVS_3", "MeVS_4", "Mock_1", "Mock_2", "Mock_3", "Mock_4"), project = "HaMM")
saveRDS(hamster_all, "./HaMM_combined.rds")

DefaultAssay(hamster_all) <- "RNA"


## integrate by hamster to remove batch effects

hamster.list <- SplitObject(hamster_all, split.by = "orig.ident")
for (i in 1:length(hamster.list)){
  hamster.list[[i]] <- SCTransform(hamster.list[[i]], verbose =T)
}

hamster.features <- SelectIntegrationFeatures(object.list = hamster.list, nfeatures = 3000)
hamster.list <- PrepSCTIntegration(object.list = hamster.list, anchor.features = hamster.features, 
                                   verbose = T)

hamster.anchors <- FindIntegrationAnchors(object.list = hamster.list, normalization.method = "SCT", 
                                          anchor.features = hamster.features, verbose = T)
hamster.integrated <- IntegrateData(anchorset = hamster.anchors, normalization.method = "SCT", 
                                    verbose = T)

## run dimensional reductions
#   PCA
hamster.integrated<- RunPCA(hamster.integrated, verbose = FALSE)
#   UMAP
hamster.integrated<- RunUMAP(hamster.integrated, dims = 1:30, verbose = FALSE)

hamster.integrated <- FindNeighbors(hamster.integrated, dims = 1:30)
hamster.integrated <- FindClusters(hamster.integrated, resolution = 0.5)


saveRDS(hamster.integrated, "./HaMM_combined_integrated.rds")

