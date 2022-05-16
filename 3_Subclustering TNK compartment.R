

library(Seurat)
library(ggplot2)
library(tidyverse)
library(biomaRt)

Idents(sars2.ham) <- "TNKcells"
sub <- subset(sars2.ham, idents = "TNKcells")


## read in data

HaMM.int <- readRDS("Path/To/Integrated_Preprocessed_Object/")

Idents(HaMM.int) <- "celltype"

## subset

tnk <- subset(HaMM.int, subset = celltype == "TNKcells")

DefaultAssay(tnk) <- "RNA"
tnk <- DietSeurat(tnk, assays = c("RNA"))

## integrate data


hamster.list <- SplitObject(tnk, split.by = "orig.ident")
for (i in 1:length(hamster.list)){
  hamster.list[[i]] <- SCTransform(hamster.list[[i]], verbose =T)
}

hamster.features <- SelectIntegrationFeatures(object.list = hamster.list, nfeatures = 3000)
hamster.list <- PrepSCTIntegration(object.list = hamster.list, anchor.features = hamster.features, 
                                   verbose = T)

hamster.anchors <- FindIntegrationAnchors(object.list = hamster.list, normalization.method = "SCT", 
                                          anchor.features = hamster.features, verbose = T)
tnk.int <- IntegrateData(anchorset = hamster.anchors, normalization.method = "SCT", 
                                    verbose = T)


## perform dimensional reduction

tnk.int<- RunPCA(tnk.int, verbose = FALSE)
tnk.int<- RunUMAP(tnk.int, dims = 1:30, verbose = FALSE)

## cluster
tnk.int <- FindNeighbors(tnk.int, dims = 1:30)
tnk.int <- FindClusters(tnk.int, resolution = 0.8)
tnk.int <- FindClusters(tnk.int, resolution = 0.6)
tnk.int <- FindClusters(tnk.int, resolution = 0.4)


## norm and scale RNA assay
DefaultAssay(tnk.int) <- "RNA"
tnk.int <- NormalizeData(tnk.int)
all.genes <- rownames(tnk.int)
tnk.int <- ScaleData(tnk.int, features = all.genes)

## check expression of marker genes for each cluster
Idents(tnk) <- "integrated_snn_res.0.6"

markers <- c("ENSMAUG00000000153", "Cd3e", "Cd3d", "S100a4", "Ccr7", "Nkg7", "Cxcr5", "Pdcd1", "Cd4", "Gzmb", "Gzma", "Prf1", "Ikzf2", "Foxp3", "Il7r", "Cxcr6", "Slamf1", "Rorc", "Crem")


avg <- AverageExpression(tnk, features = markers, return.seurat = T, assays = "RNA") 
avg

h <- DoHeatmap(avg, assay="RNA", size=5, features=markers, label=T, group.bar=F, group.bar.height = 0,
               draw.lines = F, hjust = 0)+ 
  scale_fill_gradientn(colors = c("blue", "white","red"), na.value = "white") +
  NoLegend()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0))

## assign celltype from excel table based on clustering at resolution 0.6 

ct_assign <- read.csv("assign_celltype_tnksub.csv", sep = ";")

# create named vector
ct_assign <- deframe(ct_assign)

# assign IDs based on named vector
tnk$subtype <- tnk$integrated_snn_res.0.6
tnk$subtype<- plyr::revalue(tnk$subtype, replace = ct_assign)


