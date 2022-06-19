#MeV/Alum hamster scRNAseq (PEI/Michael Mühlebach)

library(plyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)
library(Seurat)
library(HDF5Array)
library(rhdf5)
library(hdf5r)
library(DelayedArray)
library(DelayedMatrixStats)
require(gridExtra)
library(tidyr)
library("dendextend")
library(DESeq2)
library(pheatmap)
library(cowplot)
library(SeuratObject)

##################
#read integrated object from cluster

HaMM.int <- readRDS("Path/To/Integrated_Preprocessed_Object/")

#to make a more fine-grained clustering:
DefaultAssay(HaMM.int) <- 'integrated'
HaMM.int <- FindClusters(HaMM.int, resolution = 0.6)
UMAPPlot(HaMM.int, label=TRUE)
#remove unused columns from metadata object
HaMM.int@meta.data$integrated_snn_res.0.5 <- NULL

#save/read working object
saveRDS(HaMM.int, "/Applications/stuff/HaMM/HaMM.int.rds")
HaMM.int <- readRDS("/Applications/stuff/HaMM/HaMM.int.rds")

setwd("~/Documents/Articles/Corona_HaMM/")

##############
DefaultAssay(HaMM.int) <- 'RNA'
SCoV2_rawcounts <- FetchData(HaMM.int, grep("SCoV2", HaMM.int@assays$RNA@counts@Dimnames[[1]], value="TRUE"), slot="counts")
HaMM.int@meta.data$SCoV2_sum <- rowSums(SCoV2_rawcounts)
SCoV2_rawcounts$SCoV2sum <- rowSums(SCoV2_rawcounts)
HaMM.int@meta.data$SCoV2_load <- SCoV2_rawcounts$SCoV2sum/HaMM.int@meta.data$nCount_RNA*100
HaMM.int@meta.data = cbind(HaMM.int@meta.data, HaMM.int@reductions$umap@cell.embeddings)
DefaultAssay(HaMM.int) <- 'SCT'

HaMM.int@meta.data$treatment <- gsub("_[1-4].*","",HaMM.int@meta.data$orig.ident )
HaMM.int@meta.data$hamster <- gsub(".*_([1-4].*)","Ha\\1",HaMM.int@meta.data$orig.ident)

ggplot()+geom_point(data=subset(HaMM.int@meta.data, SCoV2_load<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", size=0.5, shape = 16)+
  geom_point(data=subset(HaMM.int@meta.data, SCoV2_load>0), aes(x=UMAP_1, y=UMAP_2, colour=log10(SCoV2_load)), size=0.5, shape = 16)+
  scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("log10 SCoV2 load")
ggsave(paste("log10SCoV2_load", "pdf", sep="."), useDingbats=FALSE)

#Plot with clusters
UMAPPlot(HaMM.int, label=TRUE)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("Clusters")
ggsave(paste("clusters", "pdf", sep="."), useDingbats=FALSE)

ggplot()+
  geom_point(data=HaMM.int@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=treatment), size=0.2, shape=16, alpha=0.5)+
  ggtitle("treatment")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("treatment.pdf", useDingbats=FALSE)
ggplot()+
  geom_point(data=HaMM.int@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=hamster), size=0.2, shape=16, alpha=0.5)+
  ggtitle("hamster")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("hamsters.pdf", useDingbats=FALSE)


###############
#cell type

markers <- c("Dcn", "Tagln", "C1qb", "Rtkn2", "Cd68", "Fcgr4", "Treml4", "Mrc1", "Acta2", "Pecam1", "Foxj1", "Lamp3", "Siglecf", "Marco", "H2-Ab1", "Adgre1", "Ccr2", "Cx3cr1", "Tcf4", "Irf8", "Flt3", "Retn", "Camp", "S100a8", "Arg1", "Ccr5", "Cd79b", "Ms4a1", "Nkg7", "Gzma", "Il7r", "Cd4", "Cd3e", "ENSMAUG00000000153", "Rtkn2")
for (gene in markers) {
  df <- FetchData(HaMM.int, gene)
  gene <- str_replace_all(gene, "-", "_")
  colnames(df) <- gene
  df = cbind(df, HaMM.int@reductions$umap@cell.embeddings)
  ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", shape=16, size=0.5)+
    geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), shape=16, size=0.5)+
    scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    coord_fixed(ratio=1)+
    ggtitle(gene)
  ggsave(paste(gene, "pdf", sep="."), useDingbats=FALSE)
}

avg <- AverageExpression(HaMM.int, features = markers, return.seurat = T, slot="data") 

DoHeatmap(avg, size=5, features=markers)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("cluster_heatmap.pdf")

# Annotation Cluster <-> cell type

Idents(HaMM.int) <- HaMM.int@meta.data$seurat_clusters
HaMM.int <- RenameIdents(HaMM.int, 
                         '0'='InterstitialMacrophages',
                         '1'='AlveolarMacrophages',
                         '2'='MonocyticMacrophages',
                         '3'='MonocyticMacrophages',
                         '4'='TNKcells',
                         '5'='TNKcells',
                         '6'='TNKcells',
                         '7'='Bcells',
                         '8'='TNKcells',
                         '9'='Endothelial',
                         '10'='TNKcells',
                         '11'='MyeloidDendritic',
                         '12'='Treml4+Macrophages',
                         '13'='MonocyticMacrophages',
                         '14'='Neutrophils',
                         '15'='pDC',
                         '16'='AT12',
                         '17'='TNKcells',
                         '18'='SmoothMuscle',
                         '19'='InterstitialMacrophages',
                         '20'='AlveolarMacrophages',
                         '21'='unknown',
                         '22'='Fibroblasts',
                         '23'='Treml4+Macrophages',
                         '24'='unknown',
                         '25'='unknown')
HaMM.int@meta.data$celltype <- Idents(HaMM.int)

celltypecolors = c(
  "TNKcells" = "#368F8B",
  "MonocyticMacrophages" = "#B7245C",
  "Endothelial" = "#0D3B66",
  "AlveolarMacrophages" = "#DFACC4",
  "AT12" = "#F97E44",
  "Bcells" = "#62C370",
  "Treml4+Macrophages" = "#3E2F5B",
  "Fibroblasts" = "#B2675E",
  "MyeloidDendritic" = "#4F6D7A",
  "Neutrophils" = "#0081AF",
  "SmoothMuscle" = "#644536",
  "InterstitialMacrophages" = "#B97C9D",
  "pDC" = "#7C6A0A",
  "unknown" = "#CAD2C5"
)

means <- HaMM.int@meta.data %>%
  group_by(celltype) %>%
  dplyr::summarise(mean_U1 = mean(UMAP_1), mean_U2 = mean(UMAP_2)) %>%
  filter(celltype != "unknown")

ggplot()+
  geom_point(data=HaMM.int@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=celltype), size=0.5, shape=16)+
  geom_text(data=subset(means, celltype!="Unclear"), aes(x=mean_U1, y=mean_U2, label=celltype), size=3)+
  ggtitle("celltype")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  scale_colour_manual(values = celltypecolors)
ggsave("celltypes.pdf", useDingbats=FALSE)

celltypecolorsdf <- as.data.frame(celltypecolors)
celltypecolorsdf <- celltypecolorsdf[order(row.names(celltypecolorsdf)),]
colorvector = as.vector(celltypecolorsdf)

#Make colors darker with time
expr <- list()
for (cbright in colorvector) {
  r=(col2rgb(cbright)-40)[[1]]
  r=ifelse(r<0,0,r)
  g=(col2rgb(cbright)-40)[[2]]
  g=ifelse(g<0,0,g)
  b=(col2rgb(cbright)-40)[[3]]
  b=ifelse(b<0,0,b)
  cdark = rgb(r, g, b, maxColorValue = 255)
  expr[[cbright]] <- scales::seq_gradient_pal(cbright, cdark, "Lab")(seq(0,1,length.out=3))
}
the_colors = as.vector(unlist(expr))

the_celltypes = c("AlveolarMacrophages",
                  "InterstitialMacrophages",
                  "MonocyticMacrophages",
                  "Treml4+Macrophages",
                  "Neutrophils",
                  "MyeloidDendritic",
                  "pDC",
                  "TNKcells",
                  "Bcells",
                  "AT12",
                  "Fibroblasts",
                  "Endothelial",
                  "SmoothMuscle",
                  "unknown")

##################
#some plots

#percent cell type per treatment
df1 <- cbind.data.frame(HaMM.int@meta.data$SCoV2_load,
                        HaMM.int@meta.data$celltype,
                        HaMM.int@meta.data$treatment,
                        HaMM.int@meta.data$hamster)
colnames(df1) <- c("SCoV2_load", "celltype", "treatment", "hamster")

a = df1 %>% group_by(treatment, hamster) %>% tally(name="tot") 
b = df1 %>% group_by(celltype, treatment, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('treatment', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
#write.table(c, paste("days", gene, "tsv", sep="."), sep="\t", quote=FALSE, row.names = FALSE)

tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("celltype", "treatment"))
tgc = tgc %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
legendcolors <- c("AlumS" ="gray80", "MeVS"="gray60", "Mock"="gray40")
ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")))+
  geom_bar(aes(colour=treatment), position=position_dodge(), stat="identity", size=0)+
  geom_point(data=c, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.8), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("celltypes per time point", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_colors)+
  scale_color_manual(values = legendcolors)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))

ggsave("barplot_celltypepercentage_pertimepoint.pdf", useDingbats=FALSE)

#percent virus positive by treatment
a = df1 %>% group_by(celltype, treatment, hamster) %>% tally(name="tot") 
b = df1 %>% filter(SCoV2_load>0) %>% group_by(celltype, treatment, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('celltype', 'treatment', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
#write.table(c, paste("days", gene, "tsv", sep="."), sep="\t", quote=FALSE, row.names = FALSE)
tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("celltype", "treatment"))
tgc = tgc %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
legendcolors <- c("AlumS" ="gray80", "MeVS"="gray60", "Mock"="gray40")
ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")))+
  geom_bar(aes(colour=treatment), position=position_dodge(), stat="identity", size=0)+
  geom_point(data=c, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.8), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("celltypes per time point", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_colors)+
  scale_color_manual(values = legendcolors)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))
ggsave("barplot_viruspositivepercelltype_pertimepoint.pdf", useDingbats=FALSE)



