## generate pathway enrichment results using KEGG and GO:BP databases

## generate DEGs for all celltypes vs mock

library(Seurat)
library(GSVA)
library(msigdbr)
library(tidyverse)
library(pheatmap)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ggpubr)
library(svglite)

#----------------------------------

hamster <- readRDS("path/to/preprocessed/object")

DefaultAssay(hamster) <- "RNA"
hamster <- DietSeurat(hamster, scale.data = T, assays = c("RNA"))

#create celltype treat
hamster$celltype_treat <- paste0(hamster$celltype, "_", hamster$treatment)
Idents(hamster) <- "celltype_treat"
head(hamster@meta.data)

#get human-to-mouse conversion of gene IDs

# human data base
ensembl_human <-  useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# convert to enterz id
hamster_to_human <- getBM(attributes=c("hgnc_symbol", "external_gene_name","entrezgene_accession",
                                       "entrezgene_id"),
                          filters = 'hgnc_symbol',
                          values = rownames(hamster), ###
                          mart = ensembl_human)


#### IDENTIFY DEGs BETWEEN ALUMS and MEVS VACCINATED vs MOCK, RUN PATHWAY ENRICHMENT ####

## create DEGs by celltype
## comparing AlumS vaccinated and MeVS vaccinated to Mock vaccinated
## run pathway enrichment using ClusterProfiler gseKEGG and gseGO
## visualise in bubbleplots

ct_list <- unique(hamster$celltype)

for (ct in ct_list){
  
  
  ####################### AlumS ############################
  
  print(ct)
  
  ident1 <- paste0(ct, "_AlumS")
  ident2 <- paste0(ct, "_Mock")
  
  # get DEGs
  as_vs_mock <- FindMarkers(hamster, ident.1 = ident1, ident.2 = ident2, logfc.threshold = 0 )
  as_vs_mock$gene <- rownames(as_vs_mock)
  # convert to human equivalents
  as_vs_mock["Name2"] <- toupper(as_vs_mock$gene)
  as_vs_mock_subset <- as_vs_mock %>% filter(Name2 %in% toupper(hamster_to_human$external_gene_name))
  as_vs_mock_subset <- as_vs_mock_subset %>% full_join(hamster_to_human, by=c("Name2"="external_gene_name"))
  saveRDS(as_vs_mock_subset, file = paste0("dir/vs_mock_DEGs/",
                                           ct, "_as_vs_mockDEG.rds"))
  
  ## run KEGG
  
  kegg_input <- as_vs_mock_subset %>%
    dplyr::select(c(entrezgene_id, avg_log2FC)) %>%
    na.omit () %>% # remove unmapped
    arrange(desc(avg_log2FC)) %>%
    deframe()
  
  
  # gene set enrichment from kegg lists
  as_vs_mock_kegg <- gseKEGG(geneList = kegg_input,
                             organism     = 'hsa',
                             minGSSize    = 20,
                             maxGSSize = 1000,
                             pvalueCutoff = 1,
                             verbose      = FALSE)
  
  # tidy result
  as_vs_mock_keggTidy <- as_vs_mock_kegg@result 
  as_vs_mock_keggTidy <- as_vs_mock_keggTidy %>%
    as_tibble() %>%
    #filter(p.adjust <= 0.05) %>%
    arrange(desc(NES)) 
  #as_vs_mock_keggTidy <- as_vs_mock_keggTidy[,1:(length(as_vs_mock_keggTidy)-3)]
  as_vs_mock_keggTidy$vax <- "AlumS"
  ## SAVE!
  saveRDS(as_vs_mock_keggTidy, file = paste0("dir/vs_mock_KEGG/",
                                             ct, "_as_vs_mockKEGG.rds"))
  
  ## run GO
  
  
  as_vs_mock_GO <- gseGO(geneList = kegg_input,
                         #keyType = "ENSEMBL",
                         ont = "BP",
                         minGSSize    = 20,
                         maxGSSize = 1000,
                         pvalueCutoff = 1,
                         OrgDb = org.Hs.eg.db, 
                         verbose      = FALSE)
  
  # tidy result
  as_vs_mock_GOTidy <- as_vs_mock_GO@result 
  as_vs_mock_GOTidy <- as_vs_mock_GOTidy %>%
    as_tibble() %>%
    arrange(desc(NES)) 
  #as_vs_mock_GOTidy <- as_vs_mock_GOTidy[,1:(length(as_vs_mock_GOTidy)-3)]
  as_vs_mock_GOTidy$vax <- "AlumS"
  saveRDS(as_vs_mock_GOTidy, file = paste0("dir/vs_mock_GO/",
                                           ct, "_as_vs_mockGO.rds"))
  
  
  ####################### MEVS ####################################
  ident1 <- paste0(ct, "_MeVS")
  ident2 <- paste0(ct, "_Mock")
  
  
  # get DEGs
  mevs_vs_mock <- FindMarkers(hamster, ident.1 = ident1, ident.2 = ident2, logfc.threshold = 0 )
  mevs_vs_mock$gene <- rownames(mevs_vs_mock)
  # convert to human equivalents
  mevs_vs_mock["Name2"] <- toupper(mevs_vs_mock$gene)
  mevs_vs_mock_subset <- mevs_vs_mock %>% filter(Name2 %in% toupper(hamster_to_human$external_gene_name))
  mevs_vs_mock_subset <- mevs_vs_mock_subset %>% full_join(hamster_to_human, by=c("Name2"="external_gene_name"))
  saveRDS(mevs_vs_mock_subset, file = paste0("dir/vs_mock_DEGs/",
                                             ct, "_mevs_vs_mockDEG.rds"))
  
  ## run KEGG
  
  kegg_input <- mevs_vs_mock_subset %>%
    dplyr::select(c(entrezgene_id, avg_log2FC)) %>%
    na.omit () %>% # remove unmapped
    arrange(desc(avg_log2FC)) %>%
    deframe()
  
  
  # gene set enrichment from kegg lists
  mevs_vs_mock_kegg <- gseKEGG(geneList = kegg_input,
                               organism     = 'hsa',
                               minGSSize    = 20,
                               maxGSSize = 1000,
                               pvalueCutoff = 1,
                               verbose      = FALSE)
  
  # tidy result
  mevs_vs_mock_keggTidy <- mevs_vs_mock_kegg@result 
  mevs_vs_mock_keggTidy <- mevs_vs_mock_keggTidy %>%
    as_tibble() %>%
    #filter(p.adjust <= 0.05) %>%
    arrange(desc(NES)) 
  #mevs_vs_mock_keggTidy <- mevs_vs_mock_keggTidy[,1:(length(mevs_vs_mock_keggTidy)-3)]
  mevs_vs_mock_keggTidy$vax <- "MeVS"
  ## SAVE!
  saveRDS(mevs_vs_mock_keggTidy, file = paste0("dir/vs_mock_KEGG/",
                                               ct, "_mevs_vs_mockKEGG.rds"))
  mevs_vs_mock_keggTidy
  ## run GO
  
  mevs_vs_mock_GO <- gseGO(geneList = kegg_input,
                           #keyType = "ENSEMBL",
                           ont = "BP",
                           minGSSize    = 20,
                           maxGSSize = 1000,
                           pvalueCutoff = 1,
                           OrgDb = org.Hs.eg.db, 
                           verbose      = FALSE)
  
  # tidy result
  mevs_vs_mock_GOTidy <- mevs_vs_mock_GO@result 
  mevs_vs_mock_GOTidy <- mevs_vs_mock_GOTidy %>%
    as_tibble() %>%
    arrange(desc(NES)) 
  #mevs_vs_mock_GOTidy <- mevs_vs_mock_GOTidy[,1:(length(mevs_vs_mock_GOTidy)-3)]
  mevs_vs_mock_GOTidy$vax <- "MeVS"
  # SAVE!
  saveRDS(mevs_vs_mock_GOTidy, file = paste0("dir/vs_mock_GO/",
                                             ct, "_mevs_vs_mockGO.rds"))
  
  
  ### combine results
  kegg_merged_res <- bind_rows(as_vs_mock_keggTidy, mevs_vs_mock_keggTidy)
  saveRDS(kegg_merged_res, file=paste0("dir/vs_mock_KEGG/",
                                       ct, "_mergedres_KEGG.rds"))
  
  GO_merged_res <- bind_rows(as_vs_mock_GOTidy, mevs_vs_mock_GOTidy)
  saveRDS(GO_merged_res, file=paste0("dir/vs_mock_GO/",
                                     ct, "_mergedres_GO.rds"))
  
}



#### CREATE BUBBLEPLOTS TO VISUALISE RESULTS ####

ct_list <- c("AlveolarMacrophages", "InterstitialMacrophages", "Treml4+Macrophages", "TNKcells")

#create counter
p <- 0
q <- 0


for (ct in ct_list){
  
  
  ## KEGG
  
  #set counter
  p <- p + 1
  
  merged_res <- readRDS(paste0("dir/vs_mock_KEGG/",
                               ct, "_mergedres_KEGG.rds"))
  
  
  ## calculate gene ratio using core enrichment genes
  gene_count<- merged_res %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
  merged_res <- left_join(merged_res, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
  
  # get NES scores for MeVS and AlumS
  mevs_min <- merged_res %>%
    filter(vax == "MeVS") %>%
    dplyr::select(c(ID, Description, NES)) %>%
    rename(mevs_NES = NES)
  
  
  as_min <- merged_res %>%
    filter(vax == "AlumS") %>%
    dplyr::select(c(ID, Description, NES)) %>%
    rename(as_NES = NES) 
  # identify top_pways
  top_pways <- plyr::join(mevs_min, as_min)
  top_pways[is.na(top_pways)] <- 0  
  
  top_pways_up <- top_pways %>%
    mutate(diff = abs(mevs_NES - as_NES)) %>%
    arrange(desc(diff)) %>%
    slice_max(order_by = diff, n = 30) %>%
    pull(ID)
  
  # subset by top pways and scale NES using scale() fxn to obtain zscores
  merged_res_plot <- merged_res %>%
    filter(ID %in% top_pways_up) %>%
    mutate(zscore = scale(NES)) %>% # 
    mutate(signifi = ifelse(qvalues < 0.05, 1, 0)) %>%
    mutate(stroke = ifelse(signifi == 1, 1.5, 0))
  
  plot <- merged_res_plot %>%
    ggplot(aes(y=Description, x=vax, size=GeneRatio, fill=NES, color=as.factor(signifi), stroke=stroke))+
    geom_point(shape=21)+
    #custom_theme+
    theme_classic()+
    scale_x_discrete(labels=c("Alum + S", expression(MeV[vac2]-SARS2-S(H))))+
    theme(axis.text.x = element_text( angle = 45, size = 8, hjust=1 )) +
    scale_size(range=c(2,7), name="size ratio")+
    scale_fill_gradient2(mid="white", high="#0808C2", low="#FF6306", name="z-score")+
    scale_color_manual(values=c("black", "transparent"), name="FDR", labels=c("<0.05", "ns"))+
    labs(x="", y="")+
    guides(gradient2=guide_legend(order=1), color=guide_legend(order=2), size=guide_legend(order=3)) +
    ggtitle(ct) +
    theme(plot.title = element_text(hjust = -0, size = 8)) 
  
  
  
  assign(paste0("p", p), plot)
  
  
  
  ## GO
  
  #set counter
  q <- q + 1
  
  merged_res <- readRDS(paste0("dir/vs_mock_GO/",
                               ct, "_mergedres_GO.rds"))
  
  ## calculate gene ratio using core enrichment genes
  gene_count<- merged_res %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
  merged_res <- left_join(merged_res, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
  
  
  # get NES scores for MeVS and AlumS
  mevs_min <- merged_res %>%
    filter(vax == "MeVS") %>%
    dplyr::select(c(ID, Description, NES)) %>%
    rename(mevs_NES = NES)
  
  
  as_min <- merged_res %>%
    filter(vax == "AlumS") %>%
    dplyr::select(c(ID, Description, NES)) %>%
    rename(as_NES = NES) 
  # identify top_pways
  top_pways <- plyr::join(mevs_min, as_min)
  top_pways[is.na(top_pways)] <- 0  
  
  top_pways_up <- top_pways %>%
    mutate(diff = abs(mevs_NES - as_NES)) %>%
    arrange(desc(diff)) %>%
    slice_max(order_by = diff, n = 30) %>%
    pull(ID)
  
  # subset by top pways
  merged_res_plot <- merged_res %>%
    filter(ID %in% top_pways_up)%>%
    mutate(zscore = scale(NES)) %>% # 
    mutate(signifi = ifelse(qvalues < 0.05, 1, 0)) %>%
    mutate(stroke = ifelse(signifi == 1, 1.5, 0))
  
  plot1 <- merged_res_plot %>%
    ggplot(aes(y=Description, x=vax, size=GeneRatio, fill=NES, color=as.factor(signifi), stroke=stroke))+
    geom_point(shape=21)+
    #custom_theme+
    theme_classic()+
    scale_x_discrete(labels=c("Alum + S", expression(MeV[vac2]-SARS2-S(H))))+
    theme(axis.text.x = element_text( angle = 45, size = 8, hjust=1 )) +
    scale_size(range=c(2,7), name="size ratio")+
    scale_fill_gradient2(mid="white", high="#0808C2", low="#FF6306", name="z-score")+
    scale_color_manual(values=c("black", "transparent"), name="FDR", labels=c("<0.05", "ns"))+
    labs(x="", y="")+
    guides(gradient2=guide_legend(order=1), color=guide_legend(order=2), size=guide_legend(order=3)) +
    ggtitle(ct) +
    theme(plot.title = element_text(hjust = -0, size = 8)) 
  
  assign(paste0("q", q), plot1)
  
  
  
  
}

#### arrange plots in grid and export as pdf and svg

kegg_grid <- ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="right", align = "v" )
ggsave("dir/scPathway_kegg_bubbleplots.pdf", 
       kegg_grid, width=12, height=11.5)
ggsave("dir/scPathway_kegg_bubbleplots.eps", 
       kegg_grid, width=12, height=11.5)
#######
go_grid <- ggarrange(q1, q2, q3, q4, ncol=2, nrow=2, common.legend = TRUE, legend="right", align = "v" )
ggsave("dir/scPathway_go_bubbleplots.pdf", 
       go_grid, width=14, height=11.5)
ggsave("dir/scPathway_go_bubbleplots.eps", 
       go_grid, width=14, height=11.5)

