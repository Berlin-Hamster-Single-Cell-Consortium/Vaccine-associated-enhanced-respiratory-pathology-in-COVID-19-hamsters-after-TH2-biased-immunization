#### RNAseq Bulk Expression Analysis #####
## Ebening, Mühlenbach et al. 2022 ##

#### library ####
library(tidyverse)
library(GOplot)
library(ggpubr)

#### data path in ####
eb_df <- read_csv("./eb_df.csv")
GO_df <- read_csv("./GO_df.csv")

#### data out ####
out <- "./Plots"
setwd(out)


#### theme and color ####
custom_theme <-
  list(
    theme_classic() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        text = element_text(size = 12),
        legend.position = "right",
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
      )
  )


#### generate GO plot tables #### 
## generate GO walter Gene/GO lists ## 

GOwalter <- data.frame(Category=GO_df$id, ID=paste("GO:", GO_df$`GO term`, sep=""),
                       Term=GO_df$Description, Genes=GO_df$`Detected Genes (Names)`, adj_pval=GO_df$`FDR p-value`)

GeneWalter <- data.frame(ID=eb_df$Name, logFC=eb_df$log_Fold_change, AveExpr=eb_df$`Max group means`,
                         P.Value=eb_df$`P-value`, adj.P.Val=eb_df$FDR, Category=eb_df$id)


## separate circle data generation per category/id

# "Alum+SvsOptiMEM" matches "Alum+S vs. OptiMEM eb.xlsx"
circ_AlumvsOpti <- circle_dat(filter(GOwalter, Category=="Alum+SvsOptiMEM"), filter(GeneWalter, Category=="Alum+S vs. OptiMEM eb.xlsx"))

# "MeVvac2vsOptiMEM"  
circ_MeVvac2vsOptiMEM <- circle_dat(filter(GOwalter, Category=="MeVvac2vsOptiMEM"), filter(GeneWalter, Category=="MeVvac2-SARS2-S(H) vs. OptiMEM eb.xlsx"))

# "mock_vs_infected_optiMEM" 
circ_MockvsInfOpti <- circle_dat(filter(GOwalter, Category=="mock_vs_infected_optiMEM"), filter(GeneWalter, Category=="mock_vs_infected_optiMEM eb.xlsx"))

# "mockvsAlum+S" 
circ_mockvsAlumS <- circle_dat(filter(GOwalter, Category=="mockvsAlum+S"), filter(GeneWalter, Category=="mock vs. Alum+S eb.xlsx"))

# "mockvsMeVvac2"  
circ_mockvsMeVvac2 <- circle_dat(filter(GOwalter, Category=="mockvsMeVvac2"), filter(GeneWalter, Category=="mock vs. MeVvac2-SARS2-S(H) eb.xlsx"))


# combine all GO-walter circle data into a single df
circle_df <- rbind(circ_AlumvsOpti, circ_MeVvac2vsOptiMEM, circ_mockvsAlumS, circ_MockvsInfOpti, circ_mockvsMeVvac2)
# 
# # filter for selected terms - pval and logFC
# genes_sel <- c(
#   "GO:0060337",
#   "GO:0009615",
#   "GO:0006950",
#   "GO:0034341",
#   "GO:0001666",
#   "GO:0034097",
#   "GO:0043408",
#   "GO:0012501",
#   "GO:0045597",
#   "GO:0009892",
#   "GO:0002443",
#   "GO:0006954",
#   "GO:0006955",
#   "GO:0006959",
#   "GO:0006952",
#   "GO:1990869",
#   "GO:0048245",
#   "GO:0035711",
#   "GO:0035706",
#   "GO:0035712",
#   "GO:0035707",
#   "GO:0035712"
# )


## binary pval ##
# adj_pval
circle_df["sig"] <- circle_df$adj_pval
circle_df <- circle_df %>% mutate(sig=replace(sig, adj_pval <= 0.05, "s"))
circle_df <- circle_df %>% mutate(sig=replace(sig, adj_pval > 0.05, "ns"))


## summarize by term ##
# pval for all terms the same
# logFC as mean
# one zscore per term

circ_AlumvsOpti_sum <- circ_AlumvsOpti %>%
  group_by(ID, term) %>%
  summarise(GO_zscore=mean(zscore), pval=mean(adj_pval), logFC=mean(logFC))
circ_AlumvsOpti_sum["id"] <- rep("circ_AlumvsOpti_sum", nrow(circ_AlumvsOpti_sum))

circ_MeVvac2vsOptiMEM_sum <- circ_MeVvac2vsOptiMEM %>%
  group_by(ID, term) %>%
  summarise(GO_zscore=mean(zscore), pval=mean(adj_pval), logFC=mean(logFC))
circ_MeVvac2vsOptiMEM_sum["id"] <- rep("circ_MeVvac2vsOptiMEM_sum", nrow(circ_MeVvac2vsOptiMEM_sum))

circ_MockvsInfOpti_sum <- circ_MockvsInfOpti %>%
  group_by(ID, term) %>%
  summarise(GO_zscore=mean(zscore), pval=mean(adj_pval), logFC=mean(logFC))
circ_MockvsInfOpti_sum["id"] <- rep("circ_MockvsInfOpti_sum", nrow(circ_MockvsInfOpti_sum))

circ_mockvsAlumS_sum <- circ_mockvsAlumS %>%
  group_by(ID, term) %>%
  summarise(GO_zscore=mean(zscore), pval=mean(adj_pval), logFC=mean(logFC))
circ_mockvsAlumS_sum["id"] <- rep("circ_mockvsAlumS_sum", nrow(circ_mockvsAlumS_sum))

circ_mockvsMeVvac2_sum <- circ_mockvsMeVvac2 %>%
  group_by(ID, term) %>%
  summarise(GO_zscore=mean(zscore), pval=mean(adj_pval), logFC=mean(logFC))
circ_mockvsMeVvac2_sum["id"] <- rep("circ_mockvsMeVvac2_sum", nrow(circ_mockvsMeVvac2_sum))

# combine
zscore_df <- rbind(circ_AlumvsOpti_sum, circ_MeVvac2vsOptiMEM_sum, circ_MockvsInfOpti_sum, circ_mockvsAlumS_sum, circ_mockvsMeVvac2_sum)


##### combine GO_df with circPlot #####
# circ_AlumvsOpti_sum
GO_df_AlumSvsOptiMEM <- GO_df %>% filter(id=="Alum+SvsOptiMEM")
GO_df_AlumSvsOptiMEM$`GO term` <- paste("GO:", GO_df_AlumSvsOptiMEM$`GO term`, sep = "")
AlumvsOpti_join_dot <- circ_AlumvsOpti_sum %>% full_join(GO_df_AlumSvsOptiMEM, by=c("ID"="GO term"))

# circ_MeVvac2vsOptiMEM_sum
GO_df_MeVvac2vsOptiMEM <- GO_df %>% filter(id=="MeVvac2vsOptiMEM")
GO_df_MeVvac2vsOptiMEM$`GO term` <- paste("GO:", GO_df_MeVvac2vsOptiMEM$`GO term`, sep = "")
MeVvac2vsOptiMem_join_dot <- circ_MeVvac2vsOptiMEM_sum %>% full_join(GO_df_MeVvac2vsOptiMEM, by=c("ID"="GO term"))

# circ_MockvsInfOpti_sum
GO_df_mock_vs_infected_optiMEM <- GO_df %>% filter(id=="mock_vs_infected_optiMEM")
GO_df_mock_vs_infected_optiMEM$`GO term` <- paste("GO:", GO_df_mock_vs_infected_optiMEM$`GO term`, sep = "")
MockvsInfectedoptiMEM_join_dot <- circ_MockvsInfOpti_sum %>% full_join(GO_df_mock_vs_infected_optiMEM, by=c("ID"="GO term"))

# circ_mockvsAlumS_sum
GO_df_mockvsAlum <- GO_df %>% filter(id=="mockvsAlum+S")
GO_df_mockvsAlum$`GO term` <- paste("GO:", GO_df_mockvsAlum$`GO term`, sep = "")
MockvsAlumS_join_dot <- circ_mockvsAlumS_sum %>% full_join(GO_df_mockvsAlum, by=c("ID"="GO term"))

# circ_mockvsMeVvac2_sum
GO_df_MockvsMeVvac2 <- GO_df %>% filter(id=="mockvsMeVvac2")
GO_df_MockvsMeVvac2$`GO term` <- paste("GO:", GO_df_MockvsMeVvac2$`GO term`, sep = "")
MockvsMeVvac2_join_dot <- circ_mockvsMeVvac2_sum %>% full_join(GO_df_MockvsMeVvac2, by=c("ID"="GO term"))

# bind
GO_df_dot <- rbind(AlumvsOpti_join_dot, MeVvac2vsOptiMem_join_dot, MockvsInfectedoptiMEM_join_dot, MockvsAlumS_join_dot, MockvsMeVvac2_join_dot)
GO_df_dot["size_ratio"] <- (GO_df_dot$`DE Genes`) / (GO_df_dot$`Detected Genes`)


#### dotplots ####
GO_df_dot["signifi"] <- GO_df_dot$`FDR p-value`
GO_df_dot <- GO_df_dot %>% mutate(signifi=replace(signifi, signifi > .05, 1))
GO_df_dot <- GO_df_dot %>% mutate(signifi=replace(signifi, signifi <= .05, 0))

## only to mock comparison ##
GO_df_dot <- GO_df_dot %>% filter(id.y!="Alum+SvsOptiMEM" & id.y!="MeVvac2vsOptiMEM")

genes_sel <- c("GO:0006950", # response to stress
               "GO:0006952",  # defense response
               "GO:0006955", # immune response
               "GO:0009615", # response to virus
               "GO:0006954", # inflammatory response
               "GO:0006959", # humoral immune response

               "GO:0034097", #  response to cytokine
               "GO:0060337", # type 1 inf signaling pathway
               "GO:0034341", #  response to interferon-gamma
               "GO:0043408", # regulation of MAPK cascade
               "GO:0002443", # leukocyte mediated immunity
               "GO:1990869", # cellular response to chemokine
               "GO:0048245", # eosinophil chemotaxis
               "GO:0035711", # T-helper 1 cell activation
               "GO:0035712", # T-helper 2 cell activation

               "GO:0045597", #  positive regulation of cell differentiation
               "GO:0009892", # negative regulation of metabolic process

               "GO:0001666", #  response to hypoxia
               "GO:0012501" # programmed cell death
               )

genes_sel_names <- c("programmed cell death",
                     "response to hypoxia",
                     "negative regulation of metabolic process",
                     "positive regulation of cell differentiation",
                     "T-helper 2 cell activation",
                     "T-helper 1 cell activation",
                     "eosinophil chemotaxis",
                     "cellular response to chemokine",
                     "leukocyte mediated immunity",
                     "regulation of MAPK cascade",
                     "response to interferon-gamma",
                     "type I interferon signaling pathway",
                     "response to cytokine",
                     "humoral immune response",
                     "inflammatory response",
                     "response to virus",
                     "immune response",
                     "defense response",
                     "response to stress"
                     )

GO_df_dot2 <- GO_df_dot %>% filter(ID %in% genes_sel)
GO_df_dot2 <- GO_df_dot2 %>% mutate(ID=factor(ID, levels=genes_sel))
GO_df_dot2 <- GO_df_dot2 %>% mutate(term=factor(term, levels=genes_sel_names))

# add delta stroke size
GO_df_dot2["stroke"] <- GO_df_dot2$signifi
GO_df_dot2[which(GO_df_dot2$stroke==1),]$stroke <- .75
GO_df_dot2[which(GO_df_dot2$stroke==0),]$stroke <- 1.5

##### dotplot 1 #####
dotplotr <- GO_df_dot2 %>%
  ggplot(aes(y=term, x=id.y, size=`size_ratio`, fill=`GO_zscore`, color=as.factor(signifi), stroke=stroke))+
  geom_point(shape=21)+
  custom_theme+
  scale_x_discrete(labels=c("mock_vs_infected_optiMEM"= "SARS-CoV-2", 
                            "mockvsAlum+S"= "Alum+S",
                            "mockvsMeVvac2" = "MeVac2S"))+
  scale_y_discrete(labels=c("programmed cell death", "response to hypoxia",
                            "negative regulation of metabolic process", "positive regulation of cell differentiation",
                            "T-helper 2 cell activation"= "TH 2 cell activation", "T-helper 1 cell activation"="TH 1 cell activation",
                            "eosinophil chemotaxis", "cellular response to chemokine",
                            "leukocyte mediated immunity", "regulation of MAPK cascade",
                            "response to interferon-gamma" = "response to IFN-gamma", "type I interferon signaling pathway"="type I IFN signaling pathway",
                            "response to cytokine", "humoral immune response",
                            "inflammatory response", "response to virus",
                            "immune response", "defense response", "response to stress"))+
  scale_size(range=c(2,7), name="size ratio")+
  scale_fill_gradient2(mid="white", high="#0808C2", low="#FF6306", name="z-score")+
  scale_color_manual(values=c("black", "grey20"), name="FDR", labels=c("<0.5", "ns"))+
  labs(x="", y="")+
  guides(gradient2=guide_legend(order=1), color=guide_legend(order=2), size=guide_legend(order=3))
#  ggtitle("Dotplot selected genes by cell paper")
dotplotr


##### dotplot 2 #####
dotplotl_terms <- c("response to interleukin-4",
                    "T cell mediated immunity",
                    "cytokine production involved in immune response",
                    "negative regulation of interleukin-1 beta secretion",
                    "regulation of CD40 signaling pathway",
                    "diaphragm contraction",
                    "involuntary skeletal muscle contraction",
                    "positive regulation of thymocyte migration",
                    "secretion",
                    "circulatory system process",
                    "positive regulation of tumor necrosis factor secretion",
                    "regulation of T cell tolerance induction",
                    "phagocytosis, engulfment",
                    "establishment or maintenance of apical/basal cell polarity",
                    "collagen metabolic process",
                    "NADPH oxidation",
                    "histone H3-K27 trimethylation",
                    "positive regulation of T-helper 2 cell cytokine production",
                    "regulation of viral entry into host cell")

dotplotl_df <- GO_df_dot %>%
  filter(term %in% dotplotl_terms)

dotplotl_df$signifi

# add delta stroke size
dotplotl_df["stroke"] <- dotplotl_df$signifi
dotplotl_df[which(dotplotl_df$stroke==1),]$stroke <- .75
dotplotl_df[which(dotplotl_df$stroke==0),]$stroke <- 1.5

# save data for dotplot
write.csv(dotplotl_df, "dotplotl_df_dotplot.csv")

dotplotl <- dotplotl_df %>%
  ggplot(aes(y=term, x=id.y, size=`size_ratio`, fill=`GO_zscore`, color=as.factor(signifi), stroke=stroke))+
  geom_point(shape=21)+
  custom_theme+
  scale_x_discrete(labels=c("mock_vs_infected_optiMEM"= "SARS-CoV-2", 
                            "mockvsAlum+S"= "Alum+S",
                            "mockvsMeVvac2" = "MeVac2S"))+
  scale_size(range=c(2,7), name="size ratio")+
  scale_y_discrete(labels=c("response to interleukin-4" = "response to IL-4",
                            "positive regulation of tumor necrosis factor secretion" = "positive regulation of TNF secretion",
                            "negative regulation of interleukin-1 beta secretion" = "negative regulation of IL-1 beta secretion",
                            "positive regulation of T-helper 2 cell cytokine production" = "positive regulation of TH 2 cell cytokine production"))+
  scale_fill_gradient2(mid="white", high="#0808C2", low="#FF6306", name="z-score")+
  scale_color_manual(values=c("black", "grey20"), name="FDR", labels=c("<0.5", "ns"))+
  labs(x="", y="")+
  guides(gradient2=guide_legend(order=1), color=guide_legend(order=2), size=guide_legend(order=3))
dotplotl

##### combine #####
ggarrange(dotplotr, dotplotl, nrow=1, common.legend = FALSE,  align = c("hv"))
ggsave("Dotplot_dlegend5.pdf", width = 13, height = 8)
ggsave("Dotplot_dlegend5.eps", width = 13, height = 8, device = "eps")



#### Chordplots ####
moregens <- read_table("./GenesForHeatmap.txt", col_names=FALSE)

moreTerms <- c("immune response",
               "innate immune response",
               "adaptive immune response",
               "regulation of leukocyte activation",
               "regulation of leukocyte mediated immunity",
               "regulation of macrophage activation",
               "response to cytokine",
               "leukocyte chemotaxis",
               "inflammatory response",
               "regulation of apoptotic process")

GeneWalter2 <- GeneWalter
GeneWalter2$ID <- toupper(GeneWalter$ID)

##### single chords #####
## 1
chord_MockvsInfOpti <- chord_dat(data = circ_MockvsInfOpti, genes = filter(GeneWalter2, Category=="mock_vs_infected_optiMEM eb.xlsx" & ID %in% toupper(moregens$X1)),
                                 process = moreTerms)

chord_MockvsInfOpti2 <- chord_MockvsInfOpti
chord_MockvsInfOpti2 <- cbind(chord_MockvsInfOpti2,rownum=seq(1, nrow(chord_MockvsInfOpti)))

chord_MockvsInfOpti_arranged <- chord_MockvsInfOpti2 %>% as.data.frame() %>% arrange(wt=logFC) # arrange order

chord_MockvsInfOpti <- chord_MockvsInfOpti[chord_MockvsInfOpti_arranged$rownum,]

GOChord_MockvsInfOpti <- GOChord(chord_MockvsInfOpti, space = 0.02, gene.order = "none", gene.space = 0.25, gene.size = 5,lfc.min=-10, lfc.max=10,
                                 nlfc=1)
GOChord_MockvsInfOpti


## 2
chord_MockvsAlum <- chord_dat(data = circ_mockvsAlumS, genes = filter(GeneWalter2, Category=="mock vs. Alum+S eb.xlsx" & ID %in% toupper(moregens$X1)),
                              process = moreTerms)
chord_MockvsAlum <- chord_MockvsAlum[chord_MockvsInfOpti_arranged$rownum,]

GOChord_MockvsAlum <- GOChord(chord_MockvsAlum, space = 0.02, gene.order = 'none', gene.space = 0.25, gene.size = 5,nlfc=1,lfc.min=-10, lfc.max=10)
GOChord_MockvsAlum


## 3 
chord_MockvsMeVac <- chord_dat(data = circ_mockvsMeVvac2, genes = filter(GeneWalter2, Category=="mock vs. MeVvac2-SARS2-S(H) eb.xlsx" & ID %in% toupper(moregens$X1)),
                               process = moreTerms)
# due to 0 expression IL19 is excluded - add it - meVac only
chord_MockvsMeVac2 <- rbind(chord_MockvsMeVac, "IL19"=chord_MockvsAlum["IL19",])
chord_MockvsMeVac2["IL19", "logFC"] <- 0
chord_MockvsMeVac <- chord_MockvsMeVac2[rownames(chord_MockvsInfOpti_arranged),]

GOChord_MockvsMeVac <- GOChord(chord_MockvsMeVac, space = 0.02, gene.order = 'none', gene.space = 0.25, gene.size = 5,nlfc=1,lfc.min=-10, lfc.max=10)
GOChord_MockvsMeVac


##### combine #####
ggarrange(GOChord_MockvsInfOpti+
            scale_fill_gradient2(low= "#FF6306", mid="white", high = "#0808C2"),
          GOChord_MockvsAlum+
            scale_fill_gradient2(low= "#FF6306", mid="white", high = "#0808C2"),
          GOChord_MockvsMeVac+
            scale_fill_gradient2(low= "#FF6306", mid="white", high = "#0808C2"), 
          nrow=1, common.legend = TRUE)

ggsave("GO_line_arranged.pdf",width=30, height = 11.75)
