# Vaccine-associated-enhanced-respiratory-pathology-in-COVID-19-hamsters-after-TH2-biased-immunization
This is the GitHub repository providing access to all the code used in the analysis of the bulk RNA sequencing and single cell RNA sequencing data in the manuscript "Vaccine-associated enhanced respiratory pathology in COVID-19 hamsters after TH2-biased immunization" (Ebenig _et al_., 2022). Code utilised for analysis of scRNAseq data is in the folder "scRNAseq", while the code used to analyse the bulk RNA sequencing data is in the folder "BulkRNAseq". 
For scRNAseq analysis, the code has been separated into four segments and the order in which the data has been processed is indicated by a prepended number, i.e. "1_", "2_" and so on. For bulk RNAseq, all of the analysis and the code for creation of the visualisations is contained within one R script.

## Abstract

Vaccine-associated enhanced respiratory disease (VAERD) has been a severe complication for some respiratory infections. To investigate the potential for VAERD induction in COVID-19, we evaluated two vaccine leads utilizing a severe hamster infection model: a TH1-biased measles vaccine-derived candidate and a TH2-biased alum-adjuvanted, non-stabilized Spike protein. The MeV-derived vaccine protected the animals, but the protein lead induced VAERD, which could be alleviated by dexamethasone treatment. Bulk transcriptomic analysis revealed that our protein vaccine prepared enhanced host gene dysregulation in the lung, exclusively up-regulating mRNAs encoding the eosinophil attractant CCL-11, TH2-driving IL-19, or TH2-cytokines IL-4, IL-5, and IL-13. scRNA-Seq identified lung macrophages or lymphoid cells as sources, respectively. Our findings imply VAERD is caused by the concerted action of hyperstimulated macrophages and TH2 cytokine-secreting lymphoid cells and potentially links VAERD to antibody-dependent enhancement (ADE). In summary, we identify the cytokine drivers and cellular contributors, which mediate VAERD after TH2-biased vaccination. 

![Ebenig et al  Graphical abstract revised R2](https://user-images.githubusercontent.com/61689250/179394082-74b96c8c-517a-4c79-8158-b8b2b4ed1040.jpg)

## Data availability

Raw data and the processed Seurat object for the single cell RNA sequencing analysis has been uploaded to GEO with the accession number [GSE196938](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196938). The raw data for bulk RNA sequencing analysis has been uploaded to GEO with accession number [GSE195939](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195939), while the processed count data can be downloaded from [FigShare](https://figshare.com/projects/Vaccine-associated-enhanced-respiratory-pathology-in-COVID-19-hamsters-after-TH2-biased-immunization/142469).

