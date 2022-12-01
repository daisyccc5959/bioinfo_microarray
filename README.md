# Bioinformatic Profiling of Candidate Gene Involved in Endometrial Cancer Development and Progression.

## Summary
Endometrial Cancer (EC) is a common type of cancer developed in the uterus. EC is the most common type of cancer in uterine cancers for its high incidence and mortality rate. There is a lack of precise biomarkers for early diagnosis. This study aims to identify the potential crucial genes to provide insight into a biomarker in EC cancer progression by bioinformatics analysis.  
GSE17025 and GSE39099 were downloaded from Gene Expression Omnibus (GEO) for EC gene profiling. The gene expression pattern of GSE17025 was stage I EC. 584 differentially expressed genes (DEGs) were screened with the control comparison. Following with two strategies to identify functional genes, (1) Protein-Protein Interaction (PPI) network construction with hub gene clustering and (2) gene co-expression network construction from case-control associated modules (WGCNA). The overlapping genes between those two strategies were selected as candidate genes. Further, GSE39099 was used to validate specific target genes expressed in the early stages. The gene expression pattern was from different stages of EC (1) atypical endometrial hyperplasia, (2) stages I and II EC, (3) stages III and IV EC, and the control. 440 DEGs were acquired and performed to the same method as a strategy for hub gene clusters. With the validation hub gene clusters, the target gene of the overlapping was *Matrix Metalloproteinase-2* (*MMP2*).  
MMP2 was calculated as a down-regulated gene in both GSE17025 and GSE39099 datasets. Kaplan-Meier survival rate in MMP2 demonstrated that low MMP2 was associated with shorter survival (p= 0.005) compared with high MMP2 expression in uterine corpus endometrial carcinoma patients (n=543). Those findings provided that MMP2 might be the potential target in EC progression. Furthermore, the pipeline is suitable for functional gene identification discovery.  

## Flow chart
![flow_chart](https://user-images.githubusercontent.com/117828593/205051190-187c9b7a-0614-47d7-ba7e-f3f636be89c0.png)

## Data source
  * [GSE17025](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17025)
  * [GSE39099](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39099)
  <img width="800" alt="download_txt" src="https://user-images.githubusercontent.com/117828593/205051735-03f5d0dd-9689-4518-b617-8414a6d86801.png">

## Reqiurement
* [R (version 4.2.2)](https://www.r-project.org/) and [Rstudio (IDE)](https://posit.co/downloads/)
* packages:
    1. data preprocessing:
        * `readr`
        * `ggplot2`
        * `reshape2`
        * `biocManager: GEOquery` and `hgu133plus2.db`  
    ```R
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("GEOquery")
    BiocManager::install("hgu133plus2.db")
    ```
    2. DEGs (differentially expressed genes):
        * `limma` (for GSE17025, w/ biological replicates)
        * `edgeR` (for GSE39099, w/o biological replicates)
        * `gplots`
    3. PPI network construction:
        * [STRING](https://string-db.org)  
          download 
          <img width="1073" alt="string_load_deg" src="https://user-images.githubusercontent.com/117828593/205051337-3ef68a8c-0914-45a6-b5cb-a77e517ad249.png">
          <img width="1139" alt="string_dowload_predictedoutput" src="https://user-images.githubusercontent.com/117828593/205051415-11d66300-72f3-4679-acf8-d7ff2f3a1453.png">
       
        * [Cytoscape](https://cytoscape.org)  
          app(MCODE)
    4. Gene co-expression network construction: 
       * `WGCNA`  
   ```R
   install.packages("BiocManager") 
   BiocManager::install("WGCNA") 
   ```
    5. [KMplot](https://kmplot.com/analysis/)
