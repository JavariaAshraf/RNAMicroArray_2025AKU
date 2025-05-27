# Investigating host transcriptional responses through RNA Microarray Analysis - Learn to Use R
Workshop on transcriptome sequencing analysis 
# 🧬 Transcriptomics Workshop

Welcome to the official repository for the **Transcriptomics Data Analysis Workshop**. This hands-on workshop is designed to guide participants through the essential steps of transcriptome data analysis, from raw reads to biological interpretation using modern bioinformatics tools and R programming.

---

## 📅 Workshop Details

- **Title:Investigating host transcriptional responses through RNA Microarray Analysis - Learn to Use R ** 
- **Level:** Beginner to Intermediate  
-- **Audience:** Life scientists, bioinformaticians, students, and researchers working with gene expression data  
- **Prerequisites:** Basic command-line and R knowledge preferred

---

## 🧾 Learning Objectives

By the end of this workshop, participants will be able to:

- Understand the experimental design and workflow of Microarray Transcriptome
- Conduct differential gene expression analysis
- Perform functional enrichment (GO, KEGG)
- Visualize transcriptomic data using heatmaps, volcano plots, and barplots
- Interpret biological meaning from statistical outputs

---

## 🧰 Topics Covered

1. **Preparing Gene Lists for Enrichment**
     - Converting gene symbols to Entrez IDs or other identifiers
     - Filtering up/down-regulated gene sets
     - Creating ranked gene lists (optional for GSEA)
2. **Functional Enrichment Analysis Using clusterProfiler**
   - Gene Ontology (GO) Enrichment
      - enrichGO() for BP, MF, CC
   - KEGG Pathway Analysis
      - enrichKEGG() for pathway-level insights
   - WikiPathways Enrichment
      - enrichWP() via the clusterProfiler + ReactomePA or wikipathways packages
   -Gene Set Enrichment Analysis (GSEA)
      -  gseGO(), gseKEGG() 

3. **Visualization of Results**
   - Dot plots dotplot()
   - Bar plots (barplot())
   - Enrichment maps (emapplot())
   - GO tree plots (goplot())
   - KEGG pathway diagrams (pathview package)

---
## STEP#1. Loading Packages and setting Working Directory
Load required packages in **R**:

```r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("DOSE")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("rWikiPathways")
library(rWikiPathways)
library(DOSE)
library(pathview)
library(clusterProfiler)
install.packages('devtools')
require(devtools)
library(org.Hs.eg.db)
library(enrichplot)
library(rWikiPathways)
setwd("C:/Path to Folder Workshop")
 ```

## STEP2. Gene list Preparation
To prepare your gene list, you have to convert gene symbols to Entrez IDs. Processing of genelist in **R**:

```r
# Load your csv file into dataframe
d <- read.csv("Analysis_CoVvsHealthy.csv", header = TRUE)

# Extract data from genelist
geneList <- d[, 2]
#d[, 2] means: select all rows from the second column
#Selects the first column of the same data frame "d[,1]", **as.character(...)** ensures the values are treated as character strings
names(geneList) <- as.character(d[, 1])

# Sort the gene list (e.g., fold changes)
geneList <- sort(geneList, decreasing = TRUE)

# View the top genes
head(geneList)
# View the genelist
View(d)

# Filter significant genes (e.g., |logFC| > 2)
gene <- names(geneList)[abs(geneList) > 2]
head(gene)
 ```
## STEP3. Over Representation of Analysis
GO, KEGG and wikiPathway Analysis **R**:

```r
# GO Ontology (Biological Processes)
enrich_GO <- enrichGO(gene = gene,
                      OrgDb = org.Hs.eg.db,
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05)
head(enrich_GO)
#Convert to Entrez IDs
enrich_GO <- setReadable(enrich_GO, 'org.Hs.eg.db', 'ENTREZID')
#Saving results in a cluster summary and writing it to csv file
cluster_summary <- data.frame(enrich_GO)
write.csv(cluster_summary, "Enrich-analysis-GO.csv")
# GO Ontology (Molecular Function)
enrich_MF <- enrichGO(gene = gene,
                      OrgDb = org.Hs.eg.db,
                      readable = T,
                      ont = "MF",
                      pvalueCutoff = 0.05)

head(enrich_MF)
#Dotplot
dotplot(enrich_GO)
#Enriched GO induced graph:
goplot(enrich_GO, showCategory = 10)

# KEGG Pathway Analysis
Ekegg <- enrichKEGG(gene=gene,organism='hsa',
                   pvalueCutoff = 0.05)
head(Ekegg)
#Convert to Entrez IDs
Ekegg <- setReadable(Ekegg, 'org.Hs.eg.db', 'ENTREZID')
#Saving results in a cluster summary and writing it to csv file
cluster_summary <- data.frame(Ekegg)
write.csv(cluster_summary, "Enrich-analysis-KEGG.csv")
#Barplot
barplot(Ekegg, 
        showCategory = 10, 
        title = "Enriched Pathways",
        font.size = 8)
#heatmap
heatplot(Ekegg, foldChange=geneList, showCategory=10)+ ggtitle("Enrich heatmap of KEGG Pathway")

# wiki Pathway Analysis
wikienrich <- enrichWP(names(geneList), organism = "Homo sapiens")
head(wikienrich)
wikienrich <- setReadable(wikienrich, 'org.Hs.eg.db', 'ENTREZID')
cluster_summary <- data.frame(wikienrich)
write.csv(cluster_summary, "wiki-Enrich-workshop.csv")
#Barplot
barplot(wikienrich, 
        showCategory = 10, 
        title = "Enriched Pathways Of Wiki Pathways",
        font.size = 8)
#Heatmap
heatplot(wikienrich, foldChange=geneList, showCategory=10)+ ggtitle("Enrich heatmap of wiki Pathway")
 ```
## STEP3. Gene set Enrichment Analysis
GO, KEGG and wikiPathway Analysis **R**:

```r

 ```
## 🗂️ Repository Structure

