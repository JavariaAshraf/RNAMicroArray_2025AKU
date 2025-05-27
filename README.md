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
## Preparing Gene Lists for Enrichment

To prepare your gene list, you may want to convert gene symbols to Entrez IDs or other identifiers. Here's an example in **R**:

```r
# Read in your differential expression results
d <- read.csv("Analysis_CoVvsHealthy.csv", header = TRUE)

# Extract data from genelist
geneList <- d[, 2]
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

## 🗂️ Repository Structure

