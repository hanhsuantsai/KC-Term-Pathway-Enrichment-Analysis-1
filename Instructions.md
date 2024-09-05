# 1. Introduction to KC Term Pathway Enrichment Analysis
Pathway enrichment analysis is a method used to identify significantly enriched biological pathways or functional categories from a given gene set. The Key Characteristics (KC) term pathway enrichment analysis leverages two primary databases: KEGG and Reactome, to associate a list of genes with specific pathways that may be overrepresented.

In this repository, we provide two functions using the clusterProfiler package in R:

**ORA** (Over-Representation Analysis): Identifies pathways enriched based on significantly differentially expressed genes (DEGs).
**GSEA** (Gene Set Enrichment Analysis): Detects pathways based on ranked gene lists.

Users can input their own gene lists (either as a vector of DEGs for ORA or a ranked list for GSEA), and these functions will return the enriched KC terms based on KEGG and Reactome gene sets.

# 2. Dependencies
To run the code successfully, you need to install and load the following R packages. Ensure your R setup has these installed:
```r
# Check if the required packages are installed and install them if not
required_packages <- c("clusterProfiler", "tidyverse")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load the libraries
library(clusterProfiler)
library(tidyverse)
```
# 3.Data Files
This repository contains the following essential data files:
- `KC_KEGG_genes.Rdata:` Contains KC terms and genes associated with the KC terms in the KEGG database.
- `KC_Reactome_genes.Rdata:` Contains KC terms and genes associated with the KC terms in the Reactome database.
- `KC_Terms.csv:` Contains descriptions of the KC terms.


**Loading the Data:** You need to load the .Rdata files in your R session before running the analysis.  
Here's how to do that:
- locate the path to your downloaded dataset
- specify the path where appropriate in the codes below before running.
```r
### load KC gene sets_Reactome
load(file = "path to KC_Reactome_genes.Rdata") # replace with path to your KC_Reactome_genes.Rdata

b <- list()
for (i in 1:length(KC_Reactome_genes)){
  a <- matrix(NA, length(KC_Reactome_genes[[i]]), 2)
  colnames(a) <- c("KC_term", "gene")
  a[,1] <- names(KC_Reactome_genes)[i]
  a[,2] <- KC_Reactome_genes[[i]]
  b[[i]] <- a  
}
KC_gene_set_reactome <- do.call("rbind", b) %>% as.data.frame()

### KEGG
load(file = "path to KC_KEGG_genes.Rdata") # replace with path to your KC_KEGG_genes.Rdata

b <- list()
for (i in 1:length(KC_KEGG_genes)){
  a <- matrix(NA, length(KC_KEGG_genes[[i]]), 2)
  colnames(a) <- c("KC_term", "gene")
  a[,1] <- names(KC_KEGG_genes)[i]
  a[,2] <- KC_KEGG_genes[[i]]
  b[[i]] <- a 
}
KC_gene_set_KEGG <- do.call("rbind", b) %>% as.data.frame()
```

# 4.Running the Code
**ORA (Over-Representation Analysis)**
This function performs an over-representation analysis using a list of differentially expressed genes (DEGs).

**Input Requirements:**

`data:` A vector of DEGs.  
`annotation:` A data frame with two columns: first, KC terms, and second, genes associated with those terms.  
`background:` A vector of background genes.  
`treatment:` A label for the experimental condition (character).  
`sig_level:` The significance threshold (usually 0.05).  
```r
##  clusterProfiler
clusterProfiler_ORA <- function(data, # a list of DEG
                                annotation, # two columns, first with term information and second with gene associates with the term
                                background, # list of background genes
                                treatment, # experimental condition
                                sig_level) {# significant level (fdr) for enrichment analysis
  
  KC_Terms <- read.csv(file = "KC_Terms.csv")  # replace with your path to KC_Terms.csv
  
  eTerm <- enricher(gene = data,
                    universe = background,
                    TERM2GENE=annotation)
  
  output <- eTerm@result %>%
    filter(p.adjust < sig_level) %>%
    mutate(treatment = treatment) %>%
    subset(, select = c("ID",
                        "p.adjust",
                        "geneID",
                        "Count"))
  
  colnames(output)[1] <- "KC_Term"
  
  output <- left_join(output, KC_Terms, by = "KC_Term")
  
}
```

### How to Run ORA:

1. Load your own DEG list as a vector and define your background genes.
2. Run the `clusterProfiler_ORA` function with your inputs.


### GSEA (Gene Set Enrichment Analysis)
This function runs GSEA using a ranked gene list and identifies the KC terms that are significantly enriched.

**Input Requirements:**
- `data:` A ranked gene list (vector).
- `annotation:` A data frame with two columns: KC terms and associated genes.
- `treatment:` A label for the experimental condition (character).
- `sig_level:` The significance threshold (usually 0.05).
```r
clusterProfiler_GSEA <- function(data,         # ranked gene list
                                 annotation,   # TERM2GENE annotation
                                 treatment,    # experimental condition label
                                 sig_level) {  # significance level (FDR cutoff)
  
  KC_Terms <- read.csv(file = "KC_Terms.csv")  # replace with your path to KC_Terms.csv
  
  eTerm <- GSEA(data,                          # ranked gene list
                TERM2GENE = annotation,        # gene set annotations
                pvalueCutoff = sig_level,      # cutoff for significance
                pAdjustMethod = "fdr",         # adjust p-values using FDR
                by = "fgsea")                  # method for running GSEA
  
  output <- eTerm@result %>%
    filter(p.adjust < sig_level) %>%           # filter by significance level
    mutate(treatment = treatment) %>%          # add treatment label
    subset(, select = c("ID",                  # select relevant columns
                        "p.adjust",
                        "core_enrichment",
                        "NES",
                        "treatment"))
  
  colnames(output)[1] <- "KC_Term"             # rename ID to KC_Term
  
  output <- left_join(output, KC_Terms, by = "KC_Term")  # add KC term descriptions
  
  return(output)                               # return the processed output
}
```
### How to Run GSEA:
1. Prepare a ranked gene list as input.
2. Run the clusterProfiler_GSEA function with your data and parameters.

## 4. Expected Outcome
Both ORA and GSEA will return a data frame with the following columns:
- **KC_Term:** The enriched KC term.
- **p.adjust:** Adjusted p-value (FDR).
- **geneID or core_enrichment:** The genes driving the enrichment.
- **Count or NES:** Number of genes enriched (ORA) or normalized enrichment score (GSEA).
- **treatment:** The experimental condition label.
-  **KC_Term Description:** Description of the KC Term enriched.

