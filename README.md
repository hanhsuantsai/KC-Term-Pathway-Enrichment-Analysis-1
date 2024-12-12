# KC Term Pathway Enrichment Analysis  

This repository enables pathway enrichment analysis using Key Characteristics (KC) terms with **Over-Representation Analysis (ORA)** and **Gene Set Enrichment Analysis (GSEA)**. Users can input their gene lists and identify enriched pathways leveraging KEGG and Reactome databases. Follow the instructions below to successfully perform ORA and GSEA in R.  


**Note:** The analysis script only supports **human gene symbols**. Ensure your input gene list contains valid human gene symbols.


### 1. Load Dependencies
To run the code successfully, you need to install and load the following R packages.
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
### 2. Load Gene Sets
Run the following R code to load Gene Sets:
```r
# Loading KC gene sets Reactome from GitHub repository
KC_gene_set_REACTOME <- read.csv(url("https://raw.githubusercontent.com/kingdave-hub/KC-Term-Pathway-Enrichment-Analysis/main/data/KC_gene_set_REACTOME.csv"))
head(KC_gene_set_REACTOME)

# Loading KC gene sets KEGG from GitHub repository
KC_gene_set_KEGG <- read.csv(url("https://raw.githubusercontent.com/kingdave-hub/KC-Term-Pathway-Enrichment-Analysis/main/data/KC_gene_set_KEGG.csv"))
head(KC_gene_set_KEGG)
```

### 3.Running the Code
**ORA (Over-Representation Analysis)**
This section defines the function to perform over-representation analysis using a list of differentially expressed genes (DEGs).

```r
# clusterProfiler
clusterProfiler_ORA <- function(data, # a list of DEG
                                annotation, # select either KC_gene_set_KEGG or KC_gene_set_REACTOME
                                background, # list of background genes
                                treatment, # experimental condition
                                sig_level) {# significant level (fdr) for enrichment analysis
  
KC_Terms <- read.csv(url("https://raw.githubusercontent.com/kingdave-hub/KC-Term-Pathway-Enrichment-Analysis/main/data/KC_Terms.csv"))
  
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


### GSEA (Gene Set Enrichment Analysis)
This section defines the function to perform GSEA using a list of differentially expressed genes (DEGs).

```r
clusterProfiler_GSEA <- function(data,         # ranked gene list
                                 annotation,   # select either KC_gene_set_KEGG or KC_gene_set_REACTOME
                                 treatment,    # experimental condition label
                                 sig_level) {  # significance level (FDR cutoff)
  
KC_Terms <- read.csv(url("https://raw.githubusercontent.com/kingdave-hub/KC-Term-Pathway-Enrichment-Analysis/main/data/KC_Terms.csv"))
  
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

## 4. Run the analysis with an example dataset
This step demonstrates how to input data and use the ORA/GSEA functions with an example dataset (available in the `data` folder) after completing steps 1-3.  
Use the provided example to understand the workflow, including data preparation, sorting significant genes, and calling the functions.

```r
## read the csv from url in data

data <- read.csv(url("https://raw.githubusercontent.com/kingdave-hub/KC-Term-Pathway-Enrichment-Analysis/main/data/demo.csv"))

res <-data  %>%
  filter(padj < 0.1) %>%
  filter(abs(log2FoldChange) > 1)
 
annotation <- KC_gene_set_REACTOME

background <- data$gene

treatment <- "test"

sig_level <- 1


ORA<- clusterProfiler_ORA(data =res$gene,
                          annotation = KC_gene_set_REACTOME,
                          background = background,
                          treatment = treatment,
                          sig_level = sig_level) %>%
  subset(, select = c(1,2,5,6))%>%
  mutate(Enrichment = "ORA") %>%
  mutate(Treatment = treatment)

# Save ORA results
write.csv(ORA, "results/ORA_results.csv", row.names = FALSE)


sort <- data %>%
  filter(!is.na(gene)) %>%
  filter(log2FoldChange != 0) %>%
  group_by(gene) %>%
  summarise(log2FoldChange = max(log2FoldChange)) %>% as.data.frame() %>%
  arrange(desc(log2FoldChange))

ranked_genes <- setNames(sort$log2FoldChange, sort$gene)
GSEA <- clusterProfiler_GSEA(ranked_genes, annotation,
                             treatment,
                             sig_level) %>%
  subset(, select = c(1,2,6,7)) %>%
  mutate(Enrichment = "GSEA") %>%
  mutate(Treatment = treatment)

# Save GSEA results
write.csv(GSEA, "results/GSEA_results.csv", row.names = FALSE)
```


## Citations

This project uses the following R packages:

- **tidyverse**:
   -Wickham, H., & RStudio. (2019). _tidyverse: Easily install and load the 'tidyverse'_. R package version 1.2.1. Available at: [https://CRAN.R-project.org/package=tidyverse](https://CRAN.R-project.org/package=tidyverse)

- **clusterProfiler**:
  -Yu, G., Wang, L.-G., Han, Y., & He, Q.-Y. (2012). _clusterProfiler: an R package for comparing biological themes among gene clusters_. OMICS: A Journal of Integrative Biology, 16(5), 284-287. Available at: [https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html](https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html)

## License
This project is licensed under the MIT License - see the [LICENSE.txt](LICENSE.txt) file for details.


## Contact
For questions, contact [owarekings@tamu.edu](mailto:owarekings@tamu.edu).
