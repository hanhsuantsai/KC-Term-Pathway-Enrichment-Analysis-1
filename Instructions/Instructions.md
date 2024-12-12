# KC Term Pathway Enrichment Analysis

This repository provides the necessary resources and scripts for performing pathway enrichment analysis using Key Characteristics (KC) terms. The repository supports both Over-Representation Analysis (ORA) and Gene Set Enrichment Analysis (GSEA), allowing users to input their gene lists and identify enriched KC terms. The Key Characteristics (KC) term pathway enrichment analysis leverages two primary databases: **KEGG** and **Reactome**, to associate a list of genes with specific pathways that may be overrepresented.  

In this repository, we provide two functions using the clusterProfiler package in R:  
**ORA** (Over-Representation Analysis): Identifies pathways enriched based on significantly differentially expressed genes (DEGs).  
**GSEA** (Gene Set Enrichment Analysis): Detects pathways based on ranked gene lists.

Users can input their own gene lists (either as a vector of DEGs for ORA or a ranked list for GSEA), and these functions will return the enriched KC terms based on KEGG and Reactome gene sets.


## A. Folder Structure

The repository is organized into three main folders for better management of files and scripts:

- **data/**: Contains raw and processed data files used for pathway enrichment.
  - `KC_gene_set_KEGG.csv`: Contains KC terms and genes associated with the KC terms in the KEGG database.
  - `KC_gene_set_REACTOME.csv`: Contains KC terms and genes associated with the KC terms in the Reactome database.
  - `KC_Terms.csv`: Contains descriptions of the KC terms.
  
- **scripts/**: Contains the R scripts used for the pathway enrichment analysis.
  - `ORA_function.R`: Script for running Over-Representation Analysis (ORA).
  - `GSEA_function.R`: Script for running Gene Set Enrichment Analysis (GSEA).

- **results/**: Enrichment results will be automatically saved in this directory after running the scripts.

## B. Automated Data Loading
The datasets required for the analysis are automatically loaded through the provided R code. You do not need to manually download or specify paths for the following datasets:  
- `KC_gene_set_KEGG.csv:` Contains KC terms and genes associated with the KC terms in the KEGG database.
- `KC_gene_set_REACTOME.csv:` Contains KC terms and genes associated with the KC terms in the Reactome database.
- `KC_Terms.csv:` Contains descriptions of the KC terms.


## C. Usage

### 1. Supported Gene Symbols for KC Term Pathway Enrichment Analysis
The analysis script only supports **human gene symbols**. Make sure your list of genes contains correct human gene symbols before starting the analysis.



### 2. Dependencies
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
### 3. Loading the data 
Run the following R code to load the datasets:
```r
# Loading KC gene sets Reactome from GitHub repository
KC_gene_set_REACTOME <- read.csv(url("https://raw.githubusercontent.com/kingdave-hub/KC-Term-Pathway-Enrichment-Analysis/main/data/KC_gene_set_REACTOME.csv"))

# Loading KC gene sets KEGG from GitHub repository
KC_gene_set_KEGG <- read.csv(url("https://raw.githubusercontent.com/kingdave-hub/KC-Term-Pathway-Enrichment-Analysis/main/data/KC_gene_set_KEGG.csv"))
```

### 4.Running the Code
**ORA (Over-Representation Analysis)**
This function performs an over-representation analysis using a list of differentially expressed genes (DEGs).

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

### How to run ORA:

1. **Load your DEG list:** Prepare your DEG list as a vector.
2. **Load your background gene list:** Prepare your background gene list as another vector.
3. **Load the annotation data:** The annotation data contains KC terms associated with genes from KEGG or Reactome. This is already loaded in your R.
4. **Run the ORA function:** Call the clusterProfiler_ORA function with your inputs.

### GSEA (Gene Set Enrichment Analysis)
This function runs GSEA using a ranked gene list and identifies the KC terms that are significantly enriched.

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
### How to run GSEA:
1. Prepare a ranked gene list as input.
2. **Load the annotation data:** The annotation data contains KC terms associated with genes from KEGG or Reactome. This is already loaded in your R.
3. **Run the Gsea function:** Call the clusterProfiler_GSEA function with your inputs.

## 5. Save Results
Results from the enrichment analysis will be saved into the `results` folder. This folder is created automatically if it doesn't already exist.
Run the following code to save the result:
```r
# Ensure the "results" folder exists
if (!dir.exists("results")) {
  dir.create("results")
}

# Saving ORA and GSEA results into the 'results' folder
write.csv(ora_results, file = "results/ora_results.csv", row.names = FALSE)
write.csv(gsea_results, file = "results/gsea_results.csv", row.names = FALSE)
```
Files that will be saved include:
- `ora_results.csv`
- `gsea_results.csv`

## 6. Run Example Dataset
A script named run_example_analysis.R is provided in the scripts/ folder to demonstrate how to run ORA and GSEA on example datasets in the `data/` folder.  
Before running this script, ensure you have followed the above instructions to install the necessary R packages, and load the required functions.

To execute the script, use the following command in R:

```r
## run example dataset
source(url("https://raw.githubusercontent.com/kingdave-hub/KC-Term-Pathway-Enrichment-Analysis/main/scripts/run_example_analysis.R"))
```
This will automatically load example dataset from the data/ folder, perform the analyses, and save the results in the results/ folder.

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
