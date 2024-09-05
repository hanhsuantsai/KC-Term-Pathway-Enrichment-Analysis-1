# KC Term Pathway Enrichment Analysis

This repository provides the necessary resources and scripts for performing pathway enrichment analysis using Key Characteristics (KC) terms with KEGG and Reactome gene sets. The repository supports both Over-Representation Analysis (ORA) and Gene Set Enrichment Analysis (GSEA), allowing users to input their gene lists and identify enriched KC terms.

## 1. Folder Structure

The repository is organized into three main folders for better management of files and scripts:

- **data/**: Contains raw and processed data files used for pathway enrichment.
  - `KC_KEGG_genes.Rdata`: Contains KC terms and genes associated with the KC terms in the KEGG database.
  - `KC_Reactome_genes.Rdata`: Contains KC terms and genes associated with the KC terms in the Reactome database.
  - `KC_Terms.csv`: Contains descriptions of the KC terms.
  
- **scripts/**: Contains the R scripts used for the pathway enrichment analysis.
  - `ORA_function.R`: Script for running Over-Representation Analysis (ORA).
  - `GSEA_function.R`: Script for running Gene Set Enrichment Analysis (GSEA).
  
- **docs/**: Documentation files, including the usage guide for running the scripts.
  - `Instructions.md`: Step-by-step guide on how to run the pathway enrichment analysis.


## Usage

To run the enrichment analysis, follow the steps outlined in the **Instructions.md** file located in the `docs/` folder. The file provides guidance on:

- Dependencies and installation
- Running ORA and GSEA with your gene lists

## Contact
For questions, contact [owarekings@tamu.edu](mailto:owarekings@tamu.edu).
