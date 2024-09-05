##  clusterProfiler
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