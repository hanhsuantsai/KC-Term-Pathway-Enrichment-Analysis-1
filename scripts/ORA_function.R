
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