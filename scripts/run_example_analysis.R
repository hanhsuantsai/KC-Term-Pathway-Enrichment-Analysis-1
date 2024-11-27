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


sort <- DEG %>%
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
