
## This assumes you have run the DESeq2_file_preparation script
# and that you have a counts dataframe and a coldata dataframe

library(DESeq2)
library(tidyverse)
library(cowplot)
library(biomaRt)

## Expression analysis of

counts <- as.data.frame(read.csv("./data/RAB234_counts.csv", row.names = 1))
coldata = read.csv("./data/coldata.csv")
counts = counts[,-1]

### MRPL28 KO cells ------------------------------------------------------------
# Going for the simplest design and comparing Mrpl28 KO Trpn vs Trpp

counts = counts %>%
  dplyr::select(starts_with("Mrpl"))

coldata = coldata %>% filter(celltype=="Mrpl")

# Create the full model for comparison of samples

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~condition) 

# Generate a linear model

dds$condition <- relevel(dds$condition, "Trpp")
dds <- DESeq(dds)

resultsNames(dds)

# Checking size factors and dispersion

sizeFactors(dds) # only takes into account the sequencing depth, looks ok
plotDispEsts(dds) # verifies normalization, graph looks a-ok


## Trpp vs Trpn for mrpl28ko cells treatment -----------------------------------

res <- results(dds, name = "condition_Trpn_vs_Trpp")

res_tbl <- as_tibble(res, rownames="ensembl")

# Import annotation file

ensembl_to_geneName <- read_csv("./data/Biomart_annotations_mm10.csv") %>% 
  filter(!duplicated(gene, drop = F))

# Add gene names and entrez id to your results

res_tbl <- res_tbl %>%
  left_join(ensembl_to_geneName, by = c("ensembl"="ENSMUG")) %>%
  arrange(padj) 

write.csv(res_tbl,"./data_output/Mrpl_Trpn_vs_Trpp/res_tbl.csv", row.names = T)

hist(res_tbl$padj) # distribution looks ok, if not filtering is needed

plotMA(res) # looks agood


# Save the signif genes

library("openxlsx")
library(ggrepel)

overexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange>=1) 
underexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange<=-1)
signif <- full_join(overexpr, underexpr)

list_of_datasets <- list("UP" = overexpr, "DOWN" = underexpr,"ALL_sign" = signif)
write.xlsx(list_of_datasets, file = "./data_output/Mrpl_Trpn_vs_Trpp/RAB234_sign_genes.xlsx")

# Volcano plot

res_tbl %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1, 
             label = ifelse(padj<0.05&log2FoldChange>=1|
                              padj<0.05&log2FoldChange<=-1,as.character(gene),''))) +
  scale_colour_manual(values = c("gray", "firebrick3")) +
  labs(col="Significantly expressed")+
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)+
  geom_text_repel(max.overlaps = 20,
                  box.padding = 0.25,
                  segment.color = 'grey50',
                  fontface = "italic",
                  size = 3)+
  labs(title = "Trp- vs Trp+ in Not nucleofected cells")+
  theme_bw()

ggsave("./figures/Mrpl_Trpn_vs_Trpp/volcanoplot.png", last_plot(), device = png, dpi= 500, width = 8, height = 6)

saveRDS(dds, file = "./data_output/Mrpl_Trpn_vs_Trpp/dds.rds")
saveRDS(res_tbl, file = "./data_output/Mrpl_Trpn_vs_Trpp/res_tbl.rds")
write_csv(res_tbl, file = "./data_output/Mrpl_Trpn_vs_Trpp/res_tbl.csv")
res_tbl1 <-  res_tbl %>% filter(!is.na(ENTREZID), padj<0.05)
saveRDS(res_tbl1, file = "./data_output/Mrpl_Trpn_vs_Trpp/res_tbl_signif.rds")
write_csv(res_tbl1, file = "./data_output/Mrpl_Trpn_vs_Trpp/res_tbl_signif.csv")


### NDUFA2 KO cells ------------------------------------------------------------
# Going for the simplest design and comparing Ndufa2 KO Trpn vs Trpp

counts <- as.data.frame(read.csv("./data/RAB234_counts.csv", row.names = 1))
coldata = read.csv("./data/coldata.csv")
counts = counts[,-1]

counts = counts %>%
  dplyr::select(starts_with("Nduf"))

coldata = coldata %>% filter(celltype=="Nduf")

# Create the full model for comparison of samples

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~condition) 

# Generate a linear model

dds$condition <- relevel(dds$condition, "Trpp")
dds <- DESeq(dds)

resultsNames(dds)

# Checking size factors and dispersion

sizeFactors(dds) # only takes into account the sequencing depth, looks ok
plotDispEsts(dds) # verifies normalization, graph looks a-ok


## Trpp vs Trpn for ndufa2ko cells treatment -----------------------------------

res <- results(dds, name = "condition_Trpn_vs_Trpp")

res_tbl <- as_tibble(res, rownames="ensembl")

# Import annotation file

ensembl_to_geneName <- read_csv("./data/Biomart_annotations_mm10.csv") %>% 
  filter(!duplicated(gene, drop = F))

# Add gene names and entrez id to your results

res_tbl <- res_tbl %>%
  left_join(ensembl_to_geneName, by = c("ensembl"="ENSMUG")) %>%
  arrange(padj) 

write.csv(res_tbl,"./data_output/Nduf_Trpn_vs_Trpp/res_tbl.csv", row.names = T)

hist(res_tbl$padj) # distribution looks ok, if not filtering is needed

plotMA(res) # looks agood


# Save the signif genes

library("openxlsx")
library(ggrepel)

overexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange>=1) 
underexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange<=-1)
signif <- full_join(overexpr, underexpr)

list_of_datasets <- list("UP" = overexpr, "DOWN" = underexpr,"ALL_sign" = signif)
write.xlsx(list_of_datasets, file = "./data_output/Nduf_Trpn_vs_Trpp/RAB234_sign_genes.xlsx")

# Volcano plot

res_tbl %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1, 
             label = ifelse(padj<0.05&log2FoldChange>=1|
                              padj<0.05&log2FoldChange<=-1,as.character(gene),''))) +
  scale_colour_manual(values = c("gray", "firebrick3")) +
  labs(col="Significantly expressed")+
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)+
  geom_text_repel(max.overlaps = 20,
                  box.padding = 0.25,
                  segment.color = 'grey50',
                  fontface = "italic",
                  size = 3)+
  labs(title = "Trp- vs Trp+ in Not nucleofected cells")+
  theme_bw()

ggsave("./figures/Nduf_Trpn_vs_Trpp/volcanoplot.png", last_plot(), device = png, dpi= 500, width = 8, height = 6)

saveRDS(dds, file = "./data_output/Nduf_Trpn_vs_Trpp/dds.rds")
saveRDS(res_tbl, file = "./data_output/Nduf_Trpn_vs_Trpp/res_tbl.rds")
write_csv(res_tbl, file = "./data_output/Nduf_Trpn_vs_Trpp/res_tbl.csv")
res_tbl1 <-  res_tbl %>% filter(!is.na(ENTREZID), padj<0.05)
saveRDS(res_tbl1, file = "./data_output/Nduf_Trpn_vs_Trpp/res_tbl_signif.rds")
write_csv(res_tbl1, file = "./data_output/Nduf_Trpn_vs_Trpp/res_tbl_signif.csv")

