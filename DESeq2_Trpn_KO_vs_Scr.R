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

### SCRAMBLE cells -------------------------------------------------------------
# Going for the simplest design and comparing Nduf or Mrpl to Scr in Trpp

counts = counts %>%
  dplyr::select(starts_with(c("Mrpl_Trpn","Scr_Trpn")))

coldata = coldata %>% filter(condition=="Trpn", celltype %in%c("Mrpl","Scr"))

# Create the full model for comparison of samples

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~celltype) 

# Generate a linear model

dds$celltype <- relevel(dds$celltype, "Scr")
dds <- DESeq(dds)

resultsNames(dds)

## Mrpl vs Scr in Trpp----------------------------------------------------------

res <- results(dds, name = "celltype_Mrpl_vs_Scr")
res_tbl <- as_tibble(res, rownames="ensembl")

# Import annotation file

ensembl_to_geneName <- read_csv("./data/Biomart_annotations_mm10.csv") %>% 
  filter(!duplicated(gene, drop = F))

# Add gene names and entrez id to your results

res_tbl <- res_tbl %>%
  left_join(ensembl_to_geneName, by = c("ensembl"="ENSMUG")) %>%
  arrange(padj) 

write.csv(res_tbl,"./data_output/Trpn_Mrpl_vs_Scr/res_tbl.csv", row.names = T)

hist(res_tbl$padj) # distribution looks ok

plotMA(res) # looks agood

# Save the signif genes

library("openxlsx")
library(ggrepel)

overexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange>=1) 
underexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange<=-1)
signif <- full_join(overexpr, underexpr)

list_of_datasets <- list("UP" = overexpr, "DOWN" = underexpr,"ALL_sign" = signif)
write.xlsx(list_of_datasets, file = "./data_output/Trpn_Mrpl_vs_Scr/RAB234_sign_genes.xlsx")

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
  labs(title = "Mrpl28 KO vs scramble cells in Trpn medium")+
  theme_bw()

ggsave("./figures/Trpn_Mrpl_vs_Scr/volcanoplot.png", last_plot(), device = png, dpi= 500, width = 8, height = 6)

# save the tibbles etc

saveRDS(dds, file = "./data_output/Trpn_Mrpl_vs_Scr/dds.rds")
saveRDS(res_tbl, file = "./data_output/Trpn_Mrpl_vs_Scr/res_tbl.rds")
write_csv(res_tbl, file = "./data_output/Trpn_Mrpl_vs_Scr/res_tbl.csv")
res_tbl1 <-  res_tbl %>% filter(!is.na(ENTREZID), padj<0.05)
saveRDS(res_tbl1, file = "./data_output/Trpn_Mrpl_vs_Scr/res_tbl_signif.rds")
write_csv(res_tbl1, file = "./data_output/Trpn_Mrpl_vs_Scr/res_tbl_signif.csv")



## Nduf vs Scr in Trpp----------------------------------------------------------

counts <- as.data.frame(read.csv("./data/RAB234_counts.csv", row.names = 1))
coldata = read.csv("./data/coldata.csv")
counts = counts[,-1]

counts = counts %>%
  dplyr::select(starts_with(c("Nduf_Trpn","Scr_Trpn")))

coldata = coldata %>% filter(condition=="Trpn", celltype %in%c("Nduf","Scr"))

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~celltype) 
# Generate a linear model

dds$celltype <- relevel(dds$celltype, "Scr")
dds <- DESeq(dds)

res <- results(dds, name = "celltype_Nduf_vs_Scr")
res_tbl <- as_tibble(res, rownames="ensembl")

# Import annotation file

ensembl_to_geneName <- read_csv("./data/Biomart_annotations_mm10.csv") %>% 
  filter(!duplicated(gene, drop = F))

# Add gene names and entrez id to your results

res_tbl <- res_tbl %>%
  left_join(ensembl_to_geneName, by = c("ensembl"="ENSMUG")) %>%
  arrange(padj) 

write.csv(res_tbl,"./data_output/Trpn_Nduf_vs_Scr/res_tbl.csv", row.names = T)

hist(res_tbl$padj) # distribution looks ok

plotMA(res) # Almost no difference


# Save the signif genes

library("openxlsx")
library(ggrepel)

overexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange>=1) 
underexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange<=-1)
signif <- full_join(overexpr, underexpr)

list_of_datasets <- list("UP" = overexpr, "DOWN" = underexpr,"ALL_sign" = signif)
write.xlsx(list_of_datasets, file = "./data_output/Trpn_Nduf_vs_Scr/RAB234_sign_genes.xlsx")

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
  labs(title = "Ndufa2 KO vs scramble cells in Trpn medium")+
  theme_bw()

ggsave("./figures/Trpn_Nduf_vs_Scr/volcanoplot.png", last_plot(), device = png, dpi= 500, width = 8, height = 6)

# save the tibbles etc

saveRDS(dds, file = "./data_output/Trpn_Nduf_vs_Scr/dds.rds")
saveRDS(res_tbl, file = "./data_output/Trpn_Nduf_vs_Scr/res_tbl.rds")
write_csv(res_tbl, file = "./data_output/Trpn_Nduf_vs_Scr/res_tbl.csv")
res_tbl1 <-  res_tbl %>% filter(!is.na(ENTREZID), padj<0.05)
saveRDS(res_tbl1, file = "./data_output/Trpn_Nduf_vs_Scr/res_tbl_signif.rds")
write_csv(res_tbl1, file = "./data_output/Trpn_Nduf_vs_Scr/res_tbl_signif.csv")
