
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
# Going for the simplest design and comparing Scr Trpn vs Trpp

counts = counts %>%
  dplyr::select(starts_with("Scr"))

coldata = coldata %>% filter(celltype=="Scr")

# Create the full model for comparison of samples

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~condition) 

# Generate a linear model

dds$condition <- relevel(dds$condition, "Trpn")
dds <- DESeq(dds)

resultsNames(dds)

# Checking PCA

rld <- rlogTransformation(dds)

p1 <- plotPCA(rld,intgroup="condition")
p2 <- plotPCA(rld,intgroup="celltype") 
p3 <- plotPCA(rld,intgroup=c("condition","celltype")) 

p3+
  coord_fixed(ratio = 1.5)

print(p1)  

plot_grid(p1, p2, p3, nrow = 1, align = "v")


# Checking sample similarity

sampleDists <- dist(t(assay(dds)))

library("RColorBrewer")
library(pheatmap)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(dds$condition, dds$celltype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Checking size factors and dispersion

sizeFactors(dds) # only takes into account the sequencing depth, looks ok
plotDispEsts(dds) # verifies normalization, graph looks a-ok


## Trpp vs Trpn for scramble cells treatment -----------------------------------

res <- results(dds, name = "condition_Trpp_vs_Trpn")

res_tbl <- as_tibble(res, rownames="ensembl")

# Import annotation file

ensembl_to_geneName <- read_csv("./data/Biomart_annotations_mm10.csv") %>% 
  filter(!duplicated(gene, drop = F))

# Add gene names and entrez id to your results

res_tbl <- res_tbl %>%
  left_join(ensembl_to_geneName, by = c("ensembl"="ENSMUG")) %>%
  arrange(padj) 

write.csv(res_tbl,"./data_output/Scr_Trpp_vs_Trpn/res_tbl.csv", row.names = T)

hist(res_tbl$padj) # distribution looks ok

plotMA(res) # looks agood

# Save the signif genes

library("openxlsx")
library(ggrepel)

overexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange>=1) 
underexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange<=-1)
signif <- full_join(overexpr, underexpr)

list_of_datasets <- list("UP" = overexpr, "DOWN" = underexpr,"ALL_sign" = signif)
write.xlsx(list_of_datasets, file = "./data_output/Scr_Trpp_vs_Trpn/RAB234_sign_genes.xlsx")

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
  geom_text_repel(max.overlaps = 15,
                  box.padding = 0.25,
                  segment.color = 'grey50',
                  fontface = "italic")+
  labs(title = "Trp+ vs Trp- in Scramble cells")+
  theme_bw()

ggsave("./figures/Scr_Trpp_vs_Trpn/volcanoplot.png", last_plot(), device = png, dpi= 500, width = 8, height = 6)


# save the tibbles etc

saveRDS(dds, file = "./data_output/Scr_Trpp_vs_Trpn/dds.rds")
saveRDS(res_tbl, file = "./data_output/Scr_Trpp_vs_Trpn/res_tbl.rds")
write_csv(res_tbl, file = "./data_output/Scr_Trpp_vs_Trpn/res_tbl.csv")
res_tbl1 <-  res_tbl %>% filter(!is.na(ENTREZID), padj<0.05)
saveRDS(res_tbl1, file = "./data_output/Scr_Trpp_vs_Trpn/res_tbl_signif.rds")
write_csv(res_tbl1, file = "./data_output/Scr_Trpp_vs_Trpn/res_tbl_signif.csv")

