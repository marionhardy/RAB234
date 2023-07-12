# RAB234
## Project summary

Mouse-isolated T cells treated with Trp+ or - medium for 72 hours.

- Non nucleofected cells : Nonuc
- Scramble cells : Scr
- Mrpl28 KO :  Mrpl
- Ndufa2 KO: Nduf

Biological triplicates, RNAseq.

## Analysis pipeline

Comparaisons done using DESeq2 (v1.36.0)

In tryptophan high medium:
- Nonuc vs Scr
- Mrpl vs Scr
- Nduf vs Scr

In tryptophan low medium:
- Nonuc vs Scr
- Mrpl vs Scr
- Nduf vs Scr

Other comparisons:
- Nonuc Trpn vs Trpp
- Scr Trpn vs Trpp
- Mrpl Trpn vs Trpp
- Nduf Trpn vs Trpp


## Plots produced per individual report

- Volcano plot
- GSEA
- ORA
- Cnet plots
- Volcano plots with specific signature genes highlighted

## Summary report plots
### For sample quality control

- MAplot
- Padj histogram
- PCA
- Similarity matrix

### For analysis

- Volcano plots
- Venn diagram for common gene differential expression
- GSEA and ORA with filtered terms of interest
- cnet plots
- Heatmap with z-score for metabolic signature
- Volcano plots with the genes from the signatures highlighted





















