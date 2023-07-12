# Merge .csv count files

## Assumes you have received data from Raphaël Helaers
# And you have an output like : "_ReadsPerGene.out.tab"

library(DESeq2)
library(tidyverse)
library(stringr)
library(biomaRt)
library(Hmisc)

names = c("Mrpl_Trpn_1","Mrpl_Trpn_2","Mrpl_Trpn_3",
          "Mrpl_Trpp_1","Mrpl_Trpp_2","Mrpl_Trpp_3",
          "Nduf_Trpn_1","Nduf_Trpn_2","Nduf_Trpn_3",
          "Nduf_Trpp_1","Nduf_Trpp_2","Nduf_Trpp_3",
          "Nonuc_Trpn_1","Nonuc_Trpn_2","Nonuc_Trpn_3",
          "Nonuc_Trpp_1","Nonuc_Trpp_2","Nonuc_Trpp_3",
          "Scr_Trpn_1","Scr_Trpn_2","Scr_Trpn_3",
          "Scr_Trpp_1","Scr_Trpp_2","Scr_Trpp_3")

files = list.files("./data/", pattern = ".out.tab")
files # check that the order of names and files is correct or you will 
      # assign wrong names to samples

for(i in 1:length(files)){
  x = read_tsv(paste0("./data/",files[i]), col_names = F)
  x = x[-c(1:5),-c(3,4)]
  colnames(x) = c("ensembl",paste0(names[i]))
  assign(names[i],x)
}

noquote(names) # Use that to quickly copy paste the names without the quotes

counts <- list(Mrpl_Trpn_1,  Mrpl_Trpn_2,  Mrpl_Trpn_3,
               Mrpl_Trpp_1,  Mrpl_Trpp_2,  Mrpl_Trpp_3,
               Nduf_Trpn_1,  Nduf_Trpn_2,  Nduf_Trpn_3,
               Nduf_Trpp_1,  Nduf_Trpp_2,  Nduf_Trpp_3,
               Nonuc_Trpn_1, Nonuc_Trpn_2, Nonuc_Trpn_3,
               Nonuc_Trpp_1, Nonuc_Trpp_2, Nonuc_Trpp_3,
               Scr_Trpn_1,  Scr_Trpn_2,  Scr_Trpn_3,
               Scr_Trpp_1,  Scr_Trpp_2,  Scr_Trpp_3 ) %>% 
  purrr::reduce(full_join, by = "ensembl")


strrep =
  sub(pattern = "\\.(.*)","",counts$ensembl)

counts$ensembl = strrep
counts = counts[!duplicated(counts$ensembl),]

rownames(counts) = counts$ensembl

write.csv(counts,"./data/RAB234_counts.csv", row.names = T)

# Create the coldata for the summarized experiment

coldata <- data.frame(
  celltype =c(rep(c("Mrpl","Nduf","Nonuc","Scr"),each = 6)),
  condition = as.factor(rep(c("Trpn","Trpp"), times = 4, each = 3)),
  replicate=as.factor(rep(c(1:3),8)))

rownames(coldata) <- colnames(counts)[-1]

write.csv(coldata,"./data/coldata.csv")

