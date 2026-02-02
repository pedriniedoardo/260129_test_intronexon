# AIM ---------------------------------------------------------------------
# try to read the data from intronic, exonic and all to be used used as a regular object

# libraries ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)
library(DoubletFinder)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# read in the data --------------------------------------------------------
# read in the cell-wise matrix for the intronic only reads
data_intronic <- readRDS(paste0("../out/object/sparse_matrix_intron_BSedo_seurat.rds"))
data_all <- readRDS(paste0("../out/object/sparse_matrix_all_BSedo_seurat.rds"))

# wrangling ---------------------------------------------------------------
# make a summary per cell of all the reads
df_all <- data.frame(colSums(data_all)) %>%
  setNames("count_all") %>%
  rownames_to_column("barcode")

df_intronic <- data.frame(colSums(data_intronic)) %>%
  setNames("count_intronic") %>%
  rownames_to_column("barcode")

df_all_full <- left_join(df_all,df_intronic,by ="barcode") %>%
  mutate(prop_intronic = count_intronic/count_all)

# save the summary table
df_all_full %>%
  write_tsv("../out/table/05_prop_intronic_BSedo_seurat.tsv")

# plot --------------------------------------------------------------------
df_all_full %>%
  ggplot(aes(x=prop_intronic))+geom_histogram()+theme_bw()
