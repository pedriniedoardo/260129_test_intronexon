# AIM ---------------------------------------------------------------------
# try to read in the table of gene per barcode extracted from the intronic and exon bam file.
# try to reproduce the counting produced from the official routine suggested by 10X
# this script is focussed on reading in all the barcodes

# libraries ---------------------------------------------------------------
library(data.table)
library(tidyverse)
library(Seurat)
library(ggExtra)
library(ggside)

# read in the files -------------------------------------------------------
# read in the LUT for the genes
LUT_gene_fix <- read_tsv("../data/LUT_gene_fix_Human.tsv")

# Load the data.table package
# library(data.table)

# Use fread to read the file and split it into three columns
df_test4_intron <- fread("/beegfs/scratch/ric.absinta/ric.absinta/analysis/BS_run05/test/Sample_W8_untreated_Wintron_test/summarized_feature_test_intron_fullBarcodes.txt", sep = " ", header = FALSE)
df_test4_exon <- fread("/beegfs/scratch/ric.absinta/ric.absinta/analysis/BS_run05/test/Sample_W8_untreated_Wintron_test/summarized_feature_test_exon_fullBarcodes.txt", sep = " ", header = FALSE)

# Rename columns for convenience
setnames(df_test4_intron, c("V1", "V2", "V3"), c("count", "gene", "barcode"))
setnames(df_test4_exon, c("V1", "V2", "V3"), c("count", "gene", "barcode"))

# View the first few rows to verify
head(df_test4_intron)
dim(df_test4_intron)
class(df_test4_intron$count)

head(df_test4_exon)
dim(df_test4_exon)
class(df_test4_exon$count)

# read in the full cell barcode
df_barcodes <- read_tsv("/beegfs/scratch/ric.absinta/ric.absinta/analysis/BS_run05/test/Sample_W8_untreated_Wintron_test/cell_barcodes_common_seurat.txt",
                        col_names = F) %>%
  # remove the flag from the barcodes
  mutate(barcode = str_remove_all(X1,"CB:Z:"))


dim(df_barcodes)

# read in the raw unfiltered table of counts that is the output of cellranger
mtx_raw <- Read10X("/beegfs/scratch/ric.absinta/ric.absinta/analysis/BS_run05/test/Sample_W8_untreated_Wintron/outs/raw_feature_bc_matrix/")

# test sparse matrix ------------------------------------------------------
# full matrix
test_sparse_intron <- df_test4_intron
test_sparse_exon <- df_test4_exon

# Convert GeneID and Barcode to factors for row and column indices
test_sparse_intron <- test_sparse_intron %>%
  mutate(gene = factor(gene,levels = LUT_gene_fix$gene_id2),
         barcode = as.factor(barcode))

test_sparse_exon <- test_sparse_exon %>%
  mutate(gene = factor(gene,levels = LUT_gene_fix$gene_id2),
         barcode = as.factor(barcode))

# Create the sparse matrix using sparseMatrix
# library(Matrix)
sparse_matrix_exon <- Matrix::sparseMatrix(
  i = as.integer(test_sparse_exon$gene),       # row indices from GeneID factor levels
  j = as.integer(test_sparse_exon$barcode),      # column indices from Barcode factor levels
  x = test_sparse_exon$count,                   # values to fill in
  dims = c(length(levels(test_sparse_exon$gene)), length(levels(test_sparse_exon$barcode))),
  dimnames = list(levels(test_sparse_exon$gene), levels(test_sparse_exon$barcode))
)

sparse_matrix_intron <- Matrix::sparseMatrix(
  i = as.integer(test_sparse_intron$gene),       # row indices from GeneID factor levels
  j = as.integer(test_sparse_intron$barcode),      # column indices from Barcode factor levels
  x = test_sparse_intron$count,                   # values to fill in
  dims = c(length(levels(test_sparse_intron$gene)), length(levels(test_sparse_intron$barcode))),
  dimnames = list(levels(test_sparse_intron$gene), levels(test_sparse_intron$barcode))
)

# fix the column name to match the mtx nomenclature
colnames(sparse_matrix_exon) <- str_remove_all(colnames(sparse_matrix_exon),pattern = "CB:Z:")
colnames(sparse_matrix_intron) <- str_remove_all(colnames(sparse_matrix_intron),pattern = "CB:Z:")

# fix the row names
rownames(sparse_matrix_exon) <- data.frame(sparse_row = rownames(sparse_matrix_exon)) %>%
  left_join(LUT_gene_fix,by = c("sparse_row" = "gene_id2")) %>%
  pull(gene_name2)

rownames(sparse_matrix_intron) <- data.frame(sparse_row = rownames(sparse_matrix_intron)) %>%
  left_join(LUT_gene_fix,by = c("sparse_row" = "gene_id2")) %>%
  pull(gene_name2)

# View the sparse matrix
print(sparse_matrix_exon)
dim(sparse_matrix_exon)

print(sparse_matrix_intron)
dim(sparse_matrix_intron)

# test --------------------------------------------------------------------
# compare the two matrices
# subset only the barcodes of interest from the cellranger output
mtx_raw_subset <- mtx_raw[,df_barcodes$barcode]
mtx_extractExon_subset <- sparse_matrix_exon[,df_barcodes$barcode]
mtx_extractIntron_subset <- sparse_matrix_intron[,df_barcodes$barcode]

# now the two matrices should be equal
dim(mtx_raw_subset)
dim(sparse_matrix_exon)
dim(mtx_extractExon_subset)
dim(sparse_matrix_intron)
dim(mtx_extractIntron_subset)

# the intron and exon matrices are already matched by barcodes and features therefore they can be summed
sparse_matrix_all <- mtx_extractExon_subset + mtx_extractIntron_subset

# check in they are identical
identical(mtx_raw_subset,sparse_matrix_all)

# compare the summaries
df_cellranger <- rowSums(mtx_raw_subset) %>%
  data.frame() %>%
  setNames("count_cellranger") %>%
  rownames_to_column("gene")

df_extracted <- rowSums(sparse_matrix_all) %>%
  data.frame() %>%
  setNames("count_extracted") %>%
  rownames_to_column("gene")

df_summary_tot <- left_join(df_cellranger,df_extracted,by = "gene") %>%
  mutate(delta = count_cellranger - count_extracted)

df_summary_tot %>%
  ggplot(aes(x=count_cellranger,y=count_extracted)) + geom_point(alpha = 0.1) + geom_abline(slope = 1,intercept = 0,color = "red",linetype = "dashed") + theme_bw()+theme()+
  scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  scale_x_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7))

# plot the summarised exon vs intron reads
df_extracted_intron <- rowSums(mtx_extractIntron_subset) %>%
  data.frame() %>%
  setNames("count_intron") %>%
  rownames_to_column("gene")

df_extracted_exon <- rowSums(mtx_extractExon_subset) %>%
  data.frame() %>%
  setNames("count_exon") %>%
  rownames_to_column("gene")

df_summary_tot2 <- left_join(df_extracted_exon,df_extracted_intron,by = "gene")

p <- df_summary_tot2 %>%
  ggplot(aes(x=count_intron,y=count_exon)) +
  geom_point(alpha = 0.1) +
  geom_abline(slope = 1,intercept = 0,color = "red",linetype = "dashed") +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme() +
  scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  scale_x_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  coord_fixed()+
  geom_vline(xintercept = mean(df_summary_tot2$count_intron),linetype = "dashed",col="gray")+
  geom_hline(yintercept = mean(df_summary_tot2$count_exon),linetype = "dashed",col="gray")

ggMarginal(p, type = "histogram", margins = "both", size = 5, color = "black")

# -------------------------------------------------------------------------
# save the cell matrix
saveRDS(sparse_matrix_exon,"../out/object/sparse_matrix_exon_BSedo_fullBarcodes.rds")
saveRDS(sparse_matrix_intron,"../out/object/sparse_matrix_intron_BSedo_fullBarcodes.rds")
