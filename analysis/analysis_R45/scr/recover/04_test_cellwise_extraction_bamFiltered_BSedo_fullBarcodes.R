# AIM ---------------------------------------------------------------------
# try to read in the table of gene per barcode extracted from the filtered bam file.
# try to reproduce the counting produced from the official routine suggested by 10X
# this script is focussed on reading in all the barcodes

# libraries ---------------------------------------------------------------
library(data.table)
library(tidyverse)
library(Seurat)

# read in the files -------------------------------------------------------
LUT_gene_fix <- read_tsv("../data/LUT_gene_fix_Human.tsv")

# Use fread to read the file and split it into three columns
df_test4_full <- fread("/beegfs/scratch/ric.absinta/ric.absinta/analysis/BS_run05/test/Sample_W8_untreated_Wintron_test/summarized_feature_test_all_fullBarcodes.txt", sep = " ", header = FALSE)

# Rename columns for convenience
setnames(df_test4_full, c("V1", "V2", "V3"), c("count", "gene", "barcode"))

# View the first few rows to verify
head(df_test4_full)
dim(df_test4_full)
class(df_test4_full$count)

# read in the full cell barcode
df_barcodes <- read_tsv("/beegfs/scratch/ric.absinta/ric.absinta/analysis/BS_run05/test/Sample_W8_untreated_Wintron_test/cell_barcodes_common_seurat.txt",
                        col_names = F) %>%
  # remove the flag from the barcodes
  mutate(barcode = str_remove_all(X1,"CB:Z:"))


# dim(df_barcodes)

# read in the raw unfiltered table of counts that is the output of cellranger
mtx_raw <- Read10X("/beegfs/scratch/ric.absinta/ric.absinta/analysis/BS_run05/test/Sample_W8_untreated_Wintron/outs/raw_feature_bc_matrix/")

# test sparse matrix ------------------------------------------------------
# full matrix
test_sparse <- df_test4_full

# Convert GeneID and Barcode to factors for row and column indices
test_sparse <- test_sparse %>%
  mutate(gene = factor(gene,levels = LUT_gene_fix$gene_id2),
         barcode = as.factor(barcode))

# Create the sparse matrix using sparseMatrix
# library(Matrix)
sparse_matrix <- Matrix::sparseMatrix(
  i = as.integer(test_sparse$gene),       # row indices from GeneID factor levels
  j = as.integer(test_sparse$barcode),      # column indices from Barcode factor levels
  x = test_sparse$count,                   # values to fill in
  dims = c(length(levels(test_sparse$gene)), length(levels(test_sparse$barcode))),
  dimnames = list(levels(test_sparse$gene), levels(test_sparse$barcode))
)

# fix the column name to match the mtx nomenclature
colnames(sparse_matrix) <- str_remove_all(colnames(sparse_matrix),pattern = "CB:Z:")

# fix the row names
rownames(sparse_matrix) <- data.frame(sparse_row = rownames(sparse_matrix)) %>%
  left_join(LUT_gene_fix,by = c("sparse_row" = "gene_id2")) %>%
  pull(gene_name2)

# View the sparse matrix
print(sparse_matrix)
dim(sparse_matrix)

# test --------------------------------------------------------------------
# compare the two matrices
# subset only the barcodes of interest from the cellranger output
mtx_raw_subset <- mtx_raw[,df_barcodes$barcode]
# now the matrix contains all the barcodes
mtx_extracted_subset <- sparse_matrix[,df_barcodes$barcode]

# now the two matrices should be equal
dim(mtx_raw_subset)
dim(mtx_extracted_subset)

identical(mtx_raw_subset,mtx_extracted_subset)

# compare the summaries
df_cellranger <- rowSums(mtx_raw_subset) %>%
  data.frame() %>%
  setNames("count_cellranger") %>%
  rownames_to_column("gene")

df_extracted <- rowSums(mtx_extracted_subset) %>%
  data.frame() %>%
  setNames("count_extracted") %>%
  rownames_to_column("gene")

df_summary_tot <- left_join(df_cellranger,df_extracted,by = "gene") %>%
  mutate(delta = count_cellranger - count_extracted)

df_summary_tot %>%
  ggplot(aes(x=count_cellranger,y=count_extracted)) + geom_point(alpha = 0.1) + geom_abline(slope = 1,intercept = 0,color = "red",linetype = "dashed") + theme_bw()+theme()+
  scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  scale_x_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7))

# -------------------------------------------------------------------------
# save the cell matrix
saveRDS(sparse_matrix,"../out/object/sparse_matrix_all_BSedo_fullBarcodes.rds")
