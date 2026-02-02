# AIM ---------------------------------------------------------------------
# try to read in the table of gene per barcode extracted from the filtered bam file.
# try to reproduce the counting produced from the official routine suggested by 10X

# libraries ---------------------------------------------------------------
library(data.table)
library(tidyverse)
library(Seurat)

# read in the files -------------------------------------------------------
# read in the LUT for the geneID conversion
# "/beegfs/scratch/ric.cosr/pedrini.edoardo/reference/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
# LUT_gene <- rtracklayer::import("/beegfs/scratch/ric.cosr/pedrini.edoardo/reference/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf") %>%
#   as.data.frame() %>%
#   dplyr::filter(type == "gene")
# 
# LUT_gene_fix <- LUT_gene %>%
#   # add the rowids compatible with the count matrix
#   mutate(gene_name2 = rownames(mtx_raw)) %>%
#   # make the gene_id compatible with the count matrix
#   mutate(gene_id2 = paste0("GX:Z:",gene_id)) %>%
#   # remove the column that do not bring new information
#   select(-c(transcript_id:ont))
# 
# # sanity check, the difference between gene_name and gene_name2 should be on a few genes because of the dot notation of the redoundant ids
# LUT_gene_fix %>%
#   mutate(test = gene_name == gene_name2) %>%
#   filter(test == F)
# 
# # save the object
# LUT_gene_fix %>%
#   write_tsv("../data/LUT_gene_fix_Human.tsv")
LUT_gene_fix <- read_tsv("../data/LUT_gene_fix_Human.tsv")

# Load the data.table package
# library(data.table)

# Use fread to read the file and split it into three columns
df_test4_full <- fread("/beegfs/scratch/ric.absinta/ric.absinta/analysis/BS_run05/test/Sample_W8_untreated_Wintron_test/summarized_feature_test_all_cellBarcodes.txt", sep = " ", header = FALSE)

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


dim(df_barcodes)

# read in the raw unfiltered table of counts that is the output of cellranger
mtx_raw <- Read10X("/beegfs/scratch/ric.absinta/ric.absinta/analysis/BS_run05/test/Sample_W8_untreated_Wintron/outs/raw_feature_bc_matrix/")

# wrangling ---------------------------------------------------------------
# # select just a few cells to start the investigation
# df_test4 <- df_test4_full %>%
#   filter(barcode %in% c("CB:Z:TCATGCCGTTTGAACC-1","CB:Z:AAACCCAAGCTATCCA-1"))
# 
# # add the gene informations to the subset of cells
# df_test_bamExtr <- df_test4 %>%
#   mutate(enembl_gene_id = str_remove_all(gene,pattern="GX:Z:")) %>%
#   left_join(LUT_gene_full,"enembl_gene_id") %>%
#   select(count_extracted = count,enembl_gene_id,symbol,barcode) %>%
#   mutate(cell_barcode = str_remove_all(barcode,"CB:Z:"))
# 
# # subset only the barcodes of interest from the cellranger output
# mtx_raw_subset <- mtx_raw[,c("TCATGCCGTTTGAACC-1","AAACCCAAGCTATCCA-1")]
# dim(mtx_raw_subset)
# 
# # summarised the reads per gene
# feature_allCellranger_test <- mtx_raw_subset %>%
#   data.frame() %>%
#   rownames_to_column("symbol_cellranger") %>%
#   # add the gene ids defined by the feature.csv file the order is mainteined
#   mutate(enembl_gene_id = LUT_gene_full$enembl_gene_id,
#          symbol = LUT_gene_full$symbol) %>%
#   pivot_longer(names_to = "cell_barcode",
#                values_to = "count_cellranger",-c(symbol_cellranger,enembl_gene_id,symbol)) %>%
#   mutate(cell_barcode = str_replace_all(cell_barcode,"\\.","-"))
# 
# # join the two table to confirm that the counts are the same
# test_compare_cell <- feature_allCellranger_test %>%
#   select(enembl_gene_id,cell_barcode,count_cellranger,symbol) %>%
#   left_join(df_test_bamExtr %>% select(-c(barcode,symbol)),by = c("cell_barcode","enembl_gene_id")) %>%
#   # convert to 0 the counts from genes that are not
#   mutate(count_extracted = case_when(is.na(count_extracted)~0,
#                                      T~count_extracted)) %>% 
#   mutate(delta = count_cellranger - count_extracted)
# 
# # all the genes have a delta of 0
# test_compare_cell %>%
#   filter(delta != 0)

# test sparse matrix ------------------------------------------------------
# full matrix
test_sparse <- df_test4_full

# # Convert GeneID and Barcode to factors for row and column indices
# test_sparse[, gene := as.factor(gene)]
# test_sparse[, barcode := as.factor(barcode)]
# 
# # Create the sparse matrix using sparseMatrix
# # library(Matrix)
# sparse_matrix <- Matrix::sparseMatrix(
#   i = as.integer(test_sparse$gene),       # row indices from GeneID factor levels
#   j = as.integer(test_sparse$barcode),      # column indices from Barcode factor levels
#   x = test_sparse$count,                   # values to fill in
#   dims = c(length(levels(test_sparse$gene)), length(levels(test_sparse$barcode))),
#   dimnames = list(levels(test_sparse$gene), levels(test_sparse$barcode))
# )

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

# now the two matrices should be equal
dim(mtx_raw_subset)
dim(sparse_matrix)

identical(mtx_raw_subset,sparse_matrix)

# compare the summaries
df_cellranger <- rowSums(mtx_raw_subset) %>%
  data.frame() %>%
  setNames("count_cellranger") %>%
  rownames_to_column("gene")

df_extracted <- rowSums(sparse_matrix) %>%
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
saveRDS(sparse_matrix,"../out/object/sparse_matrix_all_BSedo_seurat.rds")
