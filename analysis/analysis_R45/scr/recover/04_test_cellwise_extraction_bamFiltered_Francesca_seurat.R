# AIM ---------------------------------------------------------------------
# try to read in the table of gene per barcode extracted from the filtered bam file.
# try to reproduce the counting produced from the official routine suggested by 10X
# this script is focussed on reading in all the cell barcodes

# libraries ---------------------------------------------------------------
library(data.table)
library(tidyverse)
library(Seurat)

# read in the files -------------------------------------------------------
LUT_gene_fix <- read_tsv("../data/LUT_gene_fix_Human.tsv")

# wrangling ---------------------------------------------------------------
# locate the files from the analysis of Francesca
files <- c("K16","S1","Unsorted")

# file_id <- "K16"
# loop the processing and saving of the full matrices
lapply(files,function(file_id){
  # track the processing
  print(file_id)
  
  # determine the file to be loaded
  file <- paste0("/beegfs/scratch/ric.cosr/ric.cosr/CRtest_genova_pedrini/cellranger/cellranger7/output/",
                 file_id,
                 "/outs/summarized_feature_test_all_cellBarcodes.txt")
  file
  # Use fread to read the file and split it into three columns
  df_test4_full <- fread(file, sep = " ", header = FALSE)
  
  # Rename columns for convenience
  setnames(df_test4_full, c("V1", "V2", "V3"), c("count", "gene", "barcode"))
  
  # View the first few rows to verify
  head(df_test4_full)
  dim(df_test4_full)
  class(df_test4_full$count)
  
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
  
  # -------------------------------------------------------------------------
  # save the cell matrix
  saveRDS(sparse_matrix,paste0("../out/object/sparse_matrix_all_",file_id,"_seurat.rds"))
})


# test --------------------------------------------------------------------
# compare the extracted matrix from the calculated matrix
# read in the sparse matrices
list_sparse_matrix <- lapply(files,function(file_id){
  # print(x)
  sparse_matrix <- readRDS(paste0("../out/object/sparse_matrix_all_",file_id,"_seurat.rds"))
  return(sparse_matrix)
}) %>%
  setNames(files)

# read in the matrices from the specific pbjects
list_sparse_matrix2 <- lapply(files,function(file_id){
  # print(x)
  mtx_raw <- Read10X(paste0("/beegfs/scratch/ric.cosr/ric.cosr/CRtest_genova_pedrini/cellranger/cellranger7/output/",
                            file_id,
                            "/outs/raw_feature_bc_matrix/"))
  return(mtx_raw)
}) %>%
  setNames(files)

# read in the barcodes
list_barcodes <- lapply(files,function(file_id){
  # print(x)
  df_barcodes <- read_tsv(paste0("/beegfs/scratch/ric.cosr/ric.cosr/CRtest_genova_pedrini/cellranger/cellranger7/output/",
                                 file_id,
                                 "/outs/cell_barcodes_common_seurat.txt"),
                          col_names = F) %>%
    # remove the flag from the barcodes
    mutate(barcode = str_remove_all(X1,"CB:Z:"))
  
  return(df_barcodes)
}) %>%
  setNames(files)

# make the comparison per dataset
# mat <- list_sparse_matrix$K16
# mat2 <- list_sparse_matrix2$K16
# barc <- list_barcodes$K16
list_summary <- pmap(list(list_sparse_matrix,list_sparse_matrix2,list_barcodes),function(mat,mat2,barc){
  # subset only the barcodes of interest from the cellranger output
  mtx_raw_subset <- mat2[,barc$barcode]
  # now the matrix contains all the barcodes
  mtx_extracted_subset <- mat[,barc$barcode]
  
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
  
  # df_summary_tot %>%
  #   ggplot(aes(x=count_cellranger,y=count_extracted)) + geom_point(alpha = 0.1) + geom_abline(slope = 1,intercept = 0,color = "red",linetype = "dashed") + theme_bw()+theme()+
  #   scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  #   scale_x_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7))
  return(df_summary_tot)
})

# confirm all the deltas are different from 0
lapply(list_summary, function(x){
  x %>% filter(delta != 0)
})
