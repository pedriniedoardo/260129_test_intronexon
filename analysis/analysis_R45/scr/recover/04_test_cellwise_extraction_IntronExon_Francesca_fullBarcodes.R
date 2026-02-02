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
  file_intron <- paste0("/beegfs/scratch/ric.cosr/ric.cosr/CRtest_genova_pedrini/cellranger/cellranger7/output/",
                        file_id,
                        "/outs/summarized_feature_test_intron_fullBarcodes.txt")
  
  file_exon <- paste0("/beegfs/scratch/ric.cosr/ric.cosr/CRtest_genova_pedrini/cellranger/cellranger7/output/",
                      file_id,
                      "/outs/summarized_feature_test_exon_fullBarcodes.txt")
  
  # Use fread to read the file and split it into three columns
  df_test4_intron <- fread(file_intron, sep = " ", header = FALSE)
  df_test4_exon <- fread(file_exon, sep = " ", header = FALSE)
  
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
  
  # -------------------------------------------------------------------------
  # save the cell matrix
  saveRDS(sparse_matrix_exon,paste0("../out/object/sparse_matrix_exon_",file_id,"_fullBarcodes.rds"))
  saveRDS(sparse_matrix_intron,paste0("../out/object/sparse_matrix_intron_",file_id,"_fullBarcodes.rds"))
})


# test --------------------------------------------------------------------
# compare the extracted matrix from the calculated matrix
# read in the sparse matrices from intronic reads and exonic reads
list_sparse_matrix <- lapply(files,function(file_id){
  # print(x)
  sparse_matrix_intron <- readRDS(paste0("../out/object/sparse_matrix_intron_",file_id,"_fullBarcodes.rds"))
  sparse_matrix_exon <- readRDS(paste0("../out/object/sparse_matrix_exon_",file_id,"_fullBarcodes.rds"))
  
  df_barcodes <- read_tsv(paste0("/beegfs/scratch/ric.cosr/ric.cosr/CRtest_genova_pedrini/cellranger/cellranger7/output/",
                                 file_id,
                                 "/outs/cell_barcodes_common_seurat.txt"),
                          col_names = F) %>%
    # remove the flag from the barcodes
    mutate(barcode = str_remove_all(X1,"CB:Z:"))
  
  mtx_extractExon_subset <- sparse_matrix_exon[,df_barcodes$barcode]
  mtx_extractIntron_subset <- sparse_matrix_intron[,df_barcodes$barcode]
  
  # now the two matrices should be equal
  dim(sparse_matrix_exon)
  dim(mtx_extractExon_subset)
  dim(sparse_matrix_intron)
  dim(mtx_extractIntron_subset)
  
  # the intron and exon matrices are already matched by barcodes and features therefore they can be summed
  sparse_matrix_all <- mtx_extractExon_subset + mtx_extractIntron_subset
  
  return(sparse_matrix_all)
}) %>%
  setNames(files)

# read in the matrices from the specific pbjects
list_sparse_matrix2 <- lapply(files,function(file_id){
  # print(x)
  mtx_raw <- Read10X(paste0("/beegfs/scratch/ric.cosr/ric.cosr/CRtest_genova_pedrini/cellranger/cellranger7/output/",
                            file_id,
                            "/outs/raw_feature_bc_matrix/"))
  
  df_barcodes <- read_tsv(paste0("/beegfs/scratch/ric.cosr/ric.cosr/CRtest_genova_pedrini/cellranger/cellranger7/output/",
                                 file_id,
                                 "/outs/cell_barcodes_common_seurat.txt"),
                          col_names = F) %>%
    # remove the flag from the barcodes
    mutate(barcode = str_remove_all(X1,"CB:Z:"))
  
  mtx_raw_subset <- mtx_raw[,df_barcodes$barcode]
  
  return(mtx_raw_subset)
}) %>%
  setNames(files)

# make the comparison per dataset
# mat <- list_sparse_matrix$K16
# mat2 <- list_sparse_matrix2$K16
list_summary <- pmap(list(list_sparse_matrix,list_sparse_matrix2),function(mat,mat2){

  # now the two matrices should be equal
  dim(mat)
  dim(mat2)
  
  identical(mat,mat2)
  
  # compare the summaries
  df_cellranger <- rowSums(mat2) %>%
    data.frame() %>%
    setNames("count_cellranger") %>%
    rownames_to_column("gene")
  
  df_extracted <- rowSums(mat) %>%
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
