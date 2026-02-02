# AIM ---------------------------------------------------------------------
# try to read in the table of gene per barcode extracted from the filtered bam file.
# try to reproduce the counting produced from the official routine suggested by 10X
# this script is focussed on reading in all the cell barcodes

# function generate the gtf ifle ------------------------------------------
ReadGTFIntronExon <- function(gtf_loc = NULL){
  
  # load the required libraries
  suppressMessages(suppressWarnings({
    require(rtracklayer)
    require(tidyverse)
  }))
  
  # "../pipeline_01/output/PureIntronExon/genes.gtf.gz"
  LUT_gene_fix <- rtracklayer::import(gtf_loc) %>%
    as.data.frame() %>%
    dplyr::filter(type == "gene") %>%
    dplyr::select(gene_id,gene_name) %>%
    # make the gene names unique. this is needed for the sparse matrix generation
    mutate(gene_name2 = make.unique(gene_name)) %>%
    mutate(gene_id2 = paste0("GX:Z:",gene_id)) 
  
  # DEBUG: sanity check, the difference between gene_name and gene_name2 should be on a few genes because of the dot notation of the redoundant ids
  # LUT_gene_fix %>%
  #   mutate(test = gene_name == gene_name2) %>%
  #   filter(test == F)
  
  # gene_test <- LUT_gene_fix %>%
  #   group_by(gene_name) %>%
  #   summarise(n = n()) %>%
  #   filter(n > 1)
  # 
  # LUT_gene_fix %>%
  #   filter(gene_name %in% gene_test$gene_name) %>%
  #   arrange(gene_name)
  
  return(LUT_gene_fix)
}

# test
# LUT_gene_fix <- ReadGTFIntronExon(gtf_loc = "../pipeline_01/output/PureIntronExon/genes.gtf.gz")

# function generate the pure sparse matrix --------------------------------
# define the location of the sample.
Read10XIntronExon <- function(file = NULL, gtf_obj = NULL){
  
  require(data.table)
  require(Matrix)
  require(tidyverse)
  
  # accept only one file at a time
  
  # determine the file to be loaded
  file_intron <- paste0(file,"/summarized_feature_test_intron_fullBarcodes.txt")
  file_exon <- paste0(file,"/summarized_feature_test_exon_fullBarcodes.txt")
  
  # Use fread to read the file and split it into three columns
  df_test4_intron <- fread(file_intron, sep = " ", header = FALSE)
  df_test4_exon <- fread(file_exon, sep = " ", header = FALSE)
  
  # Rename columns for convenience
  setnames(df_test4_intron, c("V1", "V2", "V3"), c("count", "gene", "barcode"))
  setnames(df_test4_exon, c("V1", "V2", "V3"), c("count", "gene", "barcode"))
  
  # DEBUG: check the structure of the data
  # head(df_test4_intron)
  # dim(df_test4_intron)
  # class(df_test4_intron$count)
  # 
  # head(df_test4_exon)
  # dim(df_test4_exon)
  # class(df_test4_exon$count)
  
  # test sparse matrix ------------------------------------------------------
  # full matrix
  test_sparse_intron <- df_test4_intron
  test_sparse_exon <- df_test4_exon
  
  # Convert GeneID and Barcode to factors for row and column indices
  test_sparse_intron <- test_sparse_intron %>%
    mutate(gene = factor(gene,levels = gtf_obj$gene_id2),
           barcode = as.factor(barcode))
  
  test_sparse_exon <- test_sparse_exon %>%
    mutate(gene = factor(gene,levels = gtf_obj$gene_id2),
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
    left_join(gtf_obj,by = c("sparse_row" = "gene_id2")) %>%
    pull(gene_name2)
  
  rownames(sparse_matrix_intron) <- data.frame(sparse_row = rownames(sparse_matrix_intron)) %>%
    left_join(gtf_obj,by = c("sparse_row" = "gene_id2")) %>%
    pull(gene_name2)
  
  # DEBUG: View the sparse matrix
  # print(sparse_matrix_exon)
  # dim(sparse_matrix_exon)
  # 
  # print(sparse_matrix_intron)
  # dim(sparse_matrix_intron)
  
  # -------------------------------------------------------------------------
  # save the cell matrix
  # saveRDS(sparse_matrix_exon,paste0("../out/object/sparse_matrix_exon_",file_id,"_fullBarcodes.rds"))
  # saveRDS(sparse_matrix_intron,paste0("../out/object/sparse_matrix_intron_",file_id,"_fullBarcodes.rds"))
  list_pure <- list(sparse_matrix_exon = sparse_matrix_exon,
                    sparse_matrix_intron = sparse_matrix_intron)
  
  return(list_pure)
  
}

# test
# file_loc <- c("../pipeline_01/output/PureIntronExon/tinygex")
# list_pure <- Read10XIntronExon(file = file_loc, gtf_obj = LUT_gene_fix)

# proportion intronic per gene --------------------------------------------
PropIntronicPerGene <- function(sparse_matrix_exon = NULL, sparse_matrix_intron = NULL){
  
  require(Matrix)
  require(tidyverse)
  
  # sum the intronic and exonic counts per gene
  df_summary_exon <- Matrix::rowSums(sparse_matrix_exon) %>%
    data.frame() %>%
    setNames("count_exon") %>%
    rownames_to_column("gene")
  
  df_summary_intron <- Matrix::rowSums(sparse_matrix_intron) %>%
    data.frame() %>%
    setNames("count_intron") %>%
    rownames_to_column("gene")
  
  # join the two dataframes
  df_summary_tot <- left_join(df_summary_exon,df_summary_intron,by = "gene") %>%
    mutate(across(everything(), ~ replace_na(.x, 0)))
  
  # calculate the proportion of intronic reads
  df_summary_tot <- df_summary_tot %>%
    mutate(prop_intron = count_intron/(count_intron + count_exon))
  
  return(df_summary_tot)
}
# test
# df_summary_tot <- PropIntronicPerGene(sparse_matrix_exon = list_pure$sparse_matrix_exon, sparse_matrix_intron = list_pure$sparse_matrix_intron)


# proportion intronic per cell --------------------------------------------
PropIntronicPerCell <- function(sparse_matrix_exon = NULL, sparse_matrix_intron = NULL){
  
  require(Matrix)
  require(tidyverse)
  
  # sum the intronic and exonic counts per gene
  df_summary_exon <- Matrix::colSums(sparse_matrix_exon) %>%
    data.frame() %>%
    setNames("count_exon") %>%
    rownames_to_column("barcode")
  
  df_summary_intron <- Matrix::colSums(sparse_matrix_intron) %>%
    data.frame() %>%
    setNames("count_intron") %>%
    rownames_to_column("barcode")
  
  # join the two dataframes. if one cell is missin in either one or the other summary, make it as a zero
  df_summary_tot <- full_join(df_summary_exon,df_summary_intron,by = "barcode") %>%
    mutate(across(everything(), ~ replace_na(.x, 0)))
  
  # calculate the proportion of intronic reads
  df_summary_tot <- df_summary_tot %>%
    mutate(prop_intron = count_intron/(count_intron + count_exon))
  
  return(df_summary_tot)
}

# test
# df_summary_tot <- PropIntronicPerCell(sparse_matrix_exon = list_pure$sparse_matrix_exon, sparse_matrix_intron = list_pure$sparse_matrix_intron)