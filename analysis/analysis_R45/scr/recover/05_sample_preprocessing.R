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

# define the filtering parameters -----------------------------------------
# in this case the matrices already contains all the right cells for the analysis, therefore I am not doing any filtering
# featureLow_thr <- 1000
# featureHigh_thr <- 6000
# mito_thr <- 15
# label <- "01000_06000_15_V5"

# read in the data --------------------------------------------------------
# # location of all the raw matrices
# id_sample <- dir("../data/test_introns/SoupX_default_cellranger7/")
# 
# # load the LUT of the dataset
# LUT <- read_csv("../data/test_introns/LUT_samples.csv")
# 
# # load the LUT of the doublet rate estimate
# df_doublet <- read_csv("../data/dublets_rate_2023.csv")

# load the full object analyzed
sobj <- readRDS("/beegfs/scratch/ric.cosr/ric.cosr/CRtest_genova_pedrini/raw_data/seu_obj/list_datasc_fix_filter_norm_doubletSinglet_SoupX_01000_06000_15_V5.rds")$sample_untreated_Wintron
# add the cluster annotation to the final object
# sobj@meta.data <- sobj@meta.data %>%
#   mutate(cell_id = case_when(seurat_clusters %in% c(6,3,11,15,19)~"ASTRO",
#                              seurat_clusters %in% c(2,5,9,10,14,17,18)~"NEU",
#                              seurat_clusters %in% c(0)~"OPC",
#                              seurat_clusters %in% c(8)~"PROG",
#                              seurat_clusters %in% c(4)~"OLIGO",
#                              seurat_clusters %in% c(12)~"MG",
#                              seurat_clusters %in% c(1,7,16,20)~"GLIA_IMM"))
# confirm the annotation
# DimPlot(sobj,group.by = "cell_id",label = T)

# extract the metadata
meta_barcode <- sobj@meta.data %>%
  rownames_to_column("barcode") %>%
  # dplyr::select(barcode,treat,doxy,exposure,cell_id)
  dplyr::select(barcode,treat,doxy,exposure)

dim(meta_barcode)

# load the proportion of intronic reads
meta_intron <- read_tsv("../out/table/05_prop_intronic_BSedo_seurat.tsv")

dim(meta_intron)

# identify the samples
id_sample <- dir("../out/object/") %>%
  str_subset(pattern = "sparse_") %>%
  str_subset(pattern = "BSedo") %>%
  str_subset(pattern = "seurat")

# wraingling --------------------------------------------------------------
# define the sample LUT
LUT <- data.frame(sample = id_sample) %>%
  mutate(test = str_extract(sample,"all|exon|intron")) %>%
  mutate(dataset = paste0("BSedo_",test))

# merge the prop of intron and the barcode annotation
meta_full <- left_join(meta_barcode,meta_intron,by="barcode")

dim(meta_full)

# run the processing ------------------------------------------------------
# x <- "sparse_matrix_all.rds"
# do the preprocessing over all the dataset and save the objects
list_datasc <- lapply(id_sample,function(x){
  # to track the processing of the progress of the lapply
  print(x)
  
  # read in the matrix
  data <- readRDS(paste0("../out/object/",x))
  
  # crete the object
  # in this case do not apply any filtering, keep all the cells and genes to make a full comparison
  datasc <- CreateSeuratObject(counts = data, project = LUT %>%
                                 filter(sample == x) %>%
                                 pull(dataset), min.cells = 0, min.features = 0)
  
  # add the metadata
  datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^MT-")
  datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  # add also the percentage of globin. in this dataset it is not meaningful as there is no blood
  datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^HB[^(P)]")
  
  # label the cells based on the mt reads content
  datasc$mt_bin <- datasc@meta.data %>%
    mutate(test = case_when(percent.mt < 1~"low",
                            percent.mt < 10~"mid",
                            T ~ "high")) %>%
    pull(test)
  
  # update the meta 
  datasc@meta.data <- datasc@meta.data %>%
    rownames_to_column("barcode") %>%
    left_join(meta_full,by = "barcode") %>%
    column_to_rownames("barcode")
  
  # datasc$treat <- LUT %>%
  #   filter(sample == x) %>%
  #   pull(treat)
  # 
  # datasc$test_intron <- LUT %>%
  #   filter(sample == x) %>%
  #   pull(test_intron)
  # 
  # datasc$doxy <- LUT %>%
  #   filter(sample == x) %>%
  #   pull(doxy)
  # 
  # datasc$exposure <- LUT %>%
  #   filter(sample == x) %>%
  #   pull(exposure)
  # 
  # datasc$ID <- LUT %>%
  #   filter(sample == x) %>%
  #   pull(ID)
  # 
  # # add the filtering variable based on the fixed threshold
  # datasc$test <- datasc@meta.data %>%
  #   mutate(test = percent.mt < mito_thr & nFeature_RNA > featureLow_thr & nFeature_RNA < featureHigh_thr) %>% 
  #   pull(test)
  # 
  # # add the filtering variable based on the
  # stats <- cbind(log10(datasc@meta.data$nCount_RNA), log10(datasc@meta.data$nFeature_RNA),
  #                datasc@meta.data$percent.mt)
  # 
  # # add the filtering variable based on the adaptive threshold
  # # library(robustbase)
  # # library(scater)
  # stats <- cbind(log10(datasc@meta.data$nCount_RNA),
  #                log10(datasc@meta.data$nFeature_RNA),
  #                datasc@meta.data$percent.mt)
  # 
  # outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
  # multi.outlier <- isOutlier(outlying, type = "higher")
  # 
  # datasc$not_outlier <- !as.vector(multi.outlier)
  
  return(datasc)
}) %>%
  setNames(id_sample %>% str_remove(".rds"))

# confirm the class of the objects ----------------------------------------
lapply(list_datasc, function(x){
  class(x@assays$RNA)
})

# save the full metadata --------------------------------------------------
meta_total <- lapply(list_datasc, function(x){
  x@meta.data %>%
    rownames_to_column("barcode") %>%
    mutate(barcode = paste0(barcode,"|",orig.ident))
}) %>%
  bind_rows(.id = "dataset")

# sample standard preprocessing
# x <- list_datasc$sparse_matrix_all
list_datasc_norm <- lapply(list_datasc, function(x){
  
  # add the cell cycle analysis after normalization of the data
  DefaultAssay(datasc_filter) <- "RNA"
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  datasc_filter <- x %>%
    NormalizeData() %>%
    CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
    # I can scale the missing features afterwards now focus on the highly variable one for speed purposes
    ScaleData(vars.to.regress = c("percent.mt","nCount_RNA"), verbose = T) %>% 
    # run this if you want to scale all the variables
    # ScaleData(vars.to.regress = c("percent.mt.harmony","nCount_RNA.harmony","S.Score.harmony","G2M.Score.harmony"), verbose = T,features = all.genes) %>% 
    RunPCA(npcs = 30, verbose = T) %>% 
    RunUMAP(reduction = "pca", dims = 1:30,return.model = TRUE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, by = 0.1)) %>%
    identity()
  
  return(datasc_filter)
})

# save the list of individual objects
saveRDS(list_datasc_norm,"../out/object/05_list_datasc_norm.rds")

# DimPlot(datasc_filter,group.by = "cell_id")
# 
# FeaturePlot(datasc_filter,features = "prop_intronic") + scale_color_viridis_c(option = "turbo")
# FeaturePlot(datasc_filter,features = "nFeature_RNA") + scale_color_viridis_c(option = "turbo")
