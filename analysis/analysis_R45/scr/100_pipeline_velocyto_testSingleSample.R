# AIM ---------------------------------------------------------------------
# test how to work with the outputs of velocyto
# this script will process a test sample

# libraries ---------------------------------------------------------------
library(Seurat)
library(patchwork)
library(ggExtra)
library(GGally)
library(SeuratWrappers)
library(velocyto.R)

# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1000 * 1024^2)

# read in the data --------------------------------------------------------
# read in the loom file
ldat <- ReadVelocity(file = "../../data/test_connect5k_PBMC/results/velocyto/merged/connect_5k_pbmc_NGSC3_ch1_gex_1/connect_5k_pbmc_NGSC3_ch1_gex_1.loom")

# confirm all the matrices are present
lapply(ldat,function(x){
  dim(x)
})

# notice that this are already filtered for a subset of barcodes
glimpse(ldat)
ldat$spliced[c("ACTB","GAPDH","TUBB"),1000:1050]

# testing -----------------------------------------------------------------
# read in the regular matrix from the cellranger output
# library(Seurat)
mtx_raw <- Read10X("../../data/test_connect5k_PBMC/results/cellranger/merged/connect_5k_pbmc_NGSC3_ch1_gex_1/outs/filtered_feature_bc_matrix/")                             

# define the filtered barcodes to make a common comparison
cell_barcodes <- colnames(mtx_raw)

# reneme the barcodes to match the one in the origina file
ldat_fix <- lapply(ldat, function(x){
  # update the barcodes
  barcode_new <- data.frame(barcode_old = colnames(x)) %>%
    mutate(barcode_new = str_remove(barcode_old,pattern = "connect_5k_pbmc_NGSC3_ch1_gex_1:") %>% str_replace(pattern = "x$",replacement = "-1")) %>%
    pull(barcode_new)
  colnames(x) <- barcode_new

  return(x)
})

# all the barcode have the tag "-1"
data.frame(cell_barcodes = cell_barcodes) %>%
  separate(col = cell_barcodes,into = c("barcode_id","barcode_tag"),sep = "-") %>%
  group_by(barcode_tag) %>%
  summarise(n = n())

# confirm that the matrices in the loom file are already filtered
lapply(ldat_fix, function(x){
  sum(!colnames(x) %in% cell_barcodes)
})

# generate the combined matrix from the pure matrices
mtx_test <- ldat_fix$spliced + ldat_fix$unspliced + ldat_fix$ambiguous
# mtx_test <- ldat_fix$spliced + ldat_fix$unspliced

# compare the summaries
df_cellranger <- Matrix::rowSums(mtx_raw) %>%
  data.frame() %>%
  setNames("count_cellranger") %>%
  rownames_to_column("gene")

test <- Matrix::rowSums(mtx_test)
df_extracted <- data.frame(count_extracted = test) %>%
  mutate(gene = names(test))

df_summary_tot <- left_join(df_cellranger,df_extracted,by = "gene") %>%
  mutate(delta = count_cellranger - count_extracted)

# plot the difference
df_summary_tot %>%
  ggplot(aes(x=count_cellranger,y=count_extracted)) +
  geom_point(alpha = 0.1) +
  geom_abline(slope = 1,intercept = 0,color = "red",linetype = "dashed") +
  theme_bw() +
  theme() +
  scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  scale_x_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7))

# confirm all the deltas are different from 0
df_summary_tot %>%
  filter(delta != 0)

# plot the comparison of the reads from all or from exonic only
df_summary_pure <- inner_join(
  Matrix::rowSums(ldat_filter$sparse_matrix_exon) %>%
    data.frame() %>%
    setNames("pure_exonic") %>%
    rownames_to_column("gene"),
  Matrix::rowSums(ldat_filter$sparse_matrix_intron) %>%
    data.frame() %>%
    setNames("pure_intronic") %>%
    rownames_to_column("gene"),by = "gene")

df_summary_tot2 <- left_join(df_summary_tot,df_summary_pure,by = "gene")

p1 <- df_summary_tot2 %>%
  ggplot(aes(x=count_cellranger,y=pure_exonic)) +
  geom_point(alpha = 0.1) +
  geom_abline(slope = 1,intercept = 0,color = "red",linetype = "dashed") +
  theme_bw() +
  theme() +
  scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  scale_x_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  coord_fixed() +
  ggtitle("pure exonic vs all")

p2 <- df_summary_tot2 %>%
  ggplot(aes(x=count_cellranger,y=pure_intronic)) +
  geom_point(alpha = 0.1) +
  geom_abline(slope = 1,intercept = 0,color = "red",linetype = "dashed") +
  theme_bw() +
  theme() +
  scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  scale_x_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  coord_fixed() +
  ggtitle("pure intronic vs all")

p3 <- df_summary_tot2 %>%
  ggplot(aes(x=pure_exonic,y=pure_intronic)) +
  geom_point(alpha = 0.1) +
  geom_abline(slope = 1,intercept = 0,color = "red",linetype = "dashed") +
  theme_bw() +
  theme() +
  scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  scale_x_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  #geom_vline(xintercept = median(df_summary_tot2$pure_exonic),col="yellow",linetype = "dashed") +
  #geom_hline(yintercept = median(df_summary_tot2$pure_intronic),col="yellow",linetype = "dashed") +
  coord_fixed() +
  ggtitle("pure intronic vs pure exonic")

wrap_plots(list(p1,p2,p3))

# add side ditribtution for the scatter
ggMarginal(p3, type="histogram")
# ggMarginal(p3, type="density")
# ggMarginal(p3, type="boxplot")
# ggMarginal(p3, type="densigram")

# -------------------------------------------------------------------------
# test proportion per gene
df_summary_gene <- PropIntronicPerGene(sparse_matrix_exon = ldat_fix$spliced,
                                       sparse_matrix_intron = ldat_fix$unspliced)

df_summary_gene %>%
  ggplot(aes(x = prop_intron)) + geom_histogram() + theme_bw()

df_summary_gene %>%
  arrange(desc(prop_intron),desc(count_intron))

df_summary_gene %>%
  arrange(prop_intron,desc(count_exon))

df_summary_gene %>%
  summarise(tot_exon = sum(count_exon),
            tot_intron = sum(count_intron)) %>%
  mutate(prop_exon = tot_exon/(tot_exon + tot_intron),
         prop_intron = tot_intron/(tot_exon + tot_intron))

# plot the different prop per count
p1 <- df_summary_gene %>%
  filter(prop_intron == 1) %>%
  ggplot(aes(x=count_intron)) +
  geom_histogram() +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = c(0,1,5,10,20,50,100,200,500,1000,2000,5000)) +
  ggtitle("hist number of genes with prop_intron == 1")

p2 <- df_summary_gene %>%
  filter(prop_intron == 0) %>%
  ggplot(aes(x=count_exon)) +
  geom_histogram() +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = c(1,10,100,10^3,10^4,10^5,10^6,10^7)) +
  ggtitle("hist number of genes with prop_intron == 0")

p1 + p2

# see the proportion of intronic reads per feature
df_test <- rtracklayer::import(gtf_file_loc) %>%
  as.data.frame() %>%
  dplyr::filter(type == "gene")

df_summary_gene %>%
  left_join(df_test,by = c("gene" = "gene_name")) %>%
  ggplot(aes(y=width,x=prop_intron)) + geom_point(alpha = 0.1,shape = 1)+theme_bw()+
  scale_y_continuous(trans = "log10") +
  geom_smooth(method = "lm",se = F) +
  # facet_wrap(~dataset,nrow = 1) +
  theme(strip.background = element_blank())

df_summary_gene %>%
  left_join(df_test,by = c("gene" = "gene_name")) %>%
  ggplot(aes(x=gene_type,y=width)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2),alpha = 0.01,shape = 1)+theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45)) +
  scale_y_continuous(trans = "log10") +
  # facet_wrap(~dataset,nrow = 1) +
  theme(strip.background = element_blank())

df_summary_gene %>%
  left_join(df_test,by = c("gene" = "gene_name")) %>%
  ggplot(aes(x=gene_type,y=prop_intron)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2),alpha = 0.01,shape = 1)+theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45)) +
  # scale_y_log10() +
  # facet_wrap(~dataset,nrow = 1) +
  theme(strip.background = element_blank())

# -------------------------------------------------------------------------
# test proportion per cell
df_summary_cell <- PropIntronicPerCell(sparse_matrix_exon = ldat_filter$sparse_matrix_exon,
                                       sparse_matrix_intron = ldat_filter$sparse_matrix_intron)

head(df_summary_cell)

# calculate more coveriated on the full object
# crete the object
datasc <- CreateSeuratObject(counts = mtx_raw, project = "test_intron", min.cells = 20, min.features = 200)

# add more covariates
datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^MT-")
datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^HB[^(P)]")

# add the PropIntronicPerCell
datasc <- AddMetaData(object = datasc,metadata = df_summary_cell %>% column_to_rownames("barcode"))

# buld the plotting funcction
panelfun <- function(data, mapping) { 
  ggplot(data = data, mapping = mapping)+  
    # geom_smooth(method = "lm")+ 
    geom_point(alpha=0.2,size = 0.2,shape = 1) 
} 

# check the correlation bewteen different covariates in the metadata
datasc@meta.data %>%
  ggpairs(columns = 2:9,
          lower = list(continuous = wrap(panelfun))) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1,angle=45),
        axis.text.y = element_text(hjust = 1,angle=45)) +
  scale_y_continuous(trans = "log1p") +
  scale_x_continuous(trans = "log1p")

# focus only on some panels of correlations
p1 <- df_summary_gene %>%
  ggplot(aes(x = prop_intron)) + geom_histogram() + theme_bw() +
  ggtitle("hist number of genes per prop_intron")

p2 <- df_summary_gene %>%
  filter(prop_intron == 1) %>%
  ggplot(aes(x=count_intron)) +
  geom_histogram() +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = c(0,1,5,10,20,50,100,200,500,1000,2000,5000)) +
  ggtitle("hist number of genes with prop_intron == 1")

p3 <- df_summary_gene %>%
  filter(prop_intron == 0) %>%
  ggplot(aes(x=count_exon)) +
  geom_histogram() +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = c(1,10,100,10^3,10^4,10^5,10^6,10^7)) +
  ggtitle("hist number of genes with prop_intron == 0")

p4 <- datasc@meta.data %>%
  ggplot(aes(x=percent.mt,y=prop_intron)) +
  geom_point(alpha=0.5) + theme_bw() +
  ggtitle("correlation prop_intron perc.mt")

p5 <- datasc@meta.data %>%
  ggplot(aes(x=percent.ribo,y=prop_intron)) +
  geom_point(alpha=0.5) + theme_bw() +
  ggtitle("correlation prop_intron perc.ribo")


(p1|p4|p5)

# simulate the creation of the objects from pure intron pure exon  --------
# here I simulate the sample preprocessing loading the matrices.
# I will skip the QC step for this sample test

# # load the intron exon pure matrices
# ldat <- Read10XIntronExon(gtf_obj = GTF_obj,file = "../../data/test_connect5k_PBMC/results/PureIntronExon/merged/connect_5k_pbmc_NGSC3_ch1_gex_1/")
# mtx_raw <- Read10X("../../data/test_connect5k_PBMC/results/cellranger/merged/connect_5k_pbmc_NGSC3_ch1_gex_1/outs/filtered_feature_bc_matrix/")                             
# 
# # define the filtered barcodes to make a common comparison
# cell_barcodes <- colnames(mtx_raw)
# 
# # filter the pure matrices over the same set of barcodes
# ldat_filter <- lapply(ldat, function(x){
#   x[,cell_barcodes]
# })

# compile a list with all the matrices of expression
list_mtx_all_filtered <- list(all = mtx_raw,
                              pure_exon = ldat_filter$sparse_matrix_exon,
                              pure_intron = ldat_filter$sparse_matrix_intron)

# confirm the dimensitons
lapply(list_mtx_all_filtered, function(x){
  dim(x)
})

# drop the min.cells and min.features but use the whithlist of barcodes from the regular cellranger output
# use the same sets of barcodes fro all the matrices to generate the objects
# run the standard processing in a loop
list_sobj <- pmap(list(list_mtx_all_filtered,names(list_mtx_all_filtered)),function(mtx,nm){
  # track the progress
  print(nm)
  
  # generate the seurat object
  sobj_total <- CreateSeuratObject(counts = mtx,
                                   project = nm,
                                   # meta.data = data.combined.all@meta.data,
                                   min.cells = 0, min.features = 0) %>%
    # this is needed as the cell cycle scoring is done on the data slot, which would be empty
    Seurat::NormalizeData(verbose = T)
  
  # add QC metadata
  DefaultAssay(sobj_total) <- "RNA"
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  sobj_total <- CellCycleScoring(sobj_total, s.features = s.genes, g2m.features = g2m.genes)
  sobj_total$percent.mt <- PercentageFeatureSet(sobj_total, pattern = "^MT-")
  sobj_total$percent.ribo <- PercentageFeatureSet(sobj_total, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  # add also the percentage of globin. in this dataset it is not meaningful as there is no blood
  sobj_total$percent.globin <- Seurat::PercentageFeatureSet(sobj_total,pattern = "^HB[^(P)]")
  
  # rescale the data for regressing out the sources of variation do not scale all the genes.
  # if needed for some plots, I can scale them before the heatmap call. for speeding up the computation I will scale only the HVF
  sobj_total <- sobj_total %>%
    # skip the normalizatio that has been already performed at the beginning
    # Seurat::NormalizeData(verbose = T) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>%
    # I can scale the missing features afterwards now focus on the highly variable one for speed purposes
    ScaleData(vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"), verbose = T) %>% 
    # run this if you want to scale all the variables
    # ScaleData(vars.to.regress = c("percent.mt.harmony","nCount_RNA.harmony","S.Score.harmony","G2M.Score.harmony"), verbose = T,features = all.genes) %>% 
    RunPCA(npcs = 30, verbose = T) %>% 
    RunUMAP(reduction = "pca", dims = 1:30,return.model = TRUE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, by = 0.1)) %>%
    identity()
  
  return(sobj_total)
})

# save the output
saveRDS(list_sobj,"../../out/object/100_list_sobj_connect_5k_pbmc_NGSC3_ch1_gex_1.rds")

# sample plot
list_plot <- pmap(list(list_sobj,names(list_sobj)),function(x,nm){
  FeaturePlot(x,features = "nFeature_RNA",order=T) + ggtitle(nm) + scale_color_viridis_c(option = "turbo")
})

wrap_plots(list_plot)
