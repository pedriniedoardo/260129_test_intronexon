# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(skimr)
library(schard)
library(ComplexHeatmap)
library(viridis)
library(patchwork)

# custom ------------------------------------------------------------------
# define the jaccard score function
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  # might not be robust in case of duplicates entries
  # union <- length(a) + length(b) - intersection 
  union <- length(union(a, b))
  return (intersection/union)
}
# test
a <- c('potato', 'tomotto', 'chips', 'baloon')
b <- c('car', 'chips', 'bird', 'salt')

jaccard(a, b)

# read in data ------------------------------------------------------------
# read in the dataset after runnign cellbender
test_wcellb <- readRDS("/idle/ric.cosr/ric.cosr/maurizio.aurora/bop/scrnaseq_cosr_standard_workflow_cr10/results/Seurat/object/CyteTypeR_annotation_cb.Rds")
DimPlot(test_wcellb,label = T)

DimPlot(test_wcellb,group.by = "cytetype_annotation_RNA_snn_res.0.3",label = T)
FeaturePlot(test_wcellb,features = "nCount_RNA",order = T) + scale_color_viridis_c(option = "turbo", trans = "log1p",breaks = scales::breaks_log(n = 5))


# -------------------------------------------------------------------------
# plot the covariated
cov_test_01 <- c("nCount_RNA", "nFeature_RNA")
# x <- "percent.ribo"
list_plot_technical_01 <- lapply(cov_test_01,function(x){
  
  plot <- FeaturePlot(test_wcellb,features = x,order = T,
                      reduction = "umap",
                      raster = T) +
  scale_color_viridis_c(option = "turbo", trans = "log",
                        breaks = scales::breaks_log(5),
                        # breaks = c(0,10,30,60,100,1000,3000,6000,10000,150000), 
                        labels = scales::label_number(scale_cut = scales::cut_short_scale())
                        ) +
    ggtitle(x)
  return(plot)
})

cov_test_02 <- c("percent.mt", "percent.ribo")
list_plot_technical_02 <- lapply(cov_test_02,function(x){
  
  plot <- FeaturePlot(test_wcellb,features = x,order = T,
                      reduction = "umap",
                      raster = T) +
    scale_color_viridis_c(option = "turbo", trans = "sqrt",
                          breaks = scales::breaks_pretty(5)) +
    ggtitle(x)
  return(plot)
})

wrap_plots(c(list_plot_technical_01,list_plot_technical_02))
ggsave("../../out/plot/000_UMAPCluster_technical.pdf",width = 10,height = 8)

# make violin plots
list_plot_technical <- lapply(c(cov_test_01,cov_test_02), function(x){ 
  test <- VlnPlot(object = test_wcellb,features = x, group.by = "cytetype_cellOntologyTerm_RNA_snn_res.0.3",raster = T)
  return(test)
})

# make it a dataframe
# x <- list_plot[[1]]
df_violin_technical <- lapply(list_plot_technical,function(x){ 
  df <- x[[1]]$data 
  
  # extract the name of the gene 
  feature <- colnames(df)[1] 
  
  df %>% 
    mutate(feature = feature) %>% 
    setNames(c("value","ident","feature")) 
}) %>% 
  bind_rows()

df_plot_violin_technical_summary <- df_violin_technical %>%
  group_by(feature) %>%
  summarise(med_score = median(value))

# plot at maximum 500 cells per group
set.seed(123)
df_plot_points <- df_violin_technical %>% 
  group_by(ident, feature) %>%
  # Pro-Tip: use slice_sample() instead of sample_n(). 
  # If a cluster has fewer than 500 cells, sample_n() will crash, 
  # but slice_sample() will safely just take all available cells.
  slice_sample(n = 100,replace = F) %>% 
  ungroup()

df_violin_technical %>%
  ggplot(aes(x = ident, y = value)) + 
  geom_violin(scale = "width")+ 
  #geom_boxplot(outlier.shape = NA,position = position_dodge(width=0.9),width=0.05) + 
  geom_point(data = df_plot_points,position = position_jitter(width = 0.2),alpha = 0.1,size = 0.5) + 
  facet_wrap(~feature,ncol = 1,scales = "free_y") + 
  theme_bw() + 
  geom_hline(data = df_plot_violin_technical_summary,aes(yintercept = med_score),linetype="dashed",col="red") +
  scale_y_continuous(trans = "sqrt",
                     breaks = scales::breaks_pretty(5),
                     labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1,angle = 45),
        )
ggsave("../../out/plot/000_ViolinCluster_technical.pdf",width = 9,height = 10)

# -------------------------------------------------------------------------

# read in the dataset not running cellbender
test_wocellb <- readRDS("/idle/ric.cosr/ric.cosr/maurizio.aurora/bop2/scrnaseq_cosr_standard_workflow_cr10/results/Seurat/object/CyteTypeR_annotation.Rds")
DimPlot(test_wocellb,label = T)

# read in the original file
test_ref <- schard::h5ad2seurat("/idle/ric.cosr/ric.cosr/maurizio.aurora/bop/scrnaseq_cosr_standard_workflow_cr10/results/Seurat/object/0b259a0f-cb45-44f9-a97f-b5dbd8478e10.h5ad")
DimPlot(test_ref,group.by = "cell_type",label = T)

# wrangling ---------------------------------------------------------------
# explore the dataset
# cytetype_ontologyID_RNA_snn_res.0.3
test_wocellb@meta.data %>%
  select(cytetype_annotation_RNA_snn_res.0.3,cytetype_cellOntologyTerm_RNA_snn_res.0.3,cytetype_annotation_RNA_snn_res.0.3,cytetype_cellState_RNA_snn_res.0.3,cytetype_RNA_snn_res.0.3,cytetype_ontologyTerm_RNA_snn_res.0.3) %>%
  skim()

test_wcellb@meta.data %>%
  group_by(cytetype_annotation_RNA_snn_res.0.3,
           cytetype_cellOntologyTerm_RNA_snn_res.0.3,
           cytetype_cellState_RNA_snn_res.0.3,
           cytetype_RNA_snn_res.0.3,
           cytetype_ontologyTerm_RNA_snn_res.0.3) %>%
  summarise()

# -------------------------------------------------------------------------
# compare the levels of the annotations
meta_test_wocellb <- test_wocellb@meta.data %>%
  group_by(cytetype_annotation_RNA_snn_res.0.3,
           cytetype_cellOntologyTerm_RNA_snn_res.0.3,
           cytetype_cellState_RNA_snn_res.0.3,
           cytetype_RNA_snn_res.0.3,
           cytetype_ontologyTerm_RNA_snn_res.0.3) %>%
  summarise()

meta_test_wocellb %>%
  rownames_to_column("cluster") %>%
  pivot_longer(names_to = "test",values_to = "ID",cytetype_annotation_RNA_snn_res.0.3:cytetype_ontologyTerm_RNA_snn_res.0.3) %>%
  filter(str_detect(ID,pattern = "macro"))

meta_test_wcellb <- test_wcellb@meta.data %>%
  group_by(cytetype_annotation_RNA_snn_res.0.3,
           cytetype_cellOntologyTerm_RNA_snn_res.0.3,
           cytetype_cellState_RNA_snn_res.0.3,
           cytetype_RNA_snn_res.0.3,
           cytetype_ontologyTerm_RNA_snn_res.0.3) %>%
  summarise()

meta_test_wcellb %>%
  rownames_to_column("cluster") %>%
  pivot_longer(names_to = "test",values_to = "ID",cytetype_annotation_RNA_snn_res.0.3:cytetype_ontologyTerm_RNA_snn_res.0.3) %>%
  filter(str_detect(ID,pattern = "macro"))

# -------------------------------------------------------------------------
# compare the meta
meta_ref <- test_ref@meta.data %>%
  select(sample_id,
         cell_type) %>%
  rownames_to_column("barcodes") %>%
  separate(barcodes,into = c("barcode_id","id","sample"),sep = "-") %>%
  mutate(barcode_id2 = paste(barcode_id,id,sep = "-"))

meta_wocellb <- test_wocellb@meta.data %>%
  select(orig.ident,
         cytetype_annotation_RNA_snn_res.0.3,
         cytetype_cellOntologyTerm_RNA_snn_res.0.3,
         cytetype_cellState_RNA_snn_res.0.3,
         cytetype_RNA_snn_res.0.3,
         cytetype_ontologyTerm_RNA_snn_res.0.3) %>%
  rownames_to_column("barcodes") %>%
  separate(barcodes,into = c("sample","barcode_id2"),sep = "_")

meta_wcellb <- test_wcellb@meta.data %>%
  select(orig.ident,
         cytetype_annotation_RNA_snn_res.0.3,
         cytetype_cellOntologyTerm_RNA_snn_res.0.3,
         cytetype_cellState_RNA_snn_res.0.3,
         cytetype_RNA_snn_res.0.3,
         cytetype_ontologyTerm_RNA_snn_res.0.3) %>%
  rownames_to_column("barcodes") %>%
  separate(barcodes,into = c("sample","barcode_id2"),sep = "_")

# infer the sample_id
# meta_wocellb %>%
#   group_by(sample) %>%
#   summarise()

# meta_ref %>%
#   left_join(meta_wocellb,by = "barcode_id2") %>%
#   filter(!is.na(orig.ident)) %>%
#   group_by(sample_id,orig.ident) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n)) %>%
#   print(n=20)

# HuLN5_sLMQ   SRR33207205
# HuLN5_sCTRL  SRR33207206
# HuLN4_sLMQ   SRR33207207
# HuLN4_sCTRL  SRR33207208
# HuLN3_sLMQ   SRR33207209
# HuLN3_sCTRL  SRR33207211

LUT_sample <- data.frame(sample_id_our = 
             c("SRR33207205",
               "SRR33207206",
               "SRR33207207",
               "SRR33207208",
               "SRR33207209",
               "SRR33207211"),
           sample_id_ref = 
             c("HuLN5_sLMQ",
               "HuLN5_sCTRL",
               "HuLN4_sLMQ",
               "HuLN4_sCTRL",
               "HuLN3_sLMQ",
               "HuLN3_sCTRL"))

# add semple conversion as metadata
meta_wocellb_full <- meta_wocellb %>%
  left_join(LUT_sample,by = c("orig.ident" = "sample_id_our")) %>%
  mutate(barcode_pivot = paste(barcode_id2,sample_id_ref,sep = "_")) %>%
  select(barcode_pivot,sample_id_our = orig.ident,sample_id_ref,cytetype_annotation_RNA_snn_res.0.3:cytetype_ontologyTerm_RNA_snn_res.0.3)

meta_wcellb_full <- meta_wcellb %>%
  left_join(LUT_sample,by = c("orig.ident" = "sample_id_our")) %>%
  mutate(barcode_pivot = paste(barcode_id2,sample_id_ref,sep = "_")) %>%
  select(barcode_pivot,sample_id_our = orig.ident,sample_id_ref,cytetype_annotation_RNA_snn_res.0.3:cytetype_ontologyTerm_RNA_snn_res.0.3)

meta_ref_full <- meta_ref %>%
  left_join(LUT_sample,by = c("sample_id" = "sample_id_ref")) %>%
  mutate(barcode_pivot = paste(barcode_id2,sample_id,sep = "_")) %>%
  select(barcode_pivot,sample_id_our,sample_id_ref = sample_id,cell_type)

# test the label transfer from the two fixed metadata
test_compare <- left_join(meta_wocellb_full,meta_ref_full,by = c("barcode_pivot","sample_id_our","sample_id_ref")) %>%
  # add the label missing cell in case we have the barcode but they do not have it
  mutate(cell_type_fix = case_when(is.na(cell_type)~"missing",
                                   T~cell_type))

test_compare_cb <- left_join(meta_wcellb_full,meta_ref_full,by = c("barcode_pivot","sample_id_our","sample_id_ref")) %>%
  # add the label missing cell in case we have the barcode but they do not have it
  mutate(cell_type_fix = case_when(is.na(cell_type)~"missing",
                                   T~cell_type))

# check the proportios
test_compare %>%
  # group_by(cytetype_annotation_RNA_snn_res.0.3,cell_type) %>%
  group_by(cell_type) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

test_compare %>%
  # group_by(cytetype_annotation_RNA_snn_res.0.3,cell_type) %>%
  group_by(sample_id_our) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

test_compare_cb %>%
  # group_by(cytetype_annotation_RNA_snn_res.0.3,cell_type) %>%
  group_by(sample_id_our) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

test_compare_cb %>%
  # group_by(cytetype_annotation_RNA_snn_res.0.3,cell_type) %>%
  group_by(cell_type) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

test_compare_cb %>%
  # group_by(cytetype_annotation_RNA_snn_res.0.3,cell_type) %>%
  group_by(cell_type,sample_id_our) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  filter(is.na(cell_type))

test_compare %>%
  # group_by(cytetype_annotation_RNA_snn_res.0.3,cell_type) %>%
  group_by(cell_type_fix) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# test the label transfer from the two fixed metadata
test_compare2 <- meta_ref_full %>%
  filter(!is.na(sample_id_our)) %>%
  left_join(meta_wocellb_full,,by = c("barcode_pivot","sample_id_our","sample_id_ref"))

# check the proportios
test_compare2 %>%
  # group_by(cytetype_annotation_RNA_snn_res.0.3,cell_type) %>%
  group_by(cytetype_annotation_RNA_snn_res.0.3) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# -------------------------------------------------------------------------

# build the dataset for the correlatino plot cross the two annotation per cells
# loop across annotations

# vec_our_annotation <- c("cytetype_annotation_RNA_snn_res.0.3",
#                         "cytetype_cellOntologyTerm_RNA_snn_res.0.3",
#                         "cytetype_cellState_RNA_snn_res.0.3",
#                         "cytetype_RNA_snn_res.0.3",
#                         "cytetype_ontologyTerm_RNA_snn_res.0.3")

vec_our_annotation <- c("cytetype_cellOntologyTerm_RNA_snn_res.0.3")

# anno <- "cytetype_annotation_RNA_snn_res.0.3"
list_hm <- lapply(vec_our_annotation, function(anno){
  # track the progress
  print(anno)
  
  # build the crossing reference
  df_crossing <- crossing(id_ref = unique(test_compare$cell_type_fix),
                          id_query = unique(test_compare[[anno]]))
  
  # build the scatter plot
  # id_ref <- "missing"
  # id_query <- "Activated mast cell"
  df_jaccard_score <- pmap(list(id_ref = df_crossing$id_ref,
                                id_query = df_crossing$id_query), function(id_ref,id_query){
                                  
                                  # calculate the jaccard score
                                  a <- test_compare %>%
                                    rownames_to_column("barcode") %>%
                                    filter(cell_type_fix == id_ref) %>% pull(barcode_pivot)
                                  
                                  b <- test_compare %>%
                                    rownames_to_column("barcode") %>%
                                    filter(.data[[anno]] == id_query) %>% pull(barcode_pivot)
                                  
                                  jaccard_score <- jaccard(a,b)
                                  
                                  # build a data.frame
                                  df <- data.frame("id_ref" = id_ref,
                                                   "id_query" = id_query,
                                                   "jaccard_score" = jaccard_score)
                                  return(df)
                                }) %>%
    bind_rows()
  
  # check the table
  head(df_jaccard_score)
  
  # shape it as a matrix
  mat_jaccard_score <- df_jaccard_score %>%
    pivot_wider(names_from = id_ref,values_from = jaccard_score) %>%
    column_to_rownames("id_query")
  
  df_jaccard_score %>%
    filter(id_ref == "missing")
  
  test_compare %>%
    filter(cell_type_fix == "missing")
  
  # plot the matrix
  ht_02 <- Heatmap(mat_jaccard_score,
                   column_title = anno,
                   name = "Jaccard\nscore",
                   # col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
                   col = viridis::viridis(option = "turbo",n = 20),
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 8),
                   column_names_side = "bottom",
                   column_names_gp = gpar(fontsize = 8),
                   column_names_rot = 45,
                   row_dend_reorder = FALSE,
                   column_dend_reorder = FALSE,
                   row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                   column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                   show_column_names = T,
                   show_row_names = T)
  
  # pdf("../out/plot/129_heatmap_jaccard_res0.9.pdf",height = 4,width = 5)
  # hm <- draw(ht_02,heatmap_legend_side = "left",padding = unit(c(40, 2, 2, 40), "mm"))
  # dev.off()
  # return(hm)
  return(ht_02)
})

list_hm2 <- lapply(list_hm, function(x){
  # 2. Capture them as "grob" objects (graphical objects)
  hm <- grid.grabExpr(draw(x,heatmap_legend_side = "left",padding = unit(c(20, 2, 2, 30), "mm")))
  return(hm)
})
wrap_plots(list_hm2)
ggsave("../../out/plot/000_test_battle_jaccard_wocb.pdf",width = 10,height = 8)

# test on the cellbender output

# anno <- "cytetype_annotation_RNA_snn_res.0.3"
list_hm_cb <- lapply(vec_our_annotation, function(anno){
  # track the progress
  print(anno)
  
  # build the crossing reference
  df_crossing <- crossing(id_ref = unique(test_compare_cb$cell_type_fix),
                          id_query = unique(test_compare_cb[[anno]]))
  
  # build the scatter plot
  # id_ref <- "missing"
  # id_query <- "Adaptive-like NK cell"
  df_jaccard_score <- pmap(list(id_ref = df_crossing$id_ref,
                                id_query = df_crossing$id_query), function(id_ref,id_query){
                                  
                                  # calculate the jaccard score
                                  a <- test_compare_cb %>%
                                    rownames_to_column("barcode") %>%
                                    filter(cell_type_fix == id_ref) %>% pull(barcode_pivot)
                                  
                                  b <- test_compare_cb %>%
                                    rownames_to_column("barcode") %>%
                                    filter(.data[[anno]] == id_query) %>% pull(barcode_pivot)
                                  
                                  jaccard_score <- jaccard(a,b)
                                  
                                  # build a data.frame
                                  df <- data.frame("id_ref" = id_ref,
                                                   "id_query" = id_query,
                                                   "jaccard_score" = jaccard_score)
                                  return(df)
                                }) %>%
    bind_rows()
  
  # check the table
  head(df_jaccard_score)
  
  # shape it as a matrix
  mat_jaccard_score <- df_jaccard_score %>%
    pivot_wider(names_from = id_ref,values_from = jaccard_score) %>%
    column_to_rownames("id_query")
  
  df_jaccard_score %>%
    filter(id_ref == "missing")
  
  test_compare %>%
    filter(cell_type_fix == "missing")
  
  # plot the matrix
  ht_02 <- Heatmap(mat_jaccard_score,
                   name = "Jaccard\nscore",
                   column_title = anno,
                   # col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
                   col = viridis::viridis(option = "turbo",n = 20),
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 8),
                   column_names_side = "bottom",
                   column_names_gp = gpar(fontsize = 8),
                   column_names_rot = 45,
                   row_dend_reorder = FALSE,
                   column_dend_reorder = FALSE,
                   row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                   column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                   show_column_names = T,
                   show_row_names = T)
  
  # pdf("../out/plot/129_heatmap_jaccard_res0.9.pdf",height = 4,width = 5)
  # hm <- draw(ht_02,heatmap_legend_side = "left",padding = unit(c(40, 2, 2, 40), "mm"))
  # dev.off()
  # return(hm)
  return(ht_02)
})

list_hm2_cb <- lapply(list_hm_cb, function(x){
  # 2. Capture them as "grob" objects (graphical objects)
  hm <- grid.grabExpr(draw(x,heatmap_legend_side = "left",padding = unit(c(20, 2, 2, 30), "mm")))
  return(hm)
})
wrap_plots(list_hm2_cb)
ggsave("../../out/plot/000_test_battle_jaccard_wcb.pdf",width = 10,height = 8)
