# AIM ---------------------------------------------------------------------
# read in the summarized feature counts per common barcode

# libraries ---------------------------------------------------------------
library(tidyverse)
library(rtracklayer)
library(GGally)
library(patchwork)
library(Seurat)

# read in the data --------------------------------------------------------
# read in the gene annotation. Notice that all the features.tsv files contains all the annotaiton needed to convert symbols to ensemble gene id and vice-versa
# LUT_gene <- read_tsv("/beegfs/scratch/ric.absinta/ric.absinta/analysis/BS_run05/test/Sample_W8_untreated_Wintron/outs/raw_feature_bc_matrix/features.tsv.gz",col_names = F) %>%
#   dplyr::rename(enembl_gene_id = X1,symbol = X2, type = X3)
# 
# # pull more annotations from the gtf file of 
# df_gtf <- rtracklayer::import("/beegfs/scratch/ric.cosr/pedrini.edoardo/reference/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf") %>%
#   as.data.frame() %>%
#   # filter only the gene annotation
#   dplyr::filter(type %in% c("gene"))
# 
# # check what genes are missing
# LUT_gene_full <- LUT_gene %>%
#   left_join(df_gtf,by = c("enembl_gene_id" = "gene_id"))
# # dplyr::filter(is.na(biotype))
# 
# # save the table for future reference
# write_tsv(LUT_gene_full,"../data/LUT_gene_full_Human.tsv")

LUT_gene_full <- read_tsv("../data/LUT_gene_full_Human.tsv")

# define the folder for the sample of interest
list_folder <- list(K16 = "/beegfs/scratch/ric.cosr/ric.cosr/CRtest_genova_pedrini/cellranger/cellranger7/output/K16/outs/",
                    S1 = "/beegfs/scratch/ric.cosr/ric.cosr/CRtest_genova_pedrini/cellranger/cellranger7/output/S1/outs/",
                    Unsorted = "/beegfs/scratch/ric.cosr/ric.cosr/CRtest_genova_pedrini/cellranger/cellranger7/output/Unsorted/outs/",
                    BSedo = "/beegfs/scratch/ric.absinta/ric.absinta/analysis/BS_run05/test/Sample_W8_untreated_Wintron/outs/"
                    )

# read in the summarized counts
# x <- list_folder[[2]]
list_feature <- lapply(list_folder, function(x){
  feature_exon <- read_table(paste0(x,"feature_counts_possorted_genome_bam_exon_filterd_common.txt"),col_names = F) %>%
    dplyr::rename(count_exon = X1,gene = X2)
  
  feature_intron <- read_table(paste0(x,"feature_counts_possorted_genome_bam_intron_filterd_common.txt"),col_names = F) %>%
    dplyr::rename(count_intron = X1,gene = X2)
  
  # this is the feature count directly from the unfiltered bam file. this do not coincide with the sum (feature-wise) of the previous two tables because it retains also the non filtered reads
  feature_all <- read_table(paste0(x,"feature_counts_possorted_genome_bam_original_common.txt"),col_names = F) %>%
    dplyr::rename(count_all = X1,gene = X2)
  
  # read in the feature counts derived from the filtering process on the original bam file
  feature_all_filtered <- read_table(paste0(x,"feature_counts_possorted_genome_bam_all_filtered_common.txt"),col_names = F) %>%
    dplyr::rename(count_allFilter = X1,gene = X2)
  
  return(list(feature_exon = feature_exon,
              feature_intron = feature_intron,
              feature_all = feature_all,
              feature_all_filtered = feature_all_filtered))
})

saveRDS(list_feature,"../out/object/list_feature.rds")

# EDA ---------------------------------------------------------------------
# notice that in the summarised expression dataset, the zero expressiong gene are not reported
# also the number of genes per summary is different

# # x <- list_feature$BSedo
# list(feature_all = x$feature_all,
#      feature_all_filtered = x$feature_all_filtered,
#      feature_exon = x$feature_exon,
#      feature_intron = x$feature_intron) %>%
#   lapply(function(x){
#     x %>%
#       setNames(c("count","gene")) %>%
#       summarise(n_gene = n(),n_gene_zero = mean(count<1))
#   })

# generate the intronic proportion metric ---------------------------------
# list_feature <- readRDS("../out/object/list_feature.rds")
list_df_tot <- lapply(list_feature,function(x){
  
  # check the counts form the filtered and intronic + exonic
  df_test02 <- full_join(x$feature_exon,x$feature_intron,by = "gene") %>%
    full_join(x$feature_all_filtered,by = "gene") %>%
    # convert he missing counts to 0
    pivot_longer(names_to = "table",values_to = "count",-gene) %>%
    mutate(count = case_when(is.na(count)~0,
                             T~count)) %>% 
    # # confirm there are no more NAs
    # filter(is.na(count))
    pivot_wider(names_from = table,values_from = count) %>%
    # verify that the sum of intron + exon is equivalent to the allFilter entry
    mutate(count_exon_intron = count_exon + count_intron) %>%
    # add a pseudocount to account also for the zero expressing genes durign log scaling of the axis
    mutate(count_intron_adj = count_intron + 1,
           count_exon_adj = count_exon + 1,
           count_allFilter_adj = count_allFilter + 1)
  
  # define the metric as percentage of intronic reads per gene
  df_tot <- df_test02 %>%
    mutate(prop_intron = count_intron / count_allFilter) %>%
    arrange(desc(prop_intron)) %>%
    mutate(gene = str_remove_all(gene,"GX:Z:")) %>%
    left_join(LUT_gene_full,by = c("gene" = "enembl_gene_id"))
  
  return(df_tot)
})

# save teh list in an object
saveRDS(list_df_tot,"../out/object/list_intronic_metrics.rds")

# run some tests ----------------------------------------------------------
# list_feature <- readRDS("../out/object/list_feature.rds")
# test 01 confirm the table of counts before and after filtering the bam file is not the same
# to avoid skipping some genes that have zero expression in one table, do not use inner_join, but full_join.
# convert he missing counts to 0
df_01 <- lapply(list_feature,function(x){
  df_test01 <- full_join(x$feature_all,x$feature_all_filtered,by = "gene") %>%
    # convert he missing counts to 0
    pivot_longer(names_to = "table",values_to = "count",-gene) %>%
    mutate(count = case_when(is.na(count)~0,
                             T~count)) %>% 
    # # confirm there are no more NAs
    # filter(is.na(count))
    pivot_wider(names_from = table,values_from = count) %>%
    # add a pseudocount to account also for the zero expressing genes durign log scaling of the axis
    mutate(count_all_adj = count_all + 1,
           count_allFilter_adj = count_allFilter + 1)
  
  return(df_test01)
}) %>%
  bind_rows(.id = "dataset")
  
df_01 %>%
  ggplot(aes(x=count_all,y=count_allFilter)) +
  geom_point(shape=1,alpha=0.1)+
  scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  scale_x_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  theme_bw() +
  facet_wrap(~dataset) +
  theme(strip.background = element_blank()) +
  geom_abline(col="red",lty="dashed")
ggsave("../out/plot/01_scatter_allFiltered_vs_all.pdf",width = 6,height = 6)


df_01 %>%
  filter(dataset == "BSedo") %>%
  ggplot(aes(x=count_all,y=count_allFilter)) +
  geom_point(shape=1,alpha=0.1)+
  scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  scale_x_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  theme_bw() +
  facet_wrap(~dataset) +
  theme(strip.background = element_blank()) +
  geom_abline(col="red",lty="dashed")

# is there any count in all that is bigger than the one in allFIlter
df_01 %>%
  mutate(diff = count_all - count_allFilter) %>%
  filter(diff < 0)

# check the counts form the filtered and intronic + exonic
df_02 <- lapply(list_feature,function(x){
  df_test02 <- full_join(x$feature_exon,x$feature_intron,by = "gene") %>%
    full_join(x$feature_all_filtered,by = "gene") %>%
    # convert he missing counts to 0
    pivot_longer(names_to = "table",values_to = "count",-gene) %>%
    mutate(count = case_when(is.na(count)~0,
                             T~count)) %>% 
    # # confirm there are no more NAs
    # filter(is.na(count))
    pivot_wider(names_from = table,values_from = count) %>%
    # verify that the sum of intron + exon is equivalent to the allFilter entry
    mutate(count_exon_intron = count_exon + count_intron) %>%
    # add a pseudocount to account also for the zero expressing genes durign log scaling of the axis
    mutate(count_intron_adj = count_intron + 1,
           count_exon_adj = count_exon + 1,
           count_allFilter_adj = count_allFilter + 1)
  
  return(df_test02)
}) %>%
  bind_rows(.id = "dataset")

#
df_02 %>%
  mutate(delta = count_allFilter - count_exon_intron) %>%
  mutate(test = !(delta == 0)) %>%
  summarise(non_zero = sum(test))

# confirm the above with a plot
df_02 %>%
  ggplot(aes(x = count_allFilter, y = count_exon_intron)) +
  geom_point(shape = 1, alpha = 0.1) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  facet_wrap(~dataset) +
  theme(strip.background = element_blank()) +
  geom_abline(slope = 1,intercept = 0,col="red",lty = "dashed")
ggsave("../out/plot/01_scatter_allFiltered_vs_ExonAndIntron.pdf",width = 6,height = 6)

df_02 %>%
  filter(dataset == "BSedo") %>%
  ggplot(aes(x = count_allFilter, y = count_exon_intron)) +
  geom_point(shape = 1, alpha = 0.1) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  facet_wrap(~dataset) +
  theme(strip.background = element_blank()) +
  geom_abline(slope = 1,intercept = 0,col="red",lty = "dashed")

# plot the count_allFilter vs the exon only
df_02 %>%
  ggplot(aes(x = count_exon, y = count_allFilter)) +
  geom_point(shape = 1, alpha = 0.1) +
  scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  scale_x_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  theme_bw() +
  facet_wrap(~dataset) +
  theme(strip.background = element_blank()) +
  geom_abline(slope = 1,intercept = 0,col="red",lty = "dashed")
ggsave("../out/plot/01_scatter_allFiltered_vs_Exon.pdf",width = 6,height = 6)


df_02 %>%
  filter(dataset == "BSedo") %>%
  ggplot(aes(x = count_exon, y = count_allFilter)) +
  geom_point(shape = 1, alpha = 0.1) +
  scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  scale_x_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  theme_bw() +
  facet_wrap(~dataset) +
  theme(strip.background = element_blank()) +
  geom_abline(slope = 1,intercept = 0,col="red",lty = "dashed")

# confirm the that total counts count_allFilter are always greater than the only exonic counts
df_02 %>%
  mutate(diff = count_allFilter - count_exon) %>%
  filter(diff < 0)


# list_df_tot <- readRDS("../out/object/list_intronic_metrics.rds")
df_tot <- list_df_tot %>%
  bind_rows(.id = "dataset")

# general distribution of the intronic reads proportion
P0 <- df_tot %>%
  ggplot(aes(x=prop_intron)) + geom_histogram() + theme_bw() +facet_wrap(~dataset,nrow = 1) +
  theme(strip.background = element_blank())

# strendness seems to be uncorrelated with intron proportions
P1 <- df_tot %>%
  ggplot(aes(x=strand,y=prop_intron)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2),alpha = 0.1,shape = 1)+theme_bw() +facet_wrap(~dataset,nrow = 1) +
  theme(strip.background = element_blank())

# there seems top be a correlation with width
P2 <- df_tot %>%
  ggplot(aes(y=width,x=prop_intron)) + geom_point(alpha = 0.1,shape = 1)+theme_bw()+scale_y_log10()+geom_smooth(method = "lm",se = F) + facet_wrap(~dataset,nrow = 1) +
  theme(strip.background = element_blank())

# check the trend with gene type
P3 <- df_tot %>%
  ggplot(aes(x=gene_type,y=prop_intron)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2),alpha = 0.01,shape = 1)+theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45)) + facet_wrap(~dataset,nrow = 1) +
  theme(strip.background = element_blank())

# check if width and gene type are correlated
P4 <- df_tot %>%
  ggplot(aes(x=gene_type,y=width)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2),alpha = 0.01,shape = 1)+theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45)) + scale_y_log10() +
  facet_wrap(~dataset,nrow = 1) +
  theme(strip.background = element_blank())

# assamble the panel
P0/P1/P2/P3/P4
ggsave("../out/plot/01_EDA_intronProp.pdf",width = 16,height = 20)
