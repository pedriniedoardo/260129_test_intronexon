# AIM ---------------------------------------------------------------------
# try to read in the output of cellranger run W or WO intronic reads counting

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# read in the data --------------------------------------------------------
# read in the filtered matrix from the filterd bam
data_Wintron <- Read10X(data.dir = "/beegfs/scratch/ric.absinta/ric.absinta/analysis/BS_run05/test/Sample_W8_untreated_Wintron/outs/filtered_feature_bc_matrix")

# read in the reference dataset, the one derived from the original bam files
data_WOintron <- Read10X("/beegfs/scratch/ric.absinta/ric.absinta/analysis/BS_run05/test/Sample_W8_untreated_WOintron/outs/filtered_feature_bc_matrix")

# wrangling ---------------------------------------------------------------
# compare the two matrices
dim(data_Wintron)
dim(data_WOintron)

# make the comparison across the common barcodes
common_barcodes <- intersect(colnames(data_Wintron),colnames(data_WOintron))

# compare the summaries
df_Wintron <- rowSums(data_Wintron[,common_barcodes]) %>%
  data.frame() %>%
  setNames("count_Wintron") %>%
  rownames_to_column("gene")

df_WOintron <- rowSums(data_WOintron[,common_barcodes]) %>%
  data.frame() %>%
  setNames("count_WOintron") %>%
  rownames_to_column("gene")

df_summary_tot <- left_join(df_Wintron,df_WOintron,by = "gene") %>%
  mutate(delta = count_Wintron - count_WOintron)

df_summary_tot %>%
  filter(delta != 0) %>%
  dim()

# use the pseudolog transfromation to show both positive and negative values
df_summary_tot %>%
  ggplot(aes(x=delta)) +
  geom_histogram() +
  scale_x_continuous(trans = "pseudo_log",breaks = c(-10000,-1000,-100,-10,-1,0,10,100,1000,10000,100000)) +
  # scale_x_continuous(trans = "log1p",breaks = c(0,1,10,100,1000,10000,100000,10^6)) +
  scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000)) +
  theme_bw()

df_summary_tot %>%
  ggplot(aes(x=count_WOintron,y=count_Wintron)) +
  geom_point(alpha = 0.1) +
  # add the value where delta !=0
  geom_abline(slope = 1,intercept = 0,color = "red",linetype = "dashed") + theme_bw()+theme()+
  scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
  scale_x_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7))
