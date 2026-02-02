# AIM ---------------------------------------------------------------------
# try to read in the output of cellranger run on the filtered bam files translated back to fastqs


# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# read in the data --------------------------------------------------------
# read in the filtered matrix from the filterd bam
data_filt <- Read10X(data.dir = "/beegfs/scratch/ric.absinta/ric.absinta/analysis/BS_run05/test/Sample_W8_Wintron_all/outs/filtered_feature_bc_matrix/")

# read in the reference dataset, the one derived from the original bam files
data_ref <- Read10X("/beegfs/scratch/ric.absinta/ric.absinta/analysis/BS_run05/test/Sample_W8_untreated_Wintron/outs/filtered_feature_bc_matrix/")

# wrangling ---------------------------------------------------------------
# compare the two matrices
dim(data_filt)
dim(data_ref)

# the two matrices are already different in size.

# make the comparison across the common barcodes
common_barcodes <- intersect(colnames(data_filt),colnames(data_ref))

# compare the summaries
df_filt <- rowSums(data_filt[,common_barcodes]) %>%
  data.frame() %>%
  setNames("count_filt") %>%
  rownames_to_column("gene")

df_ref <- rowSums(data_ref[,common_barcodes]) %>%
  data.frame() %>%
  setNames("count_ref") %>%
  rownames_to_column("gene")

df_summary_tot <- left_join(df_filt,df_ref,by = "gene") %>%
  mutate(delta = count_ref - count_filt)

df_summary_tot %>%
  filter(delta != 0) %>%
  dim()

df_summary_tot %>%
  ggplot(aes(x=delta)) +
  geom_histogram() +
  scale_x_continuous(trans = "pseudo_log",breaks = c(-10000,-1000,-100,-10,-1,0,10,100,1000,10000,100000)) +
  scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000)) +
  theme_bw()

# df_summary_tot1 <- df_summary_tot %>%
#   filter(delta == 0)
# 
# df_summary_tot2 <- df_summary_tot %>%
#   filter(delta != 0)
# 
# df_summary_tot1 %>%
#   ggplot(aes(x=count_ref,y=count_filt)) +
#   geom_point(alpha = 0.01,col="gray") +
#   # add the value where delta !=0
#   geom_point(data = df_summary_tot2,aes(col=delta!=0),alpha=0.2) +
#   geom_abline(slope = 1,intercept = 0,color = "red",linetype = "dashed") + theme_bw()+theme()+
#   scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7)) +
#   scale_x_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,100000,10^6,10^7))
