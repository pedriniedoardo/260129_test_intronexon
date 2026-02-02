# AIM ---------------------------------------------------------------------
# given the summarized features extracted from each dataset calculate the intronic metric.
# filter out the features that have more than 50% of intronic reads

# libraries ---------------------------------------------------------------
library(tidyverse)
library(UpSetR)

# read in the data --------------------------------------------------------
list_df_tot <- readRDS("../out/object/list_intronic_metrics.rds")

# wrangling ---------------------------------------------------------------
# pull the genes that have at least 50% of intronic reads
list_feature_high_intron <- lapply(list_df_tot,function(x){
  x %>%
    filter(prop_intron > 0.5) %>%
    pull(symbol) %>%
    unique()
})

# try the upset plot version
pdf("../out/plot/02_upset_intronic.pdf",width = 6,height = 4,onefile=FALSE)
upset(fromList(list_feature_high_intron), order.by = "freq") 
dev.off()

# save the intersections
df1 <- lapply(list_feature_high_intron,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "path")
head(df1)

df2 <- data.frame(gene = unique(unlist(list_feature_high_intron)))
head(df2)

# now loop through each individual gene and pick the list of all the intersections they belong to
df_int <- lapply(df2$gene,function(x){
  # pull the name of the intersections
  intersection <- df1 %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows() %>%
  arrange(int)
head(df_int,n=20)

# confirm the data and the list are congruent
df_int %>% 
  group_by(int) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))

# save the table of intersections
df_int %>%
  write_tsv("../out/table/02_df_intersection.tsv")
