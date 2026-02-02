# AIM ---------------------------------------------------------------------
# read in the features using the testSort and regular pipeline
# the aim is to see if the output is the same

# generate the counts from the test file and compare it with the original file
feature_intron_test <- read_table("/beegfs/scratch/ric.cosr/ric.cosr/CRtest_genova_pedrini/cellranger/cellranger7/output/Unsorted/outs/feature_counts_possorted_genome_bam_intron_filterd_common_test.txt",col_names = F) %>%
  dplyr::rename(count_intron_test = X1,gene = X2)

feature_intron <- read_table("/beegfs/scratch/ric.cosr/ric.cosr/CRtest_genova_pedrini/cellranger/cellranger7/output/Unsorted/outs/feature_counts_possorted_genome_bam_intron_filterd_common.txt",col_names = F) %>%
  dplyr::rename(count_intron = X1,gene = X2)

dim(feature_intron)
dim(feature_intron_test)

# join the two datasets
# the two datasets are the same
inner_join(feature_intron,feature_intron_test,by = "gene") %>%
  mutate(diff = count_intron - count_intron_test) %>%
  filter(diff != 0)
