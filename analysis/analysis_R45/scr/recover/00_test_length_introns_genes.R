# AIM ---------------------------------------------------------------------
# test how to extract the intron lenght from the gtf file

# libraries ---------------------------------------------------------------
library(GenomicRanges)
library(rtracklayer)

# Load GTF file
gtf_file <- "/beegfs/scratch/ric.cosr/pedrini.edoardo/reference/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
gtf_data <- import(gtf_file)

# # Filter for exons and introns
# exons <- gtf_data[gtf_data$type == "exon"]
# introns <- gtf_data[gtf_data$type == "intron"]

# to understant the difference between the differend entries of the file focus on one gene only
test <- gtf_data %>%
  as.data.frame() %>%
  filter(gene_name == "ACTB")

test %>%
  group_by(gene_name,type) %>%
  summarise(sum_width = sum(width))
  
# check the math
4320+12+30+6089


# -------------------------------------------------------------------------
# Filter to get only exons
exon_test_data <- test %>%
  filter(type == "exon") %>%
  select(seqnames, start, end, strand, gene_id, transcript_id)

# Function to calculate intron lengths
calculate_intron_length <- function(exons) {
  # Sort exons by start position
  exons <- exons[order(exons$start), ]
  
  # Calculate introns between exons
  intron_starts <- exons$end[-nrow(exons)] + 1
  intron_ends <- exons$start[-1] - 1
  
  # Calculate the length of each intron
  intron_lengths <- intron_ends - intron_starts + 1
  
  # Sum all intron lengths for the transcript
  total_intron_length <- sum(intron_lengths[intron_lengths > 0])
  return(total_intron_length)
}

# Calculate total intron length for each transcript
intron_lengths
exon_test_data %>%
  group_by(gene_id, transcript_id) %>%
  summarise(n = n())
  # summarize(total_intron_length = calculate_intron_length(cur_data()))

  
test_exon <- exon_test_data %>%
  mutate(id = paste0(gene_id,"-",transcript_id)) %>%
  split(f = .$id)
  
exons <- test_exon$`ENSG00000075624-ENST00000414620`

# Sort exons by start position
exons <- exons[order(exons$start), ]
    
# Calculate introns between exons
intron_starts <- exons$end[-nrow(exons)] + 1
intron_ends <- exons$start[-1] - 1
    
# Calculate the length of each intron
intron_lengths <- intron_ends - intron_starts + 1
    
# Sum all intron lengths for the transcript
total_intron_length <- sum(intron_lengths[intron_lengths > 0])
    return(total_intron_length)
  })

  
  summarize(total_intron_length = calculate_intron_length(cur_data()))
# View the results
print(intron_lengths)

