# AIM ---------------------------------------------------------------------
# read in sample velocyto outputs

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(velocyto.R)

# read in data ------------------------------------------------------------
# http://pklab.med.harvard.edu/velocyto/mouseBM/SCG71.loom
ldat <- ReadVelocity(file = "../../data/test/SCG71.loom")

# processing --------------------------------------------------------------
# This returns a list of matrices: spliced, unspliced, ambiguous
spliced <- ldat$spliced
unspliced <- ldat$unspliced
ambiguose <- ldat$ambiguous

# If you want to create a Seurat object solely from this:
seu_velo <- as.Seurat(x = ldat)
