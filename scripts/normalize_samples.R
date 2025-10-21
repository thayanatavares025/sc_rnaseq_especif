#!/bin/R

args <- commandArgs(trailingOnly = TRUE)
srr_dir <- args[1]
seurat_object <- args[2]

# thresholds
n_features <- 2000
scale_factor <- 1000

# methods
normalization_method <- "LogNormalize"
features_selection_method <- "vst"

# n_features - genes number with more variation
normalize <- function(seurat_object, srr_dir) {
  seurat_object <- NormalizeData(
    object = seurat_object, 
    normalization.method = normalization_method,
    scale.factor = scale_factor
  )
  
  seurat_object <- FindVariableFeatures(
    seurat_object, selection.method = features_selection_method, nfeatures = n_features
  )
  
  saveRDS(seurat_object, paste0(srr_dir, '_norm_seurat_object.RDS'))
}

normalize(seurat_object=seurat_object, srr_dir=srr_dir)