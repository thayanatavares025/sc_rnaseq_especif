#!/bin/R

args <- commandArgs(trailingOnly = TRUE)

matrix_dir  <- args[1]
csv_metrics <- args[2]
sra_id      <- args[3]
sample_name <- sra_id

# thresholds
thr_estimate_n_cells      <- 300
thr_mean_reads_per_cells  <- 25000
thr_median_genes_per_cell <- 900
thr_median_umi_per_cell   <- 1000
thr_n_feature_rna_min     <- 300
thr_n_feature_rna_max     <- 7500
thr_percent_mito          <- 25
thr_n_observed_cells      <- 300

load_libraries <- function() {
  libs <- c("readr", "dplyr", "ggplot2", "Seurat")
  
  for (lib in libs) {
    if (!requireNamespace(lib, quietly = TRUE)) {
      message(paste("Installing missing package:", lib))
      install.packages(lib, repos = "https://cloud.r-project.org")
    }
    library(lib, character.only = TRUE)
  }
}

string_to_numeric <- function(x) {
  as.numeric(gsub('%', '', x))
}

check_sample_quality <- function(sample_metrics) {
  check_list <- c(FALSE, FALSE, FALSE, FALSE)
  if (sample_metrics$`Estimated Number of Cells` >= thr_estimate_n_cells)       check_list[1] <- TRUE
  if (sample_metrics$`Mean Reads per Cell`      >= thr_mean_reads_per_cells)    check_list[2] <- TRUE
  if (sample_metrics$`Median Genes per Cell`    >= thr_median_genes_per_cell)   check_list[3] <- TRUE
  if (sample_metrics$`Median UMI Counts per Cell` >= thr_median_umi_per_cell)   check_list[4] <- TRUE
  
  if (all(check_list)) {
    status_flag <- "SAMPLE:SUCCESS"
  } else {
    if (sample_metrics$`Sequencing Saturation` <= 70) {
      status_flag <- "SAMPLE:FIXABLE"
    } else {
      status_flag <- "SAMPLE:FAILURE"
    }
  }
  return(status_flag)
}

create_seurat_object <- function(matrices_dir, sample_name) {
  expression_matrix <- Read10X(data.dir = matrices_dir)
  seurat_object <- CreateSeuratObject(counts = expression_matrix, project = sample_name)
  return(seurat_object)
}

calculate_mito_content <- function(seurat_object) {
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  return(seurat_object)
}

add_meta_data <- function(seurat_object, sample_name, meta_data_file) {
  meta_data_df <- read_csv(file = meta_data_file, show_col_types = FALSE) %>%
    filter(Run == sample_name)
  
  if (nrow(meta_data_df) >= 1) {
    meta_data_df <- meta_data_df %>% select(-matches("^...$"), everything())
    seurat_object@meta.data <- cbind(seurat_object@meta.data, meta_data_df[rep(1, ncol(seurat_object)), , drop=FALSE])
  } else {
    cat("No metadata found for sample:", sample_name, "\n")
  }
  return(seurat_object)
}

get_info_seurat <- function(seurat_object, min_mito = thr_percent_mito) {
  n_cells  <- ncol(seurat_object)
  n_genes  <- nrow(seurat_object)
  mean_genes      <- mean(seurat_object$nFeature_RNA)
  mean_mito       <- mean(seurat_object$percent.mt)
  high_mito_count <- sum(seurat_object$percent.mt > min_mito)
  return(list(
    n_cells = n_cells,
    n_genes = n_genes,
    mean_genes = mean_genes,
    mean_mito = mean_mito,
    high_mito_count = high_mito_count
  ))
}

filter_cells <- function(sample_name, seurat_object, status_flag) {
  # TODO [FILTER_LOGIC]: verify whether filtering should be skipped
  # when status_flag is already "SAMPLE:FAILURE".
  n_total_cells <- ncol(seurat_object)
  
  seurat_object <- subset(
    x = seurat_object,
    subset = (nFeature_RNA >= thr_n_feature_rna_min & nFeature_RNA < thr_n_feature_rna_max)
  )
  
  n_thr_features_cells <- ncol(seurat_object)
  
  seurat_object <- subset(
    x = seurat_object,
    subset = percent.mt < thr_percent_mito
  )
  
  n_observed_cells <- ncol(seurat_object)
  
  if (status_flag == "SAMPLE:SUCCESS") {
    status_flag <- ifelse(n_observed_cells >= thr_n_observed_cells, "CELL:SUCCESS", "CELL:FAILURE")
  }

  readr::write_lines(
    status_flag,
    file = paste0(sample_name, "_status.txt")
  )
  
  return(list(
    seurat_object = seurat_object,
    status_flag = status_flag,
    n_total_cells = n_total_cells,
    n_thr_features_cells = n_thr_features_cells,
    n_observed_cells = n_observed_cells
  ))
}

process_and_save_metrics <- function(
    sample_metrics,
    sample_name, n_total_cells,
    n_thr_features_cells, n_observed_cells, status_flag
) {
  sample_metrics_upgrade <- sample_metrics %>%
    select(
      `Estimated Number of Cells`, `Mean Reads per Cell`,
      `Median Genes per Cell`, `Median UMI Counts per Cell`,
      `Sequencing Saturation`
    ) %>%
    rename(
      estimate_n_cells     = `Estimated Number of Cells`,
      mean_reads_per_cell  = `Mean Reads per Cell`,
      median_genes_per_cell= `Median Genes per Cell`,
      median_umi_per_cell  = `Median UMI Counts per Cell`,
      seq_saturation       = `Sequencing Saturation`
    ) %>%
    mutate(
      sample_id             = sample_name,
      n_total_cells         = n_total_cells,
      n_thr_features_cells  = n_thr_features_cells,
      n_observed_cells      = n_observed_cells,
      status_flag           = status_flag
    ) %>%
    select(sample_id, status_flag, estimate_n_cells, mean_reads_per_cell, median_genes_per_cell,
           median_umi_per_cell, seq_saturation, n_total_cells, n_thr_features_cells, n_observed_cells)
  
  readr::write_csv(sample_metrics_upgrade, file = paste0(sample_name, "_metrics_upgrade.csv"))
}

save_seurat_object <- function(seurat_object, sample_name, status_flag) {
  if (status_flag == "CELL:SUCCESS") {
    dir.create("objects", showWarnings = FALSE)
    saveRDS(seurat_object, file = paste0("objects/", sample_name, "_seurat_object.RDS"))
    message("Seurat object saved successfully for sample: ", sample_name)
  } else {
    message("Seurat object NOT saved for sample: ", sample_name, 
            " (status_flag = ", status_flag, ")")
  }
}

save_qc_summary <- function(sample_name, info_before_filters, info_after_filters){
  qc_summary_sample <- data.frame(
    Sample = sample_name,
    Metric = c("Number of cells", 
               "Number of genes", 
               "Mean genes per cell", 
               "Mean % mitochondrial", 
               "Cells with >10% mitochondrial"),
    Before_Filter = c(info_before_filters$n_cells, 
                      info_before_filters$n_genes, 
                      round(info_before_filters$mean_genes, 1), 
                      round(info_before_filters$mean_mito, 2), 
                      info_before_filters$high_mito_count),
    After_Filter = c(info_after_filters$n_cells, 
                     info_after_filters$n_genes, 
                     round(info_after_filters$mean_genes, 1), 
                     round(info_after_filters$mean_mito, 2), 
                     info_after_filters$high_mito_count)
  )
  write.csv(qc_summary_sample, file = paste0(sample_name, "_qc_summary.csv"), row.names = FALSE)
}



load_libraries()

sample_metrics <- readr::read_csv(file = csv_metrics, show_col_types = FALSE) %>%
  mutate_if(is.character, string_to_numeric)

status_flag <- check_sample_quality(sample_metrics = sample_metrics)

seurat_object <- create_seurat_object(matrices_dir = matrix_dir, sample_name = sample_name)
seurat_object <- calculate_mito_content(seurat_object)

# TODO [METADATA]: allow optional metadata import as args[4] and enable the line below
# seurat_object <- add_meta_data(seurat_object, sample_name, meta_file)

info_before_filters <- get_info_seurat(seurat_object)

list_result <- filter_cells(
  sample_name = sample_name,
  seurat_object = seurat_object,
  status_flag = status_flag
)

seurat_object <- list_result$seurat_object
status_flag <- list_result$status_flag

process_and_save_metrics(
  sample_metrics = sample_metrics,
  sample_name = sample_name,
  n_total_cells = list_result$n_total_cells,
  n_thr_features_cells = list_result$n_thr_features_cells,
  n_observed_cells = list_result$n_observed_cells,
  status_flag = list_result$status_flag
)

info_after_filters <- get_info_seurat(seurat_object)

save_qc_summary(
  sample_name = sample_name,
  info_before_filters = info_before_filters,
  info_after_filters = info_after_filters
)

save_seurat_object(seurat_object, sample_name, status_flag)
