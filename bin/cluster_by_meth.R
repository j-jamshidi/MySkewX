#!/usr/bin/env Rscript

library(NanoMethViz)
library(tidyverse)
library(doParallel)
library(foreach)

# Patched version of NanoMethViz:::cluster_reads that adds values_fn = mean to
# pivot_wider, preventing list-col errors when the BAM contains both 5mC and
# 5hmC modification calls for the same (read_name, pos).
cluster_reads_fixed <- function(x, chr, start, end, min_pts = 5) {
  assertthat::assert_that(
    methods::is(x, "ModBamResult"),
    assertthat::is.string(chr) || (is.factor(chr) && assertthat::is.scalar(chr)),
    assertthat::is.number(start) && assertthat::is.number(end),
    assertthat::is.number(min_pts) && min_pts >= 1
  )
  methy_data <- NanoMethViz::query_methy(x, chr, start, end)
  if (nrow(methy_data) == 0) {
    stop(glue::glue("no reads containing methylation data found in specified region"))
  }
  methy_data <- methy_data %>% dplyr::filter(.data$pos >= start & .data$pos < end)
  read_stats <- NanoMethViz:::get_read_stats(methy_data)
  max_span <- max(read_stats$span)
  keep_reads <- read_stats$read_name[read_stats$span > 0.9 * max_span]
  methy_data <- methy_data %>% dplyr::filter(.data$read_name %in% keep_reads)
  mod_mat <- methy_data %>%
    dplyr::select("read_name", "pos", "mod_prob") %>%
    dplyr::arrange(.data$pos) %>%
    tidyr::pivot_wider(names_from = "pos", values_from = "mod_prob", values_fn = mean) %>%
    NanoMethViz:::df_to_matrix()
  if (nrow(mod_mat) < min_pts) {
    stop(glue::glue("fewer reads available ({nrow(mod_mat)} reads) than minimum cluster size 'min_pts' ({min_pts})"))
  }
  mod_mat_filled <- mod_mat[order(rownames(mod_mat)), ]
  col_missingness <- NanoMethViz:::mat_col_map(mod_mat_filled, NanoMethViz:::missingness)
  mod_mat_filled <- mod_mat_filled[, col_missingness < 0.6]
  row_missingness <- NanoMethViz:::mat_row_map(mod_mat_filled, NanoMethViz:::missingness)
  mod_mat_filled <- mod_mat_filled[row_missingness < 0.3, ]
  for (i in seq_len(nrow(mod_mat_filled))) {
    mod_mat_filled[i, is.na(mod_mat_filled[i, ])] <- mean(mod_mat_filled[i, ], na.rm = TRUE)
  }
  if (nrow(mod_mat_filled) < min_pts) {
    stop(glue::glue("fewer reads available ({nrow(mod_mat_filled)} reads) than minimum cluster size 'min_pts' ({min_pts})"))
  }
  dbsc <- dbscan::hdbscan(mod_mat_filled, minPts = min_pts)
  clust_df <- data.frame(read_name = rownames(mod_mat_filled), cluster_id = dbsc$cluster)
  clust_df %>%
    dplyr::inner_join(read_stats, by = "read_name") %>%
    dplyr::arrange(.data$cluster_id) %>%
    dplyr::mutate(
      cluster_id = as.factor(.data$cluster_id),
      start = as.integer(.data$start),
      end = as.integer(.data$end),
      span = as.integer(.data$span)
    )
}

apply_cluster_reads_parallel <- function(mbr, bed, min_pts, num_cores = 4) {
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  process_row <- function(row) {
    current_row <- bed[row, ]

    row_cluster <- tryCatch({
      cluster_reads_fixed(mbr, current_row$chr, current_row$start, current_row$end, min_pts = min_pts)},
      error = function(err){
        message(paste("Error in row", row, ":", err$message))
        return(NA)
      })

    if (all(is.na(row_cluster))) {
      return(NULL)
    }

    row_cluster <- row_cluster %>% mutate(CGI_id = paste0(current_row$chr, ":", current_row$start, "-", current_row$end), chr = current_row$chr, start = current_row$start, end = current_row$end)

    row_cluster <- row_cluster %>% group_by(cluster_id) %>% mutate(avg_cluster_methylation = mean(mean))

    if (nlevels(row_cluster$cluster_id) == 2) {
      low_mC <- row_cluster %>% filter(cluster_id %in% c("1", "2")) %>% pull(avg_cluster_methylation) %>% min()
      high_mC <- row_cluster %>% filter(cluster_id %in% c("1", "2")) %>% pull(avg_cluster_methylation) %>% max()
      row_cluster <- row_cluster %>% mutate(assigned_X = case_when(avg_cluster_methylation == low_mC ~ "Xa", avg_cluster_methylation == high_mC ~ "Xi", TRUE ~ "NA"))
    } else {
      row_cluster$assigned_X <- NA
    }

    return(row_cluster)
  }

  res <- foreach(row = 1:nrow(bed), .combine = bind_rows,
                 .packages = c("tidyverse", "NanoMethViz"),
                 .export = c("cluster_reads_fixed")) %dopar% {
    process_row(row)
  }

  stopCluster(cl)

  res <- bind_rows(res)

  return(res)
}

calculate_skew_by_block <- function(clustered_reads, haplotyped_reads){
  if (is.null(clustered_reads) || nrow(clustered_reads) == 0 || !("assigned_X" %in% colnames(clustered_reads))) {
    stop("No reads were successfully clustered. Check that the BAM index is up to date (run: samtools index <bam>).")
  }
  clustered_reads <- clustered_reads %>% filter(assigned_X %in% c("Xa","Xi"))
  clustered_reads <- clustered_reads %>% distinct(read_name, .keep_all = TRUE)
  df2 <- left_join(clustered_reads,haplotyped_reads)
  df2 <- df2 %>% filter(!is.na(HP))

  counts_by_block <- df2 %>% group_by(PS, assigned_X, HP) %>% summarise(counts = n(), .groups = "drop")

  skew_by_block <- counts_by_block %>%
    unite(combi, assigned_X, HP) %>%
    mutate(combi = recode(combi, "Xa_1" = "H1_Xa", "Xa_2" = "H2_Xa", "Xi_1" = "H1_Xi", "Xi_2" = "H2_Xi")) %>%
    pivot_wider(id_cols = PS, names_from = combi, values_from = counts, values_fill = 0)

  # Ensure all 4 expected columns exist; pivot_wider omits columns for combinations
  # that never appear in the data (values_fill only fills within existing columns)
  for (col in c("H1_Xa", "H2_Xa", "H1_Xi", "H2_Xi")) {
    if (!col %in% names(skew_by_block)) skew_by_block[[col]] <- 0L
  }

  skew_by_block <- skew_by_block %>%
    mutate(H1_Xa_skew = (H1_Xa + H2_Xi) / (H1_Xa + H1_Xi + H2_Xa + H2_Xi))

  return(skew_by_block)
}


# Retrieve the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

lib <- args[1]
bam <- args[2]
BED <- read_tsv(args[3], col_names = c("chr", "start", "end"))
haplotyped_reads <- read_tsv(args[4], col_names = c("read_name", "HP", "PS"))
ncpus <- strtoi(args[5])

mbr <- ModBamResult(
    methy = ModBamFiles(
        samples = lib,
        paths = bam
    ),
    samples = data.frame(
        sample = lib,
        group = 1
    )
)

clustered_reads <- apply_cluster_reads_parallel(mbr, BED, min_pts = 5, num_cores = ncpus)
write_tsv(clustered_reads, paste0(lib, "_CGIX_clustered_reads.tsv.gz"))

skew <- calculate_skew_by_block(clustered_reads, haplotyped_reads)

write_tsv(skew, paste0(lib,"_CGIX_skew.tsv.gz"))
