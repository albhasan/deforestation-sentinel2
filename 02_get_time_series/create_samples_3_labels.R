#!/usr/bin/Rscript

# Recode the labels in the samples files from 5 to 3.

suppressMessages(library(dplyr))

base_dir <- "./data/validation"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
      stop("This script takes 1 parameter: A path to a RDS file with a SITS tibble, including its time series.",
           call. = FALSE)
}
rds_file <- args[[1]]

rds_file %>%
    ensurer::ensure_that(file.exists(.), err_desc = "File not found: ") %>%
    readRDS() %>%
    dplyr::filter(label != "Water") %>%
    dplyr::mutate(label = dplyr::recode(label,
                                        #"Forest"
                                        #"Deforestatio"
                                        #NonForest
                                        Pasture      = "NonForest",
                                        NatNonForest = "NonForest")) %>%
    ensurer::ensure_that(nrow(.) > 0, err_desc = "No rows!") %>%
    ensurer::ensure_that(length(unique(.$label)) == 3,
                         err_desc = "Label missmatch!") %>%
    saveRDS(file.path(base_dir, paste0(tools::file_path_sans_ext(basename(rds_file)),
                                       "_3l.rds")))
