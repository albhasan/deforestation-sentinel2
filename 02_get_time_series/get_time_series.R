#!/usr/bin/env Rscript
suppressMessages(library(dplyr))
suppressMessages(library(sf))
suppressMessages(library(sits))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
      stop("This script takes 3 parameters: An input file (CSV), a type of brick, and an output file (RDS).",
           call. = FALSE)
}

source("./other/util.R")

sample_csv <- args[[1]]
b_type     <- args[[2]]
out_file   <- args[[3]]

brick_path <- "./brick_sentinel2"

brick_tb <- brick_path %>%
    get_brick_md() %>%
    dplyr::filter(brick_type == b_type,
                  resolution == "10m", file_ext == "vrt") %>%
    dplyr::arrange(tile, img_date, band) %>%
    ensurer::ensure_that(nrow(.) > 0, err_desc = "No files found!") %>%
    ensurer::ensure_that(length(unique(.$img_date)) == 1,
                         err_desc = sprintf("More than one date found: %s",
                                            brick_path)) %>%
    ensurer::ensure_that(length(unique(.$tile)) == 1,
                         err_desc = sprintf("More than one tile found: %s",
                                            brick_path)) %>%
    ensurer::ensure_that(!"" %in% .$band) %>%
    ensurer::ensure_that(all(table(.$file_band) == 1),
                         err_desc = "Repeated bands found!") %>%
    ensurer::ensure_that(!any(is.na(.$band)), err_desc = "Unknown bands!") %>%
    ensurer::ensure_that(!any(is.na(.$brick_type)), err_desc = "Unknown brick!")

# Create a coverage for each brick.
cube <- sits::sits_cube(service = "BRICK", name = "sentinel-bricks",
                        satellite = "SENTINEL2", sensor = "MSI",
                        timeline = seq(unique(brick_tb$img_date),
                                       by = 10, length.out = 36) ,
                        bands = brick_tb$band,
                        files = brick_tb$file_path)

# Interpolate NAs in the time series of a sits tibble.
fill_in_ts <- function(x){
    n <- nrow(x)
    .data <- x[2:ncol(x)]
    .data[.data < -1] <- NA
    .data[.data > 1] <- NA
    res <- lapply(.data, function(x1, n){
                      if(sum(!is.na(x1)) < 2)
                          return(rep(NA, times = n))
                      return(approx(x1, n = n)$y)
                            }, n = n)
    return(dplyr::bind_cols(x[,1], res))
}

# clean the time series
sits::sits_get_data(cube, file = sample_csv) %>%
    dplyr::mutate(time_series = purrr::map(time_series, fill_in_ts)) %>%
    saveRDS(file = out_file)
