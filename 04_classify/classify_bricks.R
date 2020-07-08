#!/usr/bin/Rscript
#-------------------------------------------------------------------------------
# CLASSIFY SENTINEL-2 BRICKS
#-------------------------------------------------------------------------------
suppressMessages(library(caret))
suppressMessages(library(dplyr))
suppressMessages(library(raster))
suppressMessages(library(sits))

#---- Parameters ----

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
    stop("This script takes parameters:
         A brick type [approx|raw],
         a brick directory,
         a samples file (RDS of a sits tibble with time series),
         a comma-separated list of labels,
         a comma-separated list of bands,
         a version number,
         a base directory for storing the results.",  call. = FALSE)
}

b_type         <-                                args[[1]]
brick_dir      <-                                args[[2]]
samples_file   <-                                args[[3]]
used_labels    <- sort(unlist(stringr::str_split(args[[4]], ',')))
used_bands     <- sort(unlist(stringr::str_split(args[[5]], ',')))
version_number <-                                args[[6]]
out_base_dir   <-                                args[[7]]

stopifnot(b_type %in% c("approx", "raw"))
stopifnot(dir.exists(brick_dir))
stopifnot(file.exists(samples_file))
#stopifnot(dir.exists(out_base_dir))

#---- Setup ----

tmp_directory <- "./tmp"
dir.create(file.path(tmp_directory, "masked"))
raster::rasterOptions(tmpdir = tmp_directory)
raster::tmpDir()

out_file_template <- samples_file %>%
    basename() %>%
    tools::file_path_sans_ext()

#---- Util ----

source("./other/util.R")

# Helper for doing the classification.
classify <- function(used_bands, used_labels, brick_dir, samples_file,
                     sits_method, out_dir, version_number){
    samples_tb <- samples_file %>%
        readRDS() %>%
        select_bands(used_bands) %>%
        dplyr::filter(label %in% used_labels) %>%
        ensurer::ensure_that(nrow(.) > 0, err_desc = "Samples missing") %>%
        ensurer::ensure_that(length(unique(.$label)) == length(used_labels),
                             err_desc = "The samples are missing labels!") %>%
        ensurer::ensure_that(length(sits::sits_bands(.)) == length(used_bands),
                             err_desc = "The samples are missing bands!")
    brick_tb <- brick_dir %>%
        get_brick_md() %>%
        dplyr::filter(brick_type == b_type, resolution == "10m") %>%
        dplyr::filter(band %in% used_bands) %>%
        dplyr::arrange(tile, img_date, band) %>%
        ensurer::ensure_that(!"" %in% .$band) %>%
        ensurer::ensure_that(length(unique(.$tile)) == 1,
                             err_desc = sprintf("More than one tile found: %s",
                                                brick_dir)) %>%
        ensurer::ensure_that(length(unique(.$img_date)) == 1,
                             err_desc = sprintf("More than one date found: %s",
                                                brick_dir)) %>%
        ensurer::ensure_that(nrow(.) == length(used_bands),
                             err_desc = sprintf("Bands not found: %s",
                                                used_bands))
    stopifnot(all(match(brick_tb$band, colnames(samples_tb$time_series[[1]])[-1]) == 1:length(brick_tb$band)))
    cube <- sits::sits_cube(service = "BRICK",
                            name = "sentinel-bricks",
                            satellite = "SENTINEL2",
                            sensor = "MSI",
                            timeline = seq(unique(brick_tb$img_date),
                                           by = 10, length.out = 36),
                            bands = brick_tb$band,
                            files = brick_tb$file_path)
    write(sits::sits_bands(cube), file = file.path(out_dir, "sits_bands.txt"))
    write(used_labels, file = file.path(out_dir, "sits_labels.txt"))
    model <- samples_tb %>%
        sits::sits_train(ml_method = sits_method) %>%
        (function(x){
             saveRDS(x, file = file.path(out_dir, "model.rds"))
             return(x)
        })
    probability_map <- sits::sits_classify(data = cube, ml_model = model,
                                           multicores = 16, memsize = 4,
                                           output_dir = out_dir,
                                           version = version_number)
    classification_map <- sits::sits_label_classification(probability_map,
                                                          smoothing = "bayesian",
                                                          output_dir = out_dir,
                                                          version = version_number)
    invisible(list(probability_map, classification_map))
}

#---- Classify using Random Forest ----

out_dir <- file.path(out_base_dir,
                     b_type,
                     out_file_template,
                     paste(used_bands, collapse = '-'),
                     paste(used_labels, collapse = '-'),
                     "random-forest_1000/")
if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)
print(sprintf("Saving results to: %s", out_dir))

classify(used_bands,
         used_labels,
         brick_dir,
         samples_file,
         sits_method = sits::sits_rfor(num_trees = 1000) ,
         out_dir,
         version_number)
