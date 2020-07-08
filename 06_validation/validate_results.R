#!/usr/bin/Rscript

suppressMessages(library(dplyr))
suppressMessages(library(raster))

set.seed(666)

source("./other/util.R")

# TODO: Transform into script parameter.
result_dir <- "./results9"

my_recode <- c(Deforestatio = "Deforestation", Forest = "Forest",
               "NonForest" = "Non-Forest")
my_recode2 <- my_recode %>%
    magrittr::set_names(as.character(1:length(.)))

# Get the validation points
samples_tb <- tibble::tribble(~used_bands, ~extra_samples_file,
                "blue-bnir-green-nnir-red-swir1-swir2", "./data/validation/validation_bands.rds",
                "evi-ndmi-ndvi",                        "./data/validation/validation_indices.rds") %>%
    dplyr::mutate(samples = purrr::map(extra_samples_file, readRDS))

samples_tb$samples <- lapply(samples_tb$samples, function(x){
                                 x <- x %>%
                                     dplyr::mutate(label = dplyr::recode(label, !!!my_recode)) %>%
                                     return()
                })

# NOTE: save to pangea
samples_tb$samples[[1]] %>%
    dplyr::mutate(id = dplyr::row_number()) %>%
    dplyr::select(id, longitude, latitude, start_date, end_date, label) %>%
    readr::write_csv("./pangea/validation_dataset_bands.csv")
samples_tb$samples[[2]] %>%
    dplyr::mutate(id = dplyr::row_number()) %>%
    dplyr::select(id, longitude, latitude, start_date, end_date, label) %>%
    readr::write_csv("./pangea/validation_dataset_indices.csv")

results_tb <- result_dir %>%
    get_results() %>%
    dplyr::filter(img_type == "classification", class_type == "full") %>%
    ensurer::ensure_that(nrow(.) == 2,
                         err_desc = "Unexpected number of classification rasters.") %>%
    dplyr::mutate(label_file = file.path(dirname(file_path), "sits_labels.txt")) %>%
    ensurer::ensure_that(all(file.exists(.$label_file)),
                         err_desc = "Missing label file") %>%
    dplyr::mutate(raster_obj = purrr::map(file_path, raster::raster)) %>%
    dplyr::left_join(samples_tb, by = "used_bands") %>%
    dplyr::mutate(ref_pred = purrr::pmap(dplyr::select(., samples, raster_obj),
                                         get_ref_pred, label_vec = my_recode2))

# Print confusion matrices.
table(results_tb$ref_pred[[1]])
table(results_tb$ref_pred[[2]])

resampled <- list()
resampled[[1]] <- results_tb$ref_pred[[1]] %>%
    dplyr::group_by(predicted) %>%
    dplyr::sample_n(84) %>%
    dplyr::ungroup()
resampled[[2]] <- results_tb$ref_pred[[2]] %>%
    dplyr::group_by(predicted) %>%
    dplyr::sample_n(84) %>%
    dplyr::ungroup()
results_tb[["resampled"]] <- resampled

# Print confusion matrices.
table(results_tb$resampled[[1]])
table(results_tb$resampled[[2]])

# Print accuracy
cmt_bands   <- as.matrix(table(results_tb$resampled[[1]]))
cmt_indices <- as.matrix(table(results_tb$resampled[[2]]))
# Overall
sum(diag(cmt_bands))/sum(colSums(cmt_bands))
sum(diag(cmt_indices))/sum(colSums(cmt_indices))
# Producer
diag(cmt_bands)/rowSums(cmt_bands)
diag(cmt_indices)/rowSums(cmt_indices)
# User
diag(cmt_bands)/colSums(cmt_bands)
diag(cmt_indices)/colSums(cmt_indices)

# Compute statistics
results_tb <- results_tb %>%
    dplyr::mutate(con_mat = purrr::map(resampled, sits::sits_conf_matrix))

# Write confusion matrices to disc.
results_tb %>%
    dplyr::select(file_path, con_mat) %>%
    dplyr::mutate(out_file = tools::file_path_sans_ext(file_path),
                  out_file = stringr::str_c(out_file,
                                            "_confusion_matrix_rows_are_reference.csv",
                                            sep = '')) %>%
    dplyr::mutate(cm = purrr::map(con_mat, magrittr::extract2, "table"),
                  cm = as.matrix(cm)) %>%
    (function(x){
         for (i in 1:nrow(x)){
             print(sprintf("Writing confusion matrix to %s", x$out_file[[i]]))
            write.csv(x$cm[[i]], file =  x$out_file[[i]])
         }
    })

results_tb %>%
    dplyr::mutate(def_pa = purrr::map_dbl(con_mat, get_up_accuracy, label = "Deforestation", acc_type = "pa"),
                  def_ua = purrr::map_dbl(con_mat, get_up_accuracy, label = "Deforestation", acc_type = "ua"),
                  def_f1 = purrr::map_dbl(con_mat, get_up_accuracy, label = "Deforestation", acc_type = "F1"),
                  for_pa = purrr::map_dbl(con_mat, get_up_accuracy, label = "Forest",        acc_type = "pa"),
                  for_ua = purrr::map_dbl(con_mat, get_up_accuracy, label = "Forest",        acc_type = "ua"),
                  for_f1 = purrr::map_dbl(con_mat, get_up_accuracy, label = "Forest",        acc_type = "F1"),
                  nfo_pa = purrr::map_dbl(con_mat, get_up_accuracy, label = "Non-Forest",    acc_type = "pa"),
                  nfo_ua = purrr::map_dbl(con_mat, get_up_accuracy, label = "Non-Forest",    acc_type = "ua"),
                  nfo_F1 = purrr::map_dbl(con_mat, get_up_accuracy, label = "Non-Forest",    acc_type = "F1")) %>%
    dplyr::select(used_bands, def_pa, def_ua, for_pa, for_ua, nfo_pa, nfo_ua) %>%
    write.csv(file.path(result_dir, "accuracies.csv"))

