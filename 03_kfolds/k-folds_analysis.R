#!/usr/bin/Rscript
suppressMessages(library(dplyr))
suppressMessages(library(caret))
suppressMessages(library(sits))

source("./other/util.R")

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("This script takes parameters: An input file (RDS) of time series of samples, and an output directory.",  call. = FALSE)
}

rds_file     <- args[[1]]
out_base_dir <- args[[2]]

out_dir <- out_base_dir %>%
    file.path(tools::file_path_sans_ext(basename(rds_file)))
if (!dir.exists(out_dir))
   dir.create(out_dir, recursive = TRUE)

samples_tb <- rds_file %>%
    readRDS() %>%
    dplyr::filter(label != "Water")

# NOTE: save to PANGEA
samples_tb %>%
    dplyr::mutate(id = dplyr::row_number()) %>%
    dplyr::select(id, longitude, latitude, start_date, end_date, label) %>%
    readr::write_csv("./pangea/training_dataset.csv")

experiments <- list(all_bands = c("blue","bnir","green","nnir","red","swir1","swir2"),
                    indeces   = c("evi","ndmi","ndvi"))
warning(sprintf("Using ONLY %s combination of bands:\n%s", length(experiments),
                paste(lapply(experiments, paste, collapse = "-"),
                      collapse = "\n")))


run_kfolds <- function(bands, samples_tb){
    lapply(1:100, function(i){
        samples_tb %>%
            sits::sits_sample(n = 60) %>%
            select_bands(bands) %>%
            sits::sits_kfold_validate(folds = 10,
                                      ml_method = sits::sits_rfor(num_trees = 1000)) %>%
            sits::sits_conf_matrix() %>%
            return()
    })
}

kfold_ls <- lapply(experiments, run_kfolds, samples_tb = samples_tb)

helper_acc <- function(x, label, acc_type){
    sapply(x,  get_up_accuracy, label = label, acc_type = acc_type)
}

exp_tb <- tibble::tibble(experiment = names(experiments)) %>%
    dplyr::mutate(kfold = kfold_ls) %>%
    dplyr::mutate(def_pa = purrr::map(kfold, helper_acc, label = "Deforestatio", acc_type = "pa"),
                  def_ua = purrr::map(kfold, helper_acc, label = "Deforestatio", acc_type = "ua"),
                  for_pa = purrr::map(kfold, helper_acc, label = "Forest",       acc_type = "pa"),
                  for_ua = purrr::map(kfold, helper_acc, label = "Forest",       acc_type = "ua"),
                  nnf_pa = purrr::map(kfold, helper_acc, label = "NatNonForest", acc_type = "pa"),
                  nnf_ua = purrr::map(kfold, helper_acc, label = "NatNonForest", acc_type = "ua"),
                  nof_pa = purrr::map(kfold, helper_acc, label = "NonForest",    acc_type = "pa"),
                  nof_ua = purrr::map(kfold, helper_acc, label = "NonForest",    acc_type = "ua"),
                  pas_pa = purrr::map(kfold, helper_acc, label = "Pasture",      acc_type = "pa"),
                  pas_ua = purrr::map(kfold, helper_acc, label = "Pasture",      acc_type = "ua")) %>%
    # NOTE: Remove embedded NA columns.
    dplyr::select_if(function(x){
        if(!is.list(x))
            return(TRUE)
        res <- list()
        for(i in 1:length(x)){
            res[[i]] <- !all(is.na(x[[i]]))
        }
        return(all(res))
    })

helper_plot <- function(x, label){
    x %>%
        tibble::as_tibble() %>%
        dplyr::mutate(index = 1:dplyr::n()) %>%
        tidyr::pivot_longer(cols = !tidyselect::matches("index"),
                            names_to = "experiment") %>%
        dplyr::mutate(experiment = tools::toTitleCase(stringr::str_replace(experiment, '_', ' '))) %>%
        ggplot2::ggplot(ggplot2::aes(x = value, y = experiment)) +
        ggplot2::geom_boxplot() +
        ggplot2::xlab("Index") +
        ggplot2::ylab("Band Combination") +
        ggplot2::theme(text = ggplot2::element_text(size = 8)) +
        ggplot2::xlim(0.4, 1.0) +
        ggplot2::ggtitle(label)
}

fix_name <- function(x){
    x %>%
        stringr::str_replace('_', ' ') %>%
        stringr::str_replace(" ua", " User Accuracy") %>%
        stringr::str_replace(" pa", " Producer Accuracy") %>%
        stringr::str_replace("def ", "Deforestation ") %>%
        stringr::str_replace("for ", "Forest ") %>%
        stringr::str_replace("nnf ", "Natural Non-Forest ") %>%
        stringr::str_replace("nof ", "Non-Forest ") %>%
        stringr::str_replace("pas ", "Pasture ") %>%
        return()
}

for (name in names(exp_tb)[3:ncol(exp_tb)]) {
    fixed_name <- fix_name(name)
    helper_plot(exp_tb[[name]], fixed_name) %>%
         (function(x){
             plot(x)
             invisible(x)
         })
        ggplot2::ggsave(file = file.path(out_dir,
                                         paste0(stringr::str_replace_all(tolower(fixed_name),
                                                                  pattern = ' ',
                                                                  replacement = '_'),
                                         ".png")),
                        width = 8, height = 3, units = "cm")
}
