# Remove invalid samples of time series.
#
# @param  sits_tb A sits_tibble.
# @report report  When TRUE, not cleaning is done, just marking the offending samples.
# @return A sits_tibble.
clean_ts <- function(sits_tb, report = FALSE){
    sits_tb %>%
        sits::sits_prune() %>%
        tidyr::drop_na() %>%
        dplyr::mutate(has_na    = purrr::map_int(time_series, function(x){return(sum(is.na(x)))}),
                      has_null  = purrr::map_int(time_series, function(x){return(sum(is.null(x), na.rm = TRUE))}),
                      has_overflow  = purrr::map_int(time_series, function(x){return(sum(sum(as.matrix(x[,2:ncol(x)]) < -1, na.rm = TRUE), sum(as.matrix(x[,2:ncol(x)]) > 1, na.rm = TRUE)))}),
                      time_mean = purrr::map_dbl(time_series, function(x){return(mean(x[[1]]))}),
                      n_cols    = purrr::map_int(time_series, ncol),
                      n_rows    = purrr::map_int(time_series, nrow)) %>%
        (function(.data){
             if (report){
                 return(.data)
             }else{
                 .data <- .data %>%
                     tidyr::drop_na() %>%
                     dplyr::filter(!has_null,
                                   n_cols > 1,
                                   n_rows > 0) %>%
                     dplyr::mutate(time_series = purrr::map(time_series, function(x){
                         my_approx <- function(v) {
                             apply(v, 2,
                             function(x) {
                                 i <- tryCatch({
                                     approx(x, n = length(x))
                                 }, error = function(e) list(y = rep(0, length(x))))
                                 return(i$y)
                             })
                         }
                         data_mt <- as.matrix(x[,2:ncol(x)])
                         data_mt[data_mt <= -1] <- NA
                         data_mt[data_mt >= 1]  <- NA
                         interp_mt <- my_approx(data_mt)
                         x %>%
                             dplyr::select(Index) %>%
                             dplyr::bind_cols(tibble::as_tibble(interp_mt)) %>%
                             return()
                     })) %>%
                     dplyr::select(-has_na, -has_null, -time_mean,
                                   -overflow, -n_cols, -n_rows)
                 n_removed <- nrow(sits_tb) - nrow(.data)
                 if (n_removed > 0)
                     warning(sprintf("Removed %s invalid samples out of  %s",
                                     n_removed, nrow(sits_tb)))
                 return(.data)
             }
        }) %>%
        return()
}


# Compute all the combinations of n of the given bands.
#
# @param .data           n_bands. An integer.
# @param available_bands A character. The names of the bands.
# @return                A list of vectors.
combine_n_bands <- function(n_bands, available_bands){
    available_bands %>%
        combn(m = n_bands) %>%
        split(f = col(.)) %>%
        return()
}



# Compute vegetation indexes of Sentinel images.
#
# @param vrt_file A lenght-one character. Path to a VRT file of Sentinel-2 10 meter bands B02, B03, B04, and B08.
# @param out_file A lenght-one character. Path to a file.
# @param out_file A lenght-one character. Short name of a index.
# @return    A character. A temporal file.
compute_vi_sentinel <- function(vrt_file, out_file, index_name){
    .Deprecated("Use bash scripts instead.")
    # VRT file:
    # - B02 blue   1
    # - B03 green  2
    # - B04 red    3
    # - B08 bnir   4
    # - B11 swir1  5
    index_name <- toupper(index_name)
    cmd <- list()
    cmd[["EVI"]] <- sprintf("gdal_calc.py -A %s --A_band=4 -B %s --B_band=3 -C %s --C_band=1 --outfile=%s --calc='(10000 * 2.5 * (A - B) / (A + 6.0 * B - 7.5 * C + 10000.001)).astype(int16)' --NoDataValue=-9999 --type='Int16' --creation-option='COMPRESS=LZW' --creation-option='BIGTIFF=YES'", vrt_file, vrt_file, vrt_file, out_file)
    cmd[["NDMI"]] <- sprintf("gdal_calc.py -A %s --A_band=4 -B %s --B_band=5 --outfile=%s --calc='((A - B)/(A + B + 0.001)*10000).astype(int16)' --NoDataValue=-9999 --type='Int16' --creation-option='COMPRESS=LZW' --creation-option='BIGTIFF=YES'", vrt_file, vrt_file, out_file)
    cmd[["NDVI"]] <- sprintf("gdal_calc.py -A %s --A_band=4 -B %s --B_band=3 --outfile=%s --calc='((A - B)/(A + B + 0.001)*10000).astype(int16)' --NoDataValue=-9999 --type='Int16' --creation-option='COMPRESS=LZW' --creation-option='BIGTIFF=YES'", vrt_file, vrt_file, out_file)
    cmd[["SAVI"]] <- sprintf("gdal_calc.py -A %s --A_band=4 -B %s --B_band=3 --outfile=%s --calc='((A - B)/(A + B + 4280.0011) * (10000 + 4280)).astype(int16)' --NoDataValue=-9999 --type='Int16' --creation-option='COMPRESS=LZW' --creation-option='BIGTIFF=YES'", vrt_file, vrt_file, out_file)
    cmd[["MTVI"]]  <- sprintf("gdal_calc.py -A %s --A_band=4 -B %s --B_band=3 -C %s --C_band=2 --outfile=%s --calc='(1.2 * (1.2 * (A - C) - 2.5 * (B - C))).astype(int16)' --NoDataValue=-9999 --type='Int16' --creation-option='COMPRESS=LZW' --creation-option='BIGTIFF=YES'", vrt_file, vrt_file, vrt_file, out_file)
    cmd[["OSAVI"]] <- sprintf("gdal_calc.py -A %s --A_band=4 -B %s --B_band=3 --outfile=%s --calc='((1 + 0.16) * (A - B)/(A + B + 1600.0001)*10000).astype(int16)' --NoDataValue=-9999 --type='Int16' --creation-option='COMPRESS=LZW' --creation-option='BIGTIFF=YES'", vrt_file, vrt_file, out_file)
    cmd[["RDVI"]]  <- sprintf("gdal_calc.py -A %s --A_band=4 -B %s --B_band=3 --outfile=%s --calc='((A - B)/((A + B)/10000.001)**0.5).astype(int16)' --NoDataValue=-9999 --type='Int16' --creation-option='COMPRESS=LZW' --creation-option='BIGTIFF=YES'", vrt_file, vrt_file, out_file)
    cmd[["RDVI2"]] <- cmd[["RDVI"]]
    #cmd[["PC1RGBNIR"]] <- sprintf("gdal_calc.py -A %s --A_band=1 -B %s --B_band=2 -C %s --C_band=3 -D %s --D_band=4 --outfile=%s --calc='(0.06379167*A + 0.14536839*B + 0.04085506*C + 0.98647327*D).astype(int16)' --NoDataValue=-9999 --type='Int16' --creation-option='COMPRESS=LZW' --creation-option='BIGTIFF=YES'", vrt_file, vrt_file, vrt_file, vrt_file, out_file)
    #cmd[["PC2RGBNIR"]] <- sprintf("gdal_calc.py -A %s --A_band=1 -B %s --B_band=2 -C %s --C_band=3 -D %s --D_band=4 --outfile=%s --calc='(-0.4401453*A + -0.5317655*B + -0.7105860*C +  0.1362536*D).astype(int16)' --NoDataValue=-9999 --type='Int16' --creation-option='COMPRESS=LZW' --creation-option='BIGTIFF=YES'", vrt_file, vrt_file, vrt_file, vrt_file, out_file)
    # WRONG cmd[["EVI2"]] <- sprintf("gdal_calc.py -A %s --A_band=4 -B %s --B_band=3 --outfile=%s --calc='(24000 * (A - B)/(A + B + 10000.001)).astype(int16)' --NoDataValue=-9999 --type='Int16' --creation-option='COMPRESS=LZW' --creation-option='BIGTIFF=YES'", vrt_file, vrt_file, out_file)
# GEMI Cao et al 2009 Pinty and verstraete 1992
# n = ((2/10000 * (A^2 - B^2) + 1.5*A + 0.5*B)/(A + B + 5000.001))
# GEMI = (n*(1 - 0.25*n) - (B - 1250)/(10000.0001 - B)) * 10000
    # WRONG cmd[["GEMI"]] <- sprintf("gdal_calc.py -A %s --A_band=4 -B %s  --B_band=3 --outfile=%s --calc='((((2/10000 * (A*A - B*B) + 1.5*A + 0.5*B)/(A + B + 5000.001))*(1 - 0.25*((2/10000 * (A*A - B*B) + 1.5*A + 0.5*B)/(A + B + 5000.001))) - (B - 1250)/(10000.0001 - B)) * 10000).astype(int16)' --NoDataValue=-9999 --type='Int16' --creation-option='COMPRESS=LZW' --creation-option='BIGTIFF=YES'", vrt_file, vrt_file, out_file)
    #
    # SAVI https://github.com/sentinel-hub/custom-scripts/blob/master/sentinel-2/savi/script.js
    # NDWI https://doi.org/10.1080/01431169608948714
    #
    if(index_name %in% names(cmd) == FALSE) {
        stop(sprintf("Unknown index: %s", index_name))
    }
    res <- system(cmd[[index_name]])
    invisible(out_file)
}


# Count the number of JP2 files recursively in the given directory.
#
# @param x A path to a directory.
# @return  A character.
count_jp2 <- function(x){
    x %>%
        list.files(pattern = ".jp2$", recursive = TRUE) %>%
        length() %>%
        return()
}


# Helper function for copying images from a tibble of satellite images.
#
# @param .data   A tibble.
# @param col     A variable in the given tibble.
# @param out_dir A lenght-one character. Desitination directory.
# @return out_dir (invisible).
cp_tb_files <- function(.data, col, out_dir){
    col <- rlang::enquo(col)
    .data %>%
        dplyr::pull(!!col) %>%
        file.copy(to = out_dir, recursive = TRUE, overwrite = FALSE)
    invisible(out_dir)
}


# Helper function to ensure distinct images.
#
# @param x A tibble with at least the fields landsat_date and tier..
# @return  A tibble.
distinct_landsat <- function(x){
    x %>%
        dplyr::distinct(landsat_date, tier, .keep_all = TRUE) %>%
        return()
}


# Helper function to filter a tibble of sentinel images by a band and then pile those images.
#
# @param band_name A length-one character. Name of a Sentinel2 band.
# @param img_tb    A tibble of Sentinel-2 images.
# @param out_pattern A char
# @return
filter_and_pile <- function(band_name, img_tb, out_pattern){
    .Deprecated("Use bash scripts instead.")
    band_tb <- img_tb %>%
        dplyr::filter(band == band_name) %>%
        dplyr::arrange(acquisition)
    out_file <- paste(out_pattern,
                      stringr::str_sub(band_tb$acquisition[[1]], 1, 8),
                      band_name, paste0(unique(band_tb$resolution), ".tif"),
                      sep = '_')
    band_tb %>%
        dplyr::pull(file_path) %>%
        pile_files(out_fn = out_file)
    invisible(out_file)
}


# Get the user or producer accuracies from an accuracy object.
# @param x        An accuracy object.
# @param label    A lenght-one character. A label(or class name).
# @param acc_type A lenght-one character. The type of accuracy: Producer (pa) or User (ua)
# @return         A matrix or NA.
get_up_accuracy <- function(x, label, acc_type = "pa") {
    stopifnot(acc_type %in% c("pa", "ua", "F1"))
    accuracy_label <- acc_type
    if (acc_type == "pa")
        accuracy_label <- "Sensitivity"
    if (acc_type == "ua")
        warning("Check if Pos Pred Value is the same as User Accuracy")
        accuracy_label <- "Pos Pred Value"
    index_mt <- x %>%
        magrittr::extract2("byClass")
    stopifnot("Sensitivity" %in% colnames(index_mt))
    label <- paste0(' ', label)
    col_id <- match(accuracy_label, stringr::str_match(colnames(index_mt), accuracy_label))
    row_id <- match(toupper(label), stringr::str_match(toupper(rownames(index_mt)), toupper(label)))
    if (any(is.na(c(row_id, col_id))))
        return(NA)
    return(index_mt[row_id, col_id])
}


# Get the dimension of an image using gdal_info.
#
# @param in_file A character.
# @return        A integer or a list.
get_img_dimensions <- function(in_file) {
    stopifnot(is.atomic(in_file))
    if (is.na(in_file) || length(in_file) < 1)
        return(NA)
    if (length(in_file) == 1) {
        system2("gdalinfo", in_file, stdout = TRUE) %>%
        stringr::str_subset("Size is ") %>%
        stringr::str_split(" ") %>%
        unlist() %>%
        stringr::str_remove_all(pattern = ",") %>%
        utils::tail(2) %>%
        ensurer::ensure_that(length(.) > 1) %>%
        as.integer() %>%
        magrittr::set_names(c("pixels", "lines")) %>%
        return()
    } else {
        return(vapply(in_file, get_img_dimensions, integer(2)))
    }
}


# Get metadata of the Sentinel-2 bricks.
#
# @param in_dir          A length-one character. Path to a directory.
# @return                A tibble.
get_brick <- function(in_dir){
    .Deprecated("")
    in_dir %>%
        get_brick_md() %>%
        dplyr::arrange(tile, img_date, band) %>%
        ensurer::ensure_that(length(unique(.$img_date)) == 1,
                             err_desc = sprintf("More than one date found: %s",
                                                in_dir)) %>%
        ensurer::ensure_that(length(unique(.$tile)) == 1,
                             err_desc = sprintf("More than one tile found: %s",
                                                in_dir)) %>%
        ensurer::ensure_that(!"" %in% .$band) %>%
        return()
}


# Get the metadata of the Sentinel-2 bricks in a directory.
#
# @param in_dir A length-one character. Path to a directory.
# @return       A tibble.
get_brick_md <- function(in_dir){
    file_band_names <- tibble::tribble(
        ~file_band, ~band,
        "B02",        "blue",
        "B03",        "green",
        "B04",        "red",
        "B08",        "bnir",
        "B8A",        "nnir",
        "B11",        "swir1",
        "B12",        "swir2",
        "evi",        "evi",
        "evi2",       "evi2",
        "gemi",       "gemi",
        "mtvi",       "mtvi",
        "ndmi",       "ndmi",
        "ndvi",       "ndvi",
        "ndwi",       "ndwi",
        "osavi",      "osavi",
        "rdvi",       "rdvi",
        "savi",       "savi",
        "pc1rgbnir",  "pc1rgbnir",
        "pc2rgbnir",  "pc2rgbnir"
    )
    in_dir %>%
        list.files(pattern = ".[.](tif|vrt)$", full.names = TRUE) %>%
        tibble::enframe(name = NULL) %>%
        dplyr::rename(file_path = value) %>%
        dplyr::mutate(file_name = tools::file_path_sans_ext(basename(file_path))) %>%
        dplyr::mutate(file_ext = tools::file_ext(basename(file_path))) %>%
        tidyr::separate(col = file_name, into = c("mission", "level", "orbit",
                                                  "tile", "img_date", "file_band",
                                                  "brick_type",
                                                  "resolution"),
                        sep = '_', extra = "drop", fill = "right") %>%
        dplyr::mutate(resolution = ifelse(brick_type %in% c("10m", "20m", "60m"),
                                          brick_type, resolution),
                      brick_type = ifelse(brick_type %in% c("10m", "20m", "60m"),
                                               "", brick_type),
                      brick_type = ifelse(brick_type == "", "raw", brick_type),
                      img_date = lubridate::as_date(stringr::str_sub(img_date,
                                                                     1, 8))) %>%
        dplyr::left_join(file_band_names, by = "file_band") %>%
        return()
}


# Mask a sentinel image.
#
# @param file_path  A length-one character. Path to a file.
# @param fmask_path A length-one character. Path to a mask file.
# @param out_dir    A leghth-one character. Path to a directory.
# @return           A character. Path to a file in out_dir.
mask_sentinel <- function(file_path, fmask_path, out_dir) {
    .Deprecated("Use bash scripts instead.")
    out_file <- file_path %>%
        basename() %>%
        #tools::file_path_sans_ext() %>%
        #stringr::str_c("_masked.tif") %>%
        (function(.x){ return(file.path(out_dir, .x))})
    vrt_file <- c(file_path, fmask_path) %>%
        gdalcmdline::gdal_build_vrt(out_filename = tempfile(pattern = "gdalbuildvrt_",
                                                          fileext = ".vrt"),
                                    resolution = "highest", separate = TRUE,
                                    vrtnodata = -9999)
    cmd <- sprintf("gdal_calc.py -A %s --A_band=1 -B %s --B_band=2 --outfile=%s --calc='(numpy.where(B != 4, A, -9999)).astype(int16)' --type='Int16' --NoDataValue=-9999 --creation-option='COMPRESS=LZW'",
                   vrt_file, vrt_file, out_file)
    res <- system(cmd)
    return(out_file)
}


# Pile the given files usind gdal_merge.
#
# @param file_paths   A character. Path to files.
# @param out_fn       A length-one character. Path to the output file.
# @param gdal_format  A length-one character.
# @param no_data      A numeric.
# @param gdal_options A character.
# @return  A tibble.
pile_files <- function (file_paths, out_fn, gdal_format = "GTiff",
                        no_data = -9999, data_type = "Int16",
                        gdal_options = c("TILED=YES", "COMPRESS=LZW", "BIGTIFF=YES")){
    .Deprecated("Use bash scripts instead.")
    gdalcmdline::gdal_merge(input_files = file_paths, out_filename = out_fn,
                            ot = data_type, separate = TRUE, of = gdal_format,
                            creation_option = gdal_options, init = no_data,
                            a_nodata = no_data) %>%
    invisible()
}


# Helper function for piling Sentinel2 images.
#
# @param mission     A length-one character.
# @param level       A length-one character.
# @param orbit       A length-one character.
# @param tile        A length-one character
# @param pyear       A lenght-one integer. The PRODES year.
# @param sentinel_tb A tibble of Sentinel-2 images.
# @param out_dir     A length-one character. A path to a directory.
# @return  A tibble.
pile_sentinel_images <- function(mission, level, orbit, tile, pyear, sentinel_tb,
                        field_name,
                        out_dir){
    .Deprecated("Use bash scripts instead.")
    #field_name <- rlang::enquo(field_name)
    img_tb <- sentinel_tb %>%
        dplyr::filter(mission == mission, level == level, tile == tile,
                      pyear == pyear) %>%
        dplyr::select(!!field_name) %>%
        #dplyr::select(files_10m) %>%
        tidyr::unnest(!!field_name)
        # tidyr::unnest(files_10m)
    band_names <- img_tb %>%
        dplyr::pull(band) %>%
        unique() %>%
        sort()
    brick_path <- purrr::map_chr(band_names, filter_and_pile, img_tb = img_tb,
                                 out_pattern = file.path(out_dir,
                                                         stringr::str_c(mission, level, orbit, tile, sep = '_')))
    invisible(brick_path)
}


# Compute the PRODES year of the given date.
#
# @param x A lubridate's date object.
# @return  A integer.
prodes_year <- function(x){
    if(!lubridate::is.Date(x))
        return(NA)
    if(length(x) == 1){
        m <- lubridate::month(x)
        if(m < 8)
            return(as.integer(lubridate::year(x)))
        return(as.integer(lubridate::year(x) + 1))
    }else if(length(x) > 1){
        vapply(x, prodes_year, integer(1))
    }
}


#  Select some bands in a sits tibble using a character.
#
# @param  x A sits tibble.
# @param  bands A character. The name of some bands in the time series column of the given sits tibble.
# @return A sits tibble.
select_bands <- function(x, bands){
    bands <- c("Index", bands)
    x$time_series <- lapply(x$time_series, function(y){
                                if(all(bands %in% colnames(y)))
                                    return(y[, bands])
                                else
                                    return(NA)
                                 })
    return(x)
}


















# Split an image into chunks.
#
# @param in_file A character. Path to a file.
# @param x_size  An integer. Size of chunks in the x-direction.
# @param y_size  An integer. Size of chunks in the y-direction.
# @param out_dir A character. Path to a directory.
# @return        A tibble.
split_image <- function(in_file, xsize = 256, ysize = 256, out_dir = tempdir()){
    .Deprecated()
    # Chunk size = 64  bits * 36 images * 64 * 64 = 1.125 Mebibytes
    # Chunk size = 128 bits * 36 images * 64 * 64 = 2.25  Mebibytes
    # Chunk size = 256 bits * 36 images * 64 * 64 = 4.5   Mebibytes
    stopifnot(length(in_file) == 1)
    stopifnot(file.exists(in_file))
    stopifnot(dir.exists(out_dir))
    img_dim <- in_file %>%
        get_img_dimensions()
    pixels <- seq(from = 0, to = img_dim["pixels"], by = xsize)
    lines  <- seq(from = 0, to = img_dim["lines"],  by = ysize)
    chunk_tb <- expand.grid(pixels, lines) %>%
        tibble::as_tibble() %>%
        dplyr::rename(pixel_from = Var1, line_from = Var2) %>%
        dplyr::mutate(x_win = ifelse(img_dim["pixels"] < pixel_from + xsize, img_dim["pixels"] - pixel_from, xsize),
                      y_win = ifelse(img_dim["lines"]  < line_from  + ysize, img_dim["lines"]  - line_from,  ysize),
                      out_file = stringr::str_c(tools::file_path_sans_ext(basename(in_file)),
                                                pixel_from, line_from, sep = '_'),
                      out_file = paste0(file.path(out_dir, out_file), ".tif")) %>%
        dplyr::mutate(command = stringr::str_c("gdal_translate -q -of GTiff -srcwin",
                                               pixel_from, line_from,
                                               x_win, y_win,
                                               in_file, out_file,
                                               sep = ' ')) %>%
        dplyr::mutate(command_res = purrr::map_int(command, system)) %>%
        ensurer::ensure_that(all(dplyr::pull(., command_res) == 0),
                             err_desc = sprintf("Image chunking failed for file %s",
                                                in_file)) %>%
        dplyr::select(-command, -command_res) %>%
        return()
}



#--- Helper functions ----


# Helper for masking in parallel.
#
# @param id_row      A length-one integer. A row number in sentinel_tb.
# @param sentinel_tb A tibble describing Sentienl-2 images.
# @return            A character. Path to the masked image.
helper_mask <- function(id_row, sentinel_tb, out_dir = tempdir()){
    .Deprecated("helper_mask2")
    stopifnot(id_row %in% 1:nrow(sentinel_tb))
    in_files <- sentinel_tb %>%
        dplyr::slice(id_row) %>%
        dplyr::select(file_path, fmask_path) %>%
        unlist()
    in_files["file_path"] %>%
             mask_sentinel(fmask_path = in_files["fmask_path"],
                           out_dir = out_dir) %>%
             return()
}

# Helper for masking in parallel.
#
# @param id_row      A length-one integer. A row number in sentinel_tb.
# @param sentinel_tb A tibble describing Sentinel-2 images.
# @param var         A name of a variable (column) in sentinel_tb.
# @return            A character. Path to the masked image.
helper_mask2 <- function(id_row, img_tb, var, out_dir = tempdir()){
    .Deprecated("Use bash scripts instead.")
    stopifnot(id_row %in% 1:nrow(sentinel_tb))
    var <- dplyr::enquo(var)
    in_files <- img_tb %>%
        dplyr::slice(id_row) %>%
        dplyr::select(!!var, fmask_path) %>%
        unlist()
    in_files[1] %>%
             mask_sentinel(fmask_path = in_files["fmask_path"],
                           out_dir = out_dir) %>%
             return()
}






# DEPRECATED
helper_pile <- function(x, out_dir){
    .Deprecated()
    out_file <- x %>%
        (function(.data){
            .data %>%
                dplyr::select(mission, level, baseline, orbit, pyear, tile, band, resolution, pixel_from, line_from) %>%
                lapply(., unique) %>%
                sapply(., length) %>%
                ensurer::ensure_that(all(. == 1),
                                     err_desc = "Invalid image tibble.")
            return(.data)
        }) %>%
        dplyr::arrange(img_date) %>%
        dplyr::slice(1) %>%
        dplyr::pull(safe_path) %>%
        basename() %>%
        tools::file_path_sans_ext() %>%
        paste(unique(x$band), unique(x$pixel_from), paste0(unique(x$line_from),
                                                           ".tif"), sep = '_')
    out_file <- file.path(out_dir, out_file)
    x %>%
        dplyr::pull(file_path) %>%
        pile_files(out_fn = out_file, gdal_format = "GTiff", no_data = -9999,
                   gdal_options = c("TILED=YES", "COPY_SRC_OVERVIEWS=YES",
                                         "COMPRESS=LZW", "BIGTIFF=YES"))
    return(out_file)
}


# Helper function for piling up Sentinel-2 images.
#
# @param .x      A tibble of 36 rows and at least the columns img_date (date), mission, level, orbit, tile, acquisition, band, resolution, file_path.
# @param out_dir A length-one character. Path to a directory.
# @return        A tibble.
helper_pile2 <- function(.x, out_dir){
    .Deprecated("Use bash scripts instead.")
    file_tb <- .x %>%
        ensurer::ensure_that(nrow(.) == 36, err_desc = "Missing images!") %>%
        dplyr::arrange(img_date)
    out_file <- file_tb %>%
        dplyr::slice(1) %>%
        dplyr::select(mission, level, orbit, tile, acquisition, band, resolution) %>%
        unlist() %>%
        paste(collapse = '_') %>%
        paste0(".tif")
    out_file <- file.path(out_dir, out_file)
    file_tb %>%
        dplyr::pull(file_path) %>%
        ensurer::ensure_that(length(.) > 0,
                             err_desc = "Missing path to images!") %>%
        pile_files(out_fn = out_file)
    file_tb %>%
        dplyr::slice(1) %>%
        dplyr::select(-safe_path, -file_path, -acquisition) %>% #, -processing
        dplyr::mutate(brick_file = out_file) %>%
        return()
}


# Helper function for piling up Sentinel-2 images.
#
# @param .x      A tibble of 36 rows and at least the columns img_date (date), mission, level, orbit, tile, acquisition, band, resolution, file_path.
# @param out_dir A length-one character. Path to a directory.
# @return        A tibble.
helper_pile_masked <- function(.x, out_dir){
    .Deprecated("Use bash scripts instead.")
    file_tb <- .x %>%
        ensurer::ensure_that(nrow(.) == 36, err_desc = "Missing images!") %>%
        dplyr::arrange(img_date)
    out_file <- file_tb %>%
        dplyr::slice(1) %>%
        dplyr::select(mission, level, orbit, tile, acquisition, band, resolution) %>%
        unlist() %>%
        paste(collapse = '_') %>%
        paste0(".tif")
    out_file <- file.path(out_dir, out_file)
    file_tb %>%
        dplyr::pull(img_masked) %>%
        pile_files(out_fn = out_file)
    file_tb %>%
        dplyr::slice(1) %>%
        dplyr::select(-safe_path, -file_path, -acquisition, -processing) %>%
        dplyr::mutate(brick_file = out_file) %>%
        return()
}


# Helper function for piling up Sentinel-2 images.
#
# @param .x      A tibble of 36 rows and at least the columns img_date (date), mission, level, orbit, tile, acquisition, band, resolution, file_path.
# @param out_dir A length-one character. Path to a directory.
# @return        A tibble.
helper_pile_raw <- function(.x, out_dir){
    .Deprecated("helper_pile2")
    file_tb <- .x %>%
        ensurer::ensure_that(nrow(.) == 36, err_desc = "Missing images!") %>%
        dplyr::arrange(img_date)
    out_file <- file_tb %>%
        dplyr::slice(1) %>%
        dplyr::select(mission, level, orbit, tile, acquisition, band, resolution) %>%
        unlist() %>%
        paste(collapse = '_') %>%
        paste0(".tif")
    out_file <- file.path(out_dir, out_file)
    file_tb %>%
        dplyr::pull(file_path) %>%
        pile_files(out_fn = out_file)
    file_tb %>%
        dplyr::slice(1) %>%
        dplyr::select(-safe_path, -file_path, -acquisition, -processing) %>%
        dplyr::mutate(brick_file = out_file) %>%
        return()
}


# Helper for computing vegetation indexes.
#
# @param vrt_file A lenght-one character. Path to a VRT file of a Sentinel-2 image.
# @return         A tibble.
helper_vi <- function(vrt_file){
    .Deprecated("Use bash scripts instead.")
    #vi_names <- c("evi", "ndmi", "ndvi", "savi", "mtvi", "osavi", "rdvi")
    #vi_names <- c("evi", "ndmi", "ndvi", "savi")
    vi_names <- c("evi", "ndmi")
    out_files <- vi_names %>%
        paste0('_') %>%
        lapply(tempfile, fileext = ".tif") %>%
        magrittr::set_names(vi_names)
    for(vi in vi_names) {
        print(vi)
        compute_vi_sentinel(vrt_file, out_file = out_files[[vi]],
                            index_name = vi)
    }
    out_files %>%
        tibble::as_tibble() %>%
        return()
}


# Helper for bulding a VRT for ease vegetation index computation of Sentienl-2 images.
#
# @param B02 A lenght-one character. Path to a file of Sentinel-2 band
# @param B03 A lenght-one character. Path to a file of Sentinel-2 band.
# @param B04 A lenght-one character. Path to a file of Sentinel-2 band.
# @param B08 A lenght-one character. Path to a file of Sentinel-2 band.
# @return    A character. Path to a VRT file.
helper_vrt_vi <- function(B02, B03, B04, B08, B11){
    .Deprecated("Use bash scripts instead.")
    c(B02, B03, B04, B08, B11) %>%
        gdalcmdline::gdal_build_vrt(out_filename = tempfile(pattern = "sentinel_b2-3-4-8-11_",
                                                            fileext = ".vrt"),
                                    resolution = "highest", separate = TRUE,
                                    vrtnodata = -9999) %>%
        return()
}


# Run PCA in sample time series using a subset of the available bands and labels.
do_pca <- function(selected_bands, selected_labels, samples_tb){
    samples_tb %>%
        dplyr::filter(label %in% selected_labels) %>%
        select_bands(selected_bands) %>%
        ensurer::ensure_that(nrow(.) > 0, err_desc = "Samples are missing") %>%
        dplyr::select(label, time_series) %>%
        tidyr::unnest(time_series) %>%
        dplyr::select(-label, -Index) %>%
        stats::prcomp(center = TRUE, scale = rep(10000, length(selected_bands))) %>%
        return()
}

# Get the rotation matrix from a PCA object.
get_pca_rotation <- function(pca, components = 1:2){
    pca %>%
        magrittr::extract2("rotation") %>%
        magrittr::extract(, components) %>%
        return()
}

# Compute and extract the PCA bands in a sits tibble
get_pca_bands <- function(pca_rotation, samples_tb){
    samples_tb <- samples_tb %>%
        dplyr::mutate(time_series = purrr::map(time_series, compute_pca,
                                               pca_matrix = pca_rotation)) %>%
        sits::sits_select_bands(PC1, PC2)
}

# Filter a sits tibble.
filter_sits_tibble <- function(selected_bands, selected_labels, samples_tb){
    samples_tb %>%
        dplyr::filter(label %in% selected_labels) %>%
        select_bands(selected_bands) %>%
        ensurer::ensure_that(nrow(.) > 0, err_desc = "Samples are missing") %>%
        dplyr::select(label, time_series) %>%
        tidyr::unnest(time_series) %>%
        return()
}

# Compute the PCA given a time series tibble and a matrix of PCA rotations.
#
# @param .data      A tibble.
# @param pca_matrix A matrix.  The rotations matrix resulting from a PCA analysis using stats::prcomp
# @return           A tibble. .data plus the principal components.
compute_pca <- function(.data, pca_matrix) {
    observation_mt <- .data %>%
        dplyr::select(tidyselect::starts_with(rownames(pca_matrix))) %>%
        as.matrix()
    if (ncol(observation_mt) != nrow(pca_matrix)) {
        warning(sprintf("non-conformable arguments: %s",
                        paste(rownames(pca_matrix), collapse = " ")))
        return(NA)
    }
    .data %>%
        dplyr::bind_cols(tibble::as_tibble(observation_mt %*% pca_matrix)) %>%
        return()
}


# Pile images into a VRT.
#
# @param .data    A tibble.
# @param out_dir  A length-one charater. Path to a directory.
# @param cmd      A length-one charater. Pattern used for calling gdalbuildvrt i.e. cmd = "/usr/bin/gdalbuildvrt -separate %s %s"
# @param var_file A name of column in .data containing the paths to the files to pile.
# @return         A tibble. .data plus the principal components.
helper_pile_vrt <- function(.data, out_dir, cmd, var_file){
    .Deprecated("Use bash scripts instead.")
    var_file <- rlang::enquo(var_file)
    file_tb <- .data %>%
        ensurer::ensure_that(nrow(.) == 36, err_desc = "Missing images!") %>%
        dplyr::arrange(img_date)
    out_file <- file_tb %>%
        dplyr::slice(1) %>%
        dplyr::select(mission, level, orbit, tile, acquisition, band, resolution) %>%
        unlist() %>%
        paste(collapse = '_') %>%
        paste0(".vrt")
    out_file <- file.path(out_dir, out_file)
    file_vec <- unlist(dplyr::pull(file_tb, !!var_file))
    cmd <- sprintf(cmd,
                   out_file, paste(file_vec, collapse  = ' '))
    system(cmd)
    file_tb %>%
        dplyr::slice(1) %>%
        dplyr::select(-tidyselect::any_of("safe_path", "file_path", "acquisition", "processing")) %>%
        dplyr::mutate(brick_file = out_file) %>%
        return()
}


# Get a list of the reference and predicted labels for the given samples.
#
# @param samples     A sits tibble of sample points.
# @param raster_obj  A raster object.
# @param label_vec   A character. The labels in the given raster file.
# @return       A tibble.
get_ref_pred <- function(samples, raster_obj, label_vec){
    if(any(is.null(names(label_vec)))){
        warning("Setting names of vector of labels of the input raster.")
        label_vec <- label_vec %>%
            magrittr::set_names(1:length(.))
    }
    spdf <- sp::SpatialPointsDataFrame(coords = as.matrix(samples[,c("longitude",
                                                                     "latitude")]),
                                       data = as.data.frame(samples["label"],
                                                            stringsAsFactors = FALSE),
                                       proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>%
            ensurer::ensure_that(all(sort(label_vec) %in% sort(unique(.$label))),
                                 err_desc = sprintf("Label missmatch between the raster and the samples: %s %s",
                                                    paste(sort(label_vec), collapse = '_'),
                                                    paste(sort(unique(.$label)), collapse = '_')))
    raster_obj %>%
        raster::extract(y = spdf, sp = TRUE) %>%
        sf::st_as_sf() %>%
        sf::st_set_geometry(NULL) %>%
        tibble::tibble() %>%
        magrittr::set_colnames(c("reference", "predicted")) %>%
        dplyr::mutate(predicted = as.character(predicted)) %>%
        dplyr::mutate(predicted = dplyr::recode(predicted, !!!label_vec)) %>%
        ensurer::ensure_that(sum(is.na(.$predicted)) == 0,
                             err_desc = "missing labels!") %>%
        return()
}

# Build a tibble of metadata of the classification maps.
#
# @param in_dir A length-one character. Path to a directory with classification results.
# @return       A tibble.
get_results <- function(in_dir){
    # Check if the type of classfication is first images of the brick.
    get_class_type <- function(x){
        my_pattern = "first|last"
        if(!stringr::str_detect(x, my_pattern))
            return(NA_character_)
        stringr::str_subset(stringr::str_split(x, '/', simplify = TRUE),
                            pattern = my_pattern)
    }
    # Get n directories up in the given path.
    dir_up <- function(x, n){
        if(n < 1)
            return(x)
        dir_up(dirname(x), n - 1)
    }

    in_dir %>%
        list.files(pattern = "[.]tif", recursive = TRUE, full.names = TRUE) %>%
        tibble::enframe(name = NULL) %>%
        dplyr::rename(file_path = value) %>%
        dplyr::mutate(file_name = basename(file_path),
                      class_type = purrr::map_chr(dirname(file_path),
                                                  get_class_type),
                      file_dir   = ifelse(is.na(class_type),
                                          dirname(file_path),
                                          dir_up(dirname(file_path), 1)),
                      class_type = ifelse(is.na(class_type), "full", class_type),
                      used_method  = basename(dir_up(file_dir, 0)),
                      used_labels  = basename(dir_up(file_dir, 1)),
                      used_bands   = basename(dir_up(file_dir, 2)),
                      used_samples = basename(dir_up(file_dir, 3)),
                      brick_type   = basename(dir_up(file_dir, 4))) %>%
        dplyr::mutate(img_type = dplyr::case_when(stringr::str_detect(file_name, pattern = "postprocessing_") ~ "postprocessing",
                                                  stringr::str_detect(file_name, pattern = "_probs_class") ~ "classification",
                                                  stringr::str_detect(file_name, pattern = "_probs_")      ~ "probability")) %>%
        dplyr::select(-file_name, -file_dir) %>%
        return()
}


# Slice n of the first or last rows in the given tibble.
#
# @param x     A tibble.
# @param n     A numeric. The number of rows to slice.
# @param where A length-one charater. Indicate from where to take the rows; either the first or the last.
# @return      A tibble.
slice_n <- function(x, n = 2, where = "first"){
    if(where == "first"){
        return(dplyr::slice(x, 1:n))
    }else if (where == "last"){
        return(dplyr::slice(x, (nrow(x) - (n - 1)):nrow(x)))
    }
    stop("Unknown option!")
}


# Add coordinates as columns to an SF object.
#
# @param point_sf A sf object.
# @return         A sf object.
add_coords <- function(point_sf){
         xy <- point_sf %>%
             sf::st_coordinates() %>%
             magrittr::set_colnames(c("longitude", "latitude")) %>%
             tidyr::as_tibble()
         point_sf %>%
             dplyr::bind_cols(xy) %>%
             return()
}


# Apply postprocessing rules to join the partial classification of the brick (either using the first or last images) to the full classification of the brick.
#
# @param partial_class A length-one character. Path to a raster of a classification of the first or last images in a brick.
# @param full_class    A length-one character. Path to a raster of a classification of a full brick.
# @param rules         A length-one character. Rules to be applied to partial_class and full_class using gdal_calc syntax.
# @param partial       A length-one character. Indicate if either the first or last part of the brick is bein postprocesssed.
# @param out_dir       A length-one character. Path to a dir to store the results.
# @return              A path to the results.
postprocessing <- function(partial_class, full_class, rules,
                           partial = "first", out_dir = NULL){
    stopifnot(partial %in% c("first", "last"))
    if (is.null(out_dir)) {
        out_file <- full_class %>%
            dirname() %>%
            file.path(paste0("postprocessing_", partial, ".tif"))
    }else{
        if(!dir.exists(out_dir)){
            warning(sprintf("Creating dir %s", out_dir))
            dir.create(out_dir, recursive = TRUE)
        }
        out_file <- out_dir %>%
            file.path(paste0("postprocessing_", partial, ".tif"))
    }
    cmd <- sprintf("/usr/bin/gdal_calc.py -A %s -B %s --outfile=%s --calc='(%s).astype(int16)' --NoDataValue=-9999 --type='Int16' --creation-option='COMPRESS=LZW' --creation-option='BIGTIFF=YES'",
                   partial_class, full_class, out_file, rules)
    system(cmd)
    return(out_file)
}
