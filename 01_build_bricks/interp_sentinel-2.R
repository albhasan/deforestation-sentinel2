# Adapted from an script by Rolf Simoes - https://github.com/rolfsimoes

library(raster)
library(snow)

job_spline_fun <- function(job, data_type = "INT2S",
                           no_data = -9999, options = "COMPRESS=LZW", ...) {

  function(job) {
    my_spline <- function(v) {

      t(apply(v, 1,
              function(x) {
                na <- is.na(x)
                i <- tryCatch({
                  spline(x, n = length(x), ...)
                }, error = function(e) list(y = rep(0, length(x))))
                x[na] <- i$y[na]
                return(x)
              }))
    }

    b <- raster::brick(job$f)
    b <- raster::crop(b, raster::extent(b, job$row, job$row + job$nrows - 1, 1, raster::ncol(b)))
    out <- raster::calc(b, my_spline, datatype = data_type,
                        NAflag = no_data, options = options, format = "GTiff",
                        filename = tempfile(pattern = "job_spline_fun_", fileext = ".tif"))

    return(raster::filename(out))
  }
}

job_approx_fun <- function(out_dir = tempdir(), data_type = "INT2S",
                           no_data = -9999, options = "COMPRESS=LZW", ...) {

  function(job) {

    my_approx <- function(v) {

      t(apply(v, 1,
              function(x) {
                na <- is.na(x)
                i <- tryCatch({
                  approx(x, n = length(x), ...)
                }, error = function(e) list(y = rep(0, length(x))))
                x[na] <- i$y[na]
                return(x)
              }))
    }

    #file_name <- tempfile(pattern = "job_approx_fun_", tmpdir = out_dir, fileext = ".tif")
    b <- raster::brick(job$f)
    b <- raster::crop(b, raster::extent(b, job$row, job$row + job$nrows - 1, 1, raster::ncol(b)))
    out <- raster::calc(b, my_approx, datatype = data_type,
                        NAflag = no_data, options = options, format = "GTiff",
                        filename = tempfile(pattern = "job_approx_fun_", fileext = ".tif"))
    return(raster::filename(out))
  }
}

jobs <- function(f, fun, out_file, data_type = "INT2S",
                 no_data = -9999, options = "COMPRESS=LZW", cl = NULL) {

  b <- raster::brick(f)
  bs <- raster::blockSize(b)
  bs$n <- NULL
  jobs <- unname(do.call(mapply, args = c(list(FUN = list, SIMPLIFY = FALSE, f = f), bs)))

  if (is.null(cl))
    cl <- 5

  cl <- snow::makeSOCKcluster(cl)
  on.exit({snow::stopCluster(cl)})

  cat("Processing...")
  l <- snow::clusterApplyLB(cl, jobs, fun)

  cat("Merging result...")
  do.call(raster::merge,
          args = c(list(filename = out_file, datatype = data_type, NAflag = no_data,
                        options = options, format = "GTiff"), lapply(l, raster::brick)))
  cat(out_file)
  raster::removeTmpFiles()

  return(out_file)
}

args = commandArgs(trailingOnly = TRUE)
if ((length(args) != 3) | (!args[1] %in% c("approx", "spline")))
  stop("Usage Rscript --vanilla interp_raster.R [approx|spline] <input> <output>")

rasterOptions(format = "GTiff", maxmemory = 40*1024*1024*1024, memfrac = 0.8, progress="text", chunksize=500*1024*1024)

if (args[1] == "approx")
  jobs(args[2], fun = job_approx_fun(), out_file = args[3], cl = 10, no_data = 0)

if (args[1] == "spline")
  jobs(args[2], fun = job_spline_fun(), out_file = args[3], cl = 10, no_data = 0)

