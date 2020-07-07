#!/bin/bash

# Path to base directory with scripts.
script_dir="/home/alber/Documents/ghProjects/deforestation-sentinel2"

# Path to Sentinel images
sentinel_L1C="/home/alber/Documents/data/experiments/prodes_reproduction/papers/deforestation/data/raster/sentinel_L1C"
sentinel_L2A="/home/alber/Documents/data/experiments/prodes_reproduction/papers/deforestation/data/raster/sentinel_L2A"

# Path to brick files.
brick_dir="/disks/d3/brick_sentinel2"

#---- Build bricks ----

# Build VRTs from the bands.
"${script_dir}"/01_build_bricks/./build_vrt_base.sh bilinear "${sentinel_L2A}"/*/GRANULE/*/IMG_DATA/*/*.jp2

# Build VRTs from the FMASKs.
"${script_dir}"/01_build_bricks/./build_vrt_base.sh nearest "${sentinel_L1C}"/S2A_MSIL1C_201*/GRANULE/*/FMASK_DATA/*.tif

# NOTE: build_index.sh already uses GNU parallel to run.
parallel -j 1 "${script_dir}"/01_build_bricks/./build_index.sh ::: evi_10m ndmi_10m ndvi_10m savi_10m

# NOTE: build_index.sh already uses GNU parallel to run.
parallel -j 1 "${script_dir}"/01_build_bricks/./mask_images.sh ::: B02_10m B03_10m B04_10m B08_10m B11_10m B12_10m B8A_10m evi_10m ndmi_10m ndvi_10m savi_10m

parallel -j 8 "${script_dir}"/01_build_bricks/./build_vrt_brick.sh ::: B02_10m B03_10m B04_10m B08_10m B11_10m B12_10m B8A_10m evi_10m ndmi_10m ndvi_10m savi_10m

# NOTE: Build VRTs for masked bands and indeces
parallel -j 8 "${script_dir}"/01_build_bricks/./build_vrt_brick.sh ::: B02_masked_10m B03_masked_10m B04_masked_10m B08_masked_10m B11_masked_10m B12_masked_10m B8A_masked_10m evi_masked_10m ndmi_masked_10m ndvi_masked_10m savi_masked_10m

parallel -j 8 "${script_dir}"/01_build_bricks/./build_tif_brick.sh ::: B02_10m B03_10m B04_10m B08_10m B11_10m B12_10m B8A_10m evi_10m ndmi_10m ndvi_10m savi_10m

# NOTE: Build TIFs for masked bands and indeces
parallel -j 8 "${script_dir}"/01_build_bricks/./build_tif_brick.sh ::: B02_masked_10m B03_masked_10m B04_masked_10m B08_masked_10m B11_masked_10m B12_masked_10m B8A_masked_10m evi_masked_10m ndmi_masked_10m ndvi_masked_10m savi_masked_10m

#---- Interpolate masked bricks ----

Rscript "${script_dir}"/01_build_bricks/interp_sentinel-2.R approx "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_B02_masked_10m.tif  "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_B02_approx_10m.tif 
Rscript "${script_dir}"/01_build_bricks/interp_sentinel-2.R approx "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_B03_masked_10m.tif  "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_B03_approx_10m.tif
Rscript "${script_dir}"/01_build_bricks/interp_sentinel-2.R approx "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_B04_masked_10m.tif  "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_B04_approx_10m.tif
Rscript "${script_dir}"/01_build_bricks/interp_sentinel-2.R approx "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_B08_masked_10m.tif  "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_B08_approx_10m.tif
Rscript "${script_dir}"/01_build_bricks/interp_sentinel-2.R approx "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_B11_masked_10m.tif  "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_B11_approx_10m.tif
Rscript "${script_dir}"/01_build_bricks/interp_sentinel-2.R approx "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_B12_masked_10m.tif  "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_B12_approx_10m.tif
Rscript "${script_dir}"/01_build_bricks/interp_sentinel-2.R approx "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_B8A_masked_10m.tif  "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_B8A_approx_10m.tif
Rscript "${script_dir}"/01_build_bricks/interp_sentinel-2.R approx "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_evi_masked_10m.tif  "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_evi_approx_10m.tif
Rscript "${script_dir}"/01_build_bricks/interp_sentinel-2.R approx "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_ndmi_masked_10m.tif "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_ndmi_approx_10m.tif
Rscript "${script_dir}"/01_build_bricks/interp_sentinel-2.R approx "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_ndvi_masked_10m.tif "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_ndvi_approx_10m.tif
Rscript "${script_dir}"/01_build_bricks/interp_sentinel-2.R approx "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_savi_masked_10m.tif "${brick_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_savi_approx_10m.tif

#---- Install SITS ----

"${script_dir}"/other/install_sits.R

#---- Get time series ----

"${script_dir}"/02_get_time_series/get_time_series.R "${script_dir}"/data/validation/samples_all_bands.csv approx "${script_dir}"/data/validation/samples_A_approx.rds
"${script_dir}"/02_get_time_series/get_time_series.R "${script_dir}"/data/validation/samples_indices.csv   approx "${script_dir}"/data/validation/samples_B_approx.rds
"${script_dir}"/02_get_time_series/get_time_series.R "${script_dir}"/data/validation/samples_all_bands.csv raw    "${script_dir}"/data/validation/samples_A_raw.rds 
"${script_dir}"/02_get_time_series/get_time_series.R "${script_dir}"/data/validation/samples_indices.csv   raw    "${script_dir}"/data/validation/samples_B_raw.rds

"${script_dir}"/02_get_time_series/create_samples_3_labels.R

#---- Compute K-Folds ----

"${script_dir}"/03_kfolds/k-folds_analysis.R "${script_dir}"/data/validation/samples_B_approx_3l.rds "${script_dir}"/plot/kfold_approx
"${script_dir}"/03_kfolds/k-folds_analysis.R "${script_dir}"/data/validation/samples_B_raw_3l.rds    "${script_dir}"/plot/kfold_raw

#---- Classify bricks ----

samples_B_approx_3l="${script_dir}"/data/validation/samples_B_approx_3l.rds
three_labels="Deforestatio,Forest,NonForest"
bands="blue,bnir,green,nnir,red,swir1,swir2"
indices="evi,ndmi,ndvi"
version="009"
out_base_dir="${script_dir}"/results9

"${script_dir}"/04_classify/classify_bricks.R approx "${brick_dir}" "${samples_B_approx_3l}" "${three_labels}" "${bands}"   "${version}" "${out_base_dir}"
"${script_dir}"/04_classify/classify_bricks.R approx "${brick_dir}" "${samples_B_approx_3l}" "${three_labels}" "${indices}" "${version}" "${out_base_dir}"

#---- Post-processing ----

#---- Validation ----

"${script_dir}"/06_validation/validate_results.R 

exit 0
