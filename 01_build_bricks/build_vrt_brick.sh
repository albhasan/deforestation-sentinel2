#!/bin/bash
shopt -s nullglob

# TODO: Transform to input variable.
out_dir="./brick_sentinel2"

band="$1"
if [ "$band" == "" ] ; then
    echo "Expected 1 paramter. The name of a band."
    exit 1
fi

files=("${out_dir}"/vrt/*"${band}".vrt)
if [ ${#files[@]} -eq 0 ] ; then
    echo "Band not found!"
    exit 1
fi

vrt_files=(${out_dir}/vrt/*"${band}".vrt)
/usr/bin/gdalbuildvrt -separate -overwrite -resolution user -tr 10 10 -r bilinear -q "${out_dir}"/S2A_MSIL2A_R096_T20LKP_20180812T143751_"${band%_*}"_10m.vrt "${vrt_files[@]}"
