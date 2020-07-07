#!/bin/bash
shopt -s nullglob

# TODO: Transform to input variable.
out_dir=/disks/d3/brick_sentinel2

band="$1"
if [ "$band" == "" ] ; then
    echo "Expected 1 parameter. The name of a band."
    exit 1
fi
files=("${out_dir}"/vrt/*"${band}".vrt)
if [ ${#files[@]} -eq 0 ] ; then
    echo "ERROR: Band not found!"
    exit 1
fi
mkdir -p "${out_dir}"/tif 2> /dev/null
fmask_10m=("${out_dir}"/vrt/*Fmask4_10m.vrt)
if [ ${#fmask_10m[@]} -ne ${#files[@]} ]; then
    echo "ERROR: Number of masks don't match the number of files!"
    exit 1
fi
out_files=("${files[@]%_*}")

parameters=()
for i in "${!files[@]}" ; do
    parameters+=("${files[$i]} ${fmask_10m[$i]} ${out_dir}/tif/${out_files[$i]##*/}_masked_10m.tif")
done

calc="'(numpy.where(1.0 * B != 4.0, 1.0 * A, -9999.0)).astype(int16)'"
parallel --colsep ' ' /usr/bin/gdal_calc.py --quiet -A {1} -B {2} --outfile={3} --calc="${calc}" --NoDataValue=-9999.0 --type='Int16' --creation-option='COMPRESS=LZW' --creation-option='BIGTIFF=YES' --overwrite --quiet ::: "${parameters[@]}" 

for file in "${out_files[@]}" ; do
    /usr/bin/gdalbuildvrt -overwrite "${file}_masked_10m.vrt" "${out_dir}/tif/${file##*/}_masked_10m.tif"
done

exit 0
