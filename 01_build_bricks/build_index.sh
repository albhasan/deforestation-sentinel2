#!/bin/bash

# evi_10m ndmi_10m ndvi_10m savi_10m

shopt -s nullglob

# TODO: Transform to input variable.
out_dir="./brick_sentinel2"
index="$1"

mkdir -p "${out_dir}"/tif 2> /dev/null

B02_10m=("${out_dir}"/vrt/*B02_10m.vrt)
#B03_10m=("${out_dir}"/vrt/*B03_10m.vrt)
B04_10m=("${out_dir}"/vrt/*B04_10m.vrt)
B08_10m=("${out_dir}"/vrt/*B08_10m.vrt)
B11_10m=("${out_dir}"/vrt/*B11_10m.vrt)
#B12_10m=("${out_dir}"/vrt/*B12_10m.vrt)
#B8A_10m=("${out_dir}"/vrt/*B8A_10m.vrt)

if [ ${#B02_10m[@]} -ne ${#B04_10m[@]} ] || \
   [ ${#B02_10m[@]} -ne ${#B08_10m[@]} ] || \
   [ ${#B02_10m[@]} -ne ${#B11_10m[@]} ]; then
    echo "ERROR: Number of file bands don't match!"
    exit 1
fi

out_files=("${B02_10m[@]%B02_10m*}")

parameters=()
for i in "${!B02_10m[@]}" ; do
    #NOTE: 1 = B02_10m, 2 = B04_10m, 3 = B08_10m, 4 = B11_10m
    parameters+=("${B02_10m[$i]} ${B04_10m[$i]} ${B08_10m[$i]} ${B11_10m[$i]} ${out_dir}/tif/${out_files[$i]##*/}${index}.tif")
done

if [ "${index}" == "evi_10m" ]; then

    #NOTE: A = B08_10m, B = B04_10m, C = B02_10m
    calc="'(10000 * 2.5 * (1.0 * A - B) / (A + 6.0 * B - 7.5 * C + 10000)).astype(numpy.int16)'"
    parallel  --colsep ' ' /usr/bin/gdal_calc.py --quiet -A {3} -B {2} -C {1}  --outfile={5} --calc="${calc}" --NoDataValue=-9999 --type='Int16' --creation-option='COMPRESS=LZW' --creation-option='BIGTIFF=YES' --overwrite --quiet ::: "${parameters[@]}" 

elif [ "${index}" == "ndmi_10m" ]; then

    #NOTE: A = B08_10m, B = B11_10m
    calc="'((1.0 * A - B)/(1.0 * A + B) * 10000).astype(numpy.int16)'"
    parallel  --colsep ' ' /usr/bin/gdal_calc.py --quiet -A {3} -B {4} --outfile={5} --calc="${calc}" --NoDataValue=-9999 --type='Int16' --creation-option='COMPRESS=LZW' --creation-option='BIGTIFF=YES' --overwrite --quiet ::: "${parameters[@]}" 

elif [ "${index}" == "ndvi_10m" ]; then

    #NOTE: A = B08_10m, B = B04_10m
    calc="'((1.0 * A - B)/(1.0 * A + B) * 10000).astype(numpy.int16)'"
    parallel  --colsep ' ' /usr/bin/gdal_calc.py --quiet -A {3} -B {2} --outfile={5} --calc="${calc}" --NoDataValue=-9999 --type='Int16' --creation-option='COMPRESS=LZW' --creation-option='BIGTIFF=YES' --overwrite --quiet ::: "${parameters[@]}" 

elif [ "${index}" == "savi_10m" ]; then

    #NOTE: A = B08_10m, B = B04_10m
    calc="'((1.0 * A - B)/(1.0 * A + B + 4280.0) * (10000 + 4280)).astype(numpy.int16)'"
    parallel  --colsep ' ' /usr/bin/gdal_calc.py --quiet -A {3} -B {2} --outfile={5} --calc="${calc}" --NoDataValue=-9999 --type='Int16' --creation-option='COMPRESS=LZW' --creation-option='BIGTIFF=YES' --overwrite --quiet ::: "${parameters[@]}" 

else 
    echo "Unknown index!"
    exit 1
fi

for file in "${out_files[@]}" ; do
    /usr/bin/gdalbuildvrt -overwrite "${file}${index}.vrt" "${out_dir}/tif/${file##*/}${index}.tif"
done

exit 0
