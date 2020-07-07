#!/bin/bash
shopt -s nullglob

# TODO: Transform to input variable.
out_dir=/disks/d3/brick_sentinel2/vrt

test=$(gdalbuildvrt --help | grep -E "(\{|,)$1(\}|,)")
if [ "${test}" == "" ]; then
    echo "ERROR: Invalid interpolation!"
    exit 1
fi

interpolation="$1"
shift
files=("$@")
if [ "${#files[@]}" -eq 0 ]; then
    echo "ERROR: Missing images!"
    exit 1
fi

mkdir -p "${out_dir}" 2> /dev/null

for file in "${files[@]}" ; do
    out_file="${file##*/}"
    out_file="${out_file%.*}"
    if [ -f "${out_dir}/${out_file%_*m}"_10m.vrt ]; then  
        continue
    fi
    /usr/bin/gdalbuildvrt -resolution user -tr 10 10 -r "${interpolation}" -overwrite -q "${out_dir}/${out_file%_*m}"_10m.vrt "${file}"
done

exit 0
