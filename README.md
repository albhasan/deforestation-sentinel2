[![DOI](https://zenodo.org/badge/277586366.svg)](https://zenodo.org/badge/latestdoi/277586366)

# deforestation-sentinel2

Scripts for creating the classification maps described in the paper entitled "Detecting Tropical Deforestation using Time Series of Sentinel-2A Images" authored by Alber Sanchez, Michelle C. A. Picoli, Rolf Simoes, Gilberto Camara, Pedro Andrade, and Karine Reis.

## Instalation

The scripts can be installed using the following command:

```bash
git clone https://github.com/albhasan/deforestation-sentinel2
```

## Contents

```bash
.
├── 01_build_bricks
│   ├── build_index.sh                     Compute vegetation indices. 
│   ├── build_tif_brick.sh                 Convert inputs to GeoTif.
│   ├── build_vrt_base.sh                  Create GDAL virtual rasters of 10x10 meter resolution.
│   ├── build_vrt_brick.sh                 Ensemble inputs into image bricks.
│   ├── interp_sentinel-2.R                Interpolate NA pixels along the time dimension.
│   └── mask_images.sh                     Replace clouds with NA pixels.
├── 02_get_time_series
│   ├── create_samples_3_labels.R          Recode sample labels, and save them as RDS.
│   └── get_time_series.R                  Retrieve time series for sample points.
├── 03_kfolds
│   └── k-folds_analysis.R                 Test combination of attributes using K-fold.
├── 04_classify
│   └── classify_bricks.R                  Classify image bricks using Random Forest.
├── 05_post-processing
├── 06_validation
│   └── validate_results.R                 Validate the classifications.
├── brick_sentinel2
│   ├── tif
│   └── vrt
├── data
│   └── validation
│       ├── samples_A_approx.rds           Intermediary file.
│       ├── samples_all_bands.csv          Sample points (WGS84).
│       ├── samples_A_raw.rds              Intermediary file.
│       ├── samples_B_approx_3l.rds        Intermediary file. 
│       ├── samples_B_approx.rds           Intermediary file.
│       ├── samples_B_raw_3l.rds           Intermediary file.
│       ├── samples_B_raw.rds              Intermediary file.
│       └── samples_indices.csv            Sample points (WGS84).
├── LICENSE
├── other
│   ├── config.yml                         Configuration file of the sits package.
│   ├── install_sits.R                     Installation script of the sits package.
│   └── util.R                             Utilitary functions.
├── pangea
├── plot
│   ├── kfold_approx
│   └── kfold_raw
├── README.md                              This file.
├── results9
├── run_experiment.sh                      Main script for reproducibility.
└── tmp
```
