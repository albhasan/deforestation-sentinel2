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
│   ├── raster
│   │   ├── sentinel_L1C
│   │   └── sentinel_L2A
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

## Input images

The Sentinel-2A images below are expected in the directory `./data/raster/sentinel_L1C`. The results of processing them with Sen2Cor are expected in the directory `./data/raster/sentinel_L2A`

```bash
S2A_MSIL1C_20180812T143751_N0206_R096_T20LKP_20180812T182110.SAFE
S2A_MSIL1C_20180822T143751_N0206_R096_T20LKP_20180822T192224.SAFE
S2A_MSIL1C_20180901T143741_N0206_R096_T20LKP_20180901T182920.SAFE
S2A_MSIL1C_20180911T143741_N0206_R096_T20LKP_20180911T181809.SAFE
S2A_MSIL1C_20180921T143741_N0206_R096_T20LKP_20180921T181736.SAFE
S2A_MSIL1C_20181001T143741_N0206_R096_T20LKP_20181001T163418.SAFE
S2A_MSIL1C_20181011T143751_N0206_R096_T20LKP_20181011T181532.SAFE
S2A_MSIL1C_20181021T143751_N0206_R096_T20LKP_20181021T180732.SAFE
S2A_MSIL1C_20181031T143751_N0206_R096_T20LKP_20181103T085637.SAFE
S2A_MSIL1C_20181110T143751_N0207_R096_T20LKP_20181110T181533.SAFE
S2A_MSIL1C_20181120T143751_N0207_R096_T20LKP_20181120T162623.SAFE
S2A_MSIL1C_20181130T143741_N0207_R096_T20LKP_20181130T162039.SAFE
S2A_MSIL1C_20181210T143741_N0207_R096_T20LKP_20181210T165359.SAFE
S2A_MSIL1C_20181220T143741_N0207_R096_T20LKP_20181220T161934.SAFE
S2A_MSIL1C_20181230T143751_N0207_R096_T20LKP_20181230T162049.SAFE
S2A_MSIL1C_20190109T143751_N0207_R096_T20LKP_20190109T162038.SAFE
S2A_MSIL1C_20190119T143751_N0207_R096_T20LKP_20190119T162012.SAFE
S2A_MSIL1C_20190129T143751_N0207_R096_T20LKP_20190225T132350.SAFE
S2A_MSIL1C_20190208T143751_N0207_R096_T20LKP_20190208T180253.SAFE
S2A_MSIL1C_20190218T143751_N0207_R096_T20LKP_20190218T175945.SAFE
S2A_MSIL1C_20190228T143751_N0207_R096_T20LKP_20190228T180016.SAFE
S2A_MSIL1C_20190310T143751_N0207_R096_T20LKP_20190310T180430.SAFE
S2A_MSIL1C_20190320T143751_N0207_R096_T20LKP_20190320T194256.SAFE
S2A_MSIL1C_20190330T143751_N0207_R096_T20LKP_20190330T180249.SAFE
S2A_MSIL1C_20190409T143751_N0207_R096_T20LKP_20190409T180242.SAFE
S2A_MSIL1C_20190419T143801_N0207_R096_T20LKP_20190419T180222.SAFE
S2A_MSIL1C_20190429T143751_N0207_R096_T20LKP_20190429T192231.SAFE
S2A_MSIL1C_20190509T143751_N0207_R096_T20LKP_20190509T180139.SAFE
S2A_MSIL1C_20190519T143751_N0207_R096_T20LKP_20190519T162104.SAFE
S2A_MSIL1C_20190529T143751_N0207_R096_T20LKP_20190529T192149.SAFE
S2A_MSIL1C_20190608T143751_N0207_R096_T20LKP_20190608T161936.SAFE
S2A_MSIL1C_20190618T143751_N0207_R096_T20LKP_20190618T192212.SAFE
S2A_MSIL1C_20190628T143751_N0207_R096_T20LKP_20190628T192236.SAFE
S2A_MSIL1C_20190708T143751_N0208_R096_T20LKP_20190708T162020.SAFE
S2A_MSIL1C_20190718T143751_N0208_R096_T20LKP_20190718T161930.SAFE
S2A_MSIL1C_20190728T143751_N0208_R096_T20LKP_20190728T180109.SAFE
```
