# ForestForTrees

<!-- badges: start -->

[![R-CMD-check](https://github.com/ForestScaling/ForestForTrees/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ForestScaling/ForestForTrees/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

`ForestForTrees` is an R package that operationalizes a Bayesian scaling framework for estimating full forest size-abundance distributions from remotely sensed canopy data. Starting from crown segmentation outputs (DBH estimates and LAI), the package reconstructs the complete 10–50 cm DBH distribution — including understory trees hidden by canopy occlusion — using a truncated Pareto model fit via Stan/MCMC. The package is sensor- and algorithm-agnostic: it works with any upstream segmentation method (LiDAR or RGB, DeepForest, watershed, etc.) and requires only a numeric vector of estimated DBHs and a site-level LAI value.

The package implements the methods published in Eichenwald et al. (2025), *Global Ecology and Biogeography*.

## Key Features

- **Generalized pipeline**: Automated MCMC sampling, uncertainty propagation, and observational window detection in a sequential, reproducible workflow compatible with any upstream remote sensing pipeline
- **Gridded mapping**: `map_grid_estimates()` joins plot-level estimates of α and N_tot to spatial polygon layers, enabling landscape-scale raster mapping of forest structural parameters without additional scripting
- **Cross-biome performance**: Validated across 23 NEON sites spanning subtropical to alpine forest types, recovering both α (R² = 0.54) and total abundance (R² = 0.61) with no systematic bias and without site-specific reconfiguration
- **Flexible Bayesian modeling**: Pareto model with LAI and breakpoint corrections for canopy occlusion
- **Uncertainty propagation**: Full posterior uncertainty from α propagates through to N_tot estimates

## Included Data

- **`harv_data_sample`**: A curated sample dataset derived from the Harvard Forest CTFS-ForestGEO plot (HF253) and NEON Airborne Observation Platform (AOP). Includes tree crown coordinates, canopy structural attributes, derived DBH estimates, and hectare plot IDs. Intended for demonstration, testing, and the Harvard Forest workflow vignette.

- **`harvardshapefile`**: A sample sf object of delineated tree crown polygons covering a larger area of the Harvard Forest ForestGEO site. UTM-referenced, with attributes including crown perimeter, area, maximum canopy height, and hectare plot ID (`IDhectbest`). Used in the gridded mapping demonstration.

- **`neon_data`**: An sf collection of 21,956 individual tree canopy polygons derived from remote sensing data across 23 NEON sites. Includes DBH estimates, crown geometry, and site metadata. Used in the cross-biome benchmark vignette.

- **`neon_truemeasuredalpha`**: Field-measured Pareto α values for NEON sites, used as ground truth for package validation.

- **`alpha_priors`**: Random forest-derived prior means and standard deviations for α at each NEON site, trained on FIA data.

- **`Leaf_area_index`**: Site-level LAI values for NEON sites derived from MODIS products.

## Vignettes

The package includes two fully reproducible vignettes:

1. **Estimating Forest Size-Abundance Parameters** (`ForestForTrees: Estimating Forest Size-Abundance Parameters`): A step-by-step workflow demonstration using Harvard Forest data. Covers the full pipeline from DBH input to α and N_tot estimation, multi-plot looping, parallelization, and gridded spatial mapping using `map_grid_estimates()`.

2. **Cross-Biome Benchmark and Performance Evaluation** (`ForestForTrees: Cross-Biome Benchmark and Performance Evaluation`): Applies the pipeline across 23 NEON sites to validate parameter recovery and report runtime benchmarks. Demonstrates that the automated workflow converges reliably across diverse forest types without manual reconfiguration.

## Purpose

This package makes previously published methods for estimating forest size-abundance distributions accessible to researchers and practitioners. It is particularly useful for:

- Forest ecologists analyzing structural changes across plots or regions
- Researchers integrating remote sensing data with field measurements
- Large-scale forest monitoring and carbon accounting workflows
- Teaching and demonstration of Bayesian forest structural modeling

## Installation

```r
# Install from GitHub
devtools::install_github("ForestScaling/ForestForTrees")
```

## Basic Usage

```r
library(ForestForTrees)
library(dplyr)

data("harv_data_sample")

# Single plot
df <- harv_data_sample %>% filter(IDhectbest == 1, !is.na(dbh)) %>% select(dbh)

kde_out   <- potential_break(df)
trunc_out <- truncate_filter(kde_out)
alpha_fit <- fit_alpha_model(trunc_out, LAI = 5.4, prior_mean = 1.4, prior_sd = 0.3)
ntot_fit  <- estimate_total_trees(alpha_fit)
```

## License

MIT License
