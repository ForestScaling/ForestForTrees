#' Harvard Forest Sample Dataset
#'
#' A sample dataset derived from the Harvard Forest CTFS-ForestGEO Mapped Forest Plot
#' (HF253) and supplemented with canopy height measurements from the NEON Airborne
#' Observation Platform (AOP). This dataset provides tree-level information used in
#' the \pkg{ScalingFromSkypackage} for examples and testing.
#'
#' The dataset contains a subset of the full Harvard Forest plot (see Source), including
#' bounding box coordinates, canopy structural attributes, and derived estimates of
#' diameter at breast height (DBH).
#'
#' @format A data frame with 5903 observations and 23 variables:
#' \describe{
#'   \item{xmin, ymin, xmax, ymax}{Bounding box coordinates for each tree crown (numeric).}
#'   \item{label}{Object classification label (character, e.g. "Tree").}
#'   \item{score}{Confidence score for the tree detection (numeric).}
#'   \item{image_path}{Path to the remote sensing image tile (character).}
#'   \item{left, right, top, bottom}{Image coordinates (numeric).}
#'   \item{Area, Perimeter}{Polygon area and perimeter of the crown segment (numeric).}
#'   \item{OBJECTID, Join_Count, TARGET_FID, Id}{Identifiers and join indices (numeric/integer).}
#'   \item{IDhectbest}{Hectare identifier (numeric).}
#'   \item{Shape_Leng, Shape_Area}{Shape geometry attributes (numeric).}
#'   \item{Max_Height}{Maximum canopy height from NEON AOP canopy height model (numeric, meters).}
#'   \item{Diameter}{Estimated canopy diameter (numeric, meters).}
#'   \item{dbh}{Estimated diameter at breast height (DBH, numeric, centimeters).}
#' }
#'
#' @details
#' Canopy height values (\code{Max_Height}) were extracted from the NEON Canopy Height
#' Model (CHM) using zonal statistics applied to Harvard Forest crown segments.
#'
#' Canopy diameter (\code{Diameter}) was derived from crown geometry following:
#' \deqn{Diameter = 0.5 \times \sqrt{Perimeter^2 - 8 \times Area}}
#'
#' Diameter at breast height (\code{dbh}) was estimated from canopy height and canopy
#' diameter using the allometric equation implemented in \pkg{itcSegment}:
#' \deqn{dbh = itcSegment::dbh(H = Max\_Height, CA = Diameter)}
#'
#' IDhectbest is not in the original ForestGEO data, and is instead the result of a
#' fishnet shapefile with 1 hectare squares overlaid on top of the original plot.
#' The original plot is much larger than 1 hectare.
#'
#' This dataset is intended for demonstration and testing purposes within the
#' \pkg{ScalingFromSky} package. It includes only a subset of the original Harvard Forest
#' dataset and should not be used as a replacement for the full CTFS-ForestGEO plot data.
#'
#' @source
#' Harvard Forest CTFS-ForestGEO Mapped Forest Plot (HF253):
#' \url{https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=hf253}
#' NEON Airborne Observation Platform (AOP): \url{https://www.neonscience.org/data-collection/airborne-remote-sensing}
#'
#' @examples
#' data(harv_data_sample)
#' head(harv_data_sample)
#'
#' plot(harv_data_sample$Max_Height, harv_data_sample$dbh,
#'      xlab = "Canopy Height (m)", ylab = "DBH (cm)",
#'      main = "Height vs DBH in Harvard Forest sample")
#'
#' @docType data
#' @name harv_data_sample
NULL

#' NEON Site Alpha Parameter Estimations (Benchmark)
#'
#' A dataset containing estimated alpha parameters and their associated standard
#' deviations for 25 National Ecological Observatory Network (NEON) sites.
#'
#' These estimations were calculated using the full available dataset for each
#' site. This dataset serves as a benchmark or "ground truth" for comparison
#' against the interpolated alpha values generated in the package vignettes.
#'
#' @format A data frame with 25 rows and 3 variables:
#' \describe{
#'   \item{siteID}{A 4-letter code identifying the NEON site (e.g., "BART", "HARV").}
#'   \item{alpha}{The estimated alpha parameter calculated using the full site dataset.}
#'   \item{alpha_sd}{The standard deviation (uncertainty) of the alpha estimation.}
#' }
#'
#' @details
#' These values represent the "true" measured state of the alpha parameter
#' when all available data is utilized. In the package vignettes, these are
#' used to validate the accuracy of interpolation methods.
#'
#' @source Data provided by the National Ecological Observatory Network (NEON).
#' @examples
#' # Access the data
#' data(neon_truemeasuredalpha)
#' @docType data
#' @name "neon_truemeasuredalpha"
NULL

#' Remote Sensing Tree Canopy Polygons (NEON)
#'
#' A Simple Feature (sf) collection containing 21,956 individual tree canopy
#' polygons derived from remote sensing data across various NEON sites.
#'
#' @format A simple feature collection (sf) with 21,956 features and 15 fields:
#' \describe{
#'   \item{FID}{Internal feature identifier.}
#'   \item{siteID}{4-character NEON site code (e.g., "HARV", "BART").}
#'   \item{plotID}{Specific NEON plot identifier (e.g., "HARV_022").}
#'   \item{area_sq_m}{Area of the canopy polygon in square meters.}
#'   \item{left, bottom, right, top}{Bounding box coordinates for the individual polygon.}
#'   \item{score}{Quality score for the specific segmentation method used to delineate the canopy.}
#'   \item{perim_m}{Perimeter of the canopy polygon in meters.}
#'   \item{utmZone}{The UTM projection zone where the measurement was originally taken.}
#'   \item{Height}{Height of the tree in meters.}
#'   \item{Diameter}{Crown diameter of the tree in meters.}
#'   \item{dbh}{Diameter at breast height of the tree in centimeters.}
#'   \item{geometry}{The sfc_POLYGON column containing the spatial boundaries (WGS 84).}
#' }
#'
#' @details
#' This dataset provides high-resolution individual tree crown (ITC) data.
#' The polygons are intended to be used in conjunction with site-level
#' alpha estimations for spatial analysis and scaling exercises in the
#' package vignettes.
#'
#' @source Data provided by the National Ecological Observatory Network (NEON).
#' @examples
#' library(sf)
#' data(neon_data)
#' @docType data
#' @name "neon_data"
NULL

#' Random Forest Priors for Alpha (FIA-derived)
#'
#' A dataset containing prior estimations for the alpha parameter at NEON sites.
#' These priors were generated using a Random Forest model trained on Forest
#' Inventory and Analysis (FIA) data.
#'
#' @format A data frame with 27 rows (one per NEON site) and 3 variables:
#' \describe{
#'   \item{siteID}{4-character NEON site code (e.g., "BART", "HARV").}
#'   \item{prior_mean}{The mean alpha value predicted by the FIA-trained Random Forest model.}
#'   \item{prior_sd}{The standard deviation (uncertainty) associated with the prior prediction.}
#' }
#'
#' @details
#' These priors provide a baseline expectation for alpha based on regional forest
#' characteristics captured by the FIA program. In a Bayesian framework or
#' comparative analysis, these can be used as the starting point before
#' incorporating site-specific remote sensing data (`neon_data`) or
#' validation against ground truth (`neon_truemeasuredalpha`).
#'
#' @source Model derived from USDA Forest Service Forest Inventory and Analysis (FIA) Program data.
#' @examples
#' data(alpha_priors)
#' @docType data
#' @name "alpha_priors"
NULL

#' Site-Level Leaf Area Index (LAI)
#'
#' A dataset containing the Leaf Area Index (LAI) values for 36 research sites,
#' typically associated with the National Ecological Observatory Network (NEON).
#' LAI is a dimensionless measure of the one-sided green leaf area per unit ground
#' surface area.
#'
#' @format A data frame with 36 rows and 2 variables:
#' \describe{
#'   \item{site}{Character, 4-letter unique site identifier (e.g., "ABBY", "BART").}
#'   \item{Leaf_area_index}{Numeric, the calculated Leaf Area Index (m^2/m^2).
#'   Values typically range from 0 (bare ground) to >5 (dense forest).}
#' }
#'
#' @details
#' The Leaf Area Index is a critical variable in ecological modeling as it
#' influences the amount of light intercepted by the canopy and the rate of
#' transpiration.
#'
#' data(Leaf_area_index)
#' #' @docType data
#' @name "Leaf_area_index"
