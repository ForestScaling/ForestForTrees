#' Crown segmentation polygons and canopy height data for Harvard Forest (ForestGEO)
#'
#' A spatial dataset containing tree crown segmentation polygons from a subset of the
#' Harvard Forest ForestGEO plot. Each polygon represents a segmented tree crown derived
#' from high-resolution remote sensing imagery, with corresponding canopy height values
#' extracted from a lidar-based Canopy Height Model (CHM).
#'
#' This dataset provides an example of the type of spatial input data used in
#' \pkg{ScalingFromSkypackage} for linking image-based crown detections with structural
#' and demographic estimates.
#'
#' @format A simple feature (\code{sf}) object with 6 features and 12 attributes:
#' \describe{
#'   \item{score}{Prediction confidence score (unitless, 0–1).}
#'   \item{image_path}{File name of the original remote sensing image used for crown detection.}
#'   \item{left, right, top, bottom}{UTM coordinates of the polygon extent.}
#'   \item{Area, Perimeter}{Polygon area and perimeter (in square meters and meters).}
#'   \item{IDhectbest}{Plot segment or identifier within the ForestGEO study area.}
#'   \item{Shape_Leng, Shape_Area}{Geometric length and area attributes from the original shapefile.}
#'   \item{Max_Height}{Maximum canopy height (in meters) extracted from the corresponding CHM raster.}
#'   \item{geometry}{Polygon geometry column in UTM Zone 18N (EPSG: 32618).}
#' }
#'
#' @details
#' The polygons in \code{harvardshapefile} represent crown segments generated from
#' high-resolution RGB imagery and spatially aligned with the Harvard Forest ForestGEO plot.
#' Each segment is matched with a lidar-derived canopy height estimate, allowing
#' for direct comparison between image-based crown structure and field-based
#' forest inventory measurements.
#'
#' This dataset is intended for demonstration and testing purposes within
#' \pkg{ScalingFromSkypackage}.
#'
#' @source
#' Derived from Harvard Forest ForestGEO data and NEON lidar products
#' (Canopy Height Model, product DP3.30015.001).
#'
#' @examples
#' data(harvardshapefile)
#'
#' # View structure
#' str(harvardshapefile)
#'
#' # Plot canopy height by crown polygon
#' library(sf)
#' plot(harvardshapefile["Max_Height"])
#' @docType data
#' @name harvardshapefile
NULL
