#' Logarithmic binning function
#'
#' This function performs logarithmic binning of values along the x-axis,
#' optionally aggregating y-values (e.g., biomass) into the bins.
#'
#' @param x Numeric vector of values to be binned in the x dimension (e.g., diameter).
#' @param y Optional numeric vector of values to be binned in the y dimension
#'   (defaults to \code{NULL}). Typically something like biomass production.
#' @param n Integer, number of bins.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{bin_midpoint} Midpoint of bins on x-axis (arithmetic mean in linear space).
#'   \item \code{bin_value} Height of bins, either counts per unit width (if \code{y=NULL})
#'         or sum(y)/width.
#'   \item \code{bin_count} Number of individuals in bins.
#'   \item \code{bin_min} Lower edge of bins on x-axis.
#'   \item \code{bin_max} Upper edge of bins on x-axis.
#' }
#'
#' @details
#' This implementation is adapted from the \emph{forestscaling} package
#' (not available on CRAN). Please cite the original package if you use this
#' function in research outputs.
#'
#' @examples
#' set.seed(123)
#' x <- runif(100, 1, 100)
#' y <- rpois(100, lambda = 5)
#' logbin(x, y, n = 5)
#'
#' @export
logbin <- function(x, y = NULL, n) {
  logx <- log10(x)
  bin_edges <- seq(min(logx), max(logx), length.out = n + 1)
  logxbin <- rep(NA, length(logx))
  b <- bin_edges
  b[length(b)] <- b[length(b)] + 1
  for (i in 1:length(logx)) {
    logxbin[i] <- sum(logx[i] >= b)
  }
  bin_midpoints <- numeric(n)
  for (i in 1:n) {
    bin_midpoints[i] <- mean(10^(bin_edges[i:(i + 1)]))
  }
  bin_widths <- diff(10^bin_edges)
  bin_factor <- factor(logxbin, levels = 1:n)
  bin_counts <- table(bin_factor)
  if (!is.null(y)) {
    rawy <- tapply(y, bin_factor, sum)
    rawy[is.na(rawy)] <- 0
    bin_values <- as.numeric(rawy / bin_widths)
  } else {
    bin_values <- as.numeric(bin_counts / bin_widths)
  }
  return(data.frame(
    bin_midpoint = bin_midpoints,
    bin_value = bin_values,
    bin_count = as.numeric(bin_counts),
    bin_min = 10^bin_edges[1:n],
    bin_max = 10^bin_edges[2:(n + 1)]
  ))
}
