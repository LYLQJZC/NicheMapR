#' Global microclimate model using local raster climate inputs
#'
#' Convenience wrapper around \code{micro_global()} for supplying climate rasters
#' from a local directory.
#'
#' @param spatial Path to a local directory containing \code{global_climate.nc} and,
#'   when \code{runmoist = 0}, \code{soilw.mon.ltm.v2.nc}.
#' @param ... Additional arguments passed directly to \code{micro_global()}.
#'
#' @return Output from \code{micro_global()}.
#' @export
micro_world <- function(spatial, ...) {
  micro_global(spatial = spatial, ...)
}
