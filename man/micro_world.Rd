% Generated manually
\name{micro_world}
\alias{micro_world}
\title{Global microclimate model using local raster climate inputs}
\usage{
micro_world(spatial, ...)
}
\arguments{
\item{spatial}{Path to a local directory containing \code{global_climate.nc} and, when
\code{runmoist = 0}, \code{soilw.mon.ltm.v2.nc}.}

\item{...}{Additional arguments passed to \code{micro_global()}.}
}
\value{
Output list returned by \code{micro_global()}.
}
\description{
Convenience wrapper around \code{micro_global()} for supplying climate rasters
from a local directory.
}
