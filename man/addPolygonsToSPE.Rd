% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polygons.R
\name{addPolygonsToSPE}
\alias{addPolygonsToSPE}
\title{addPolygonsToSPE}
\usage{
addPolygonsToSPE(spe, polygons)
}
\arguments{
\item{spe}{A `SpatialExperiment` object to which polygons will be added.}

\item{polygons}{An `sf` object containing the polygon data.}
}
\value{
The `SpatialExperiment` object with polygons added to the `colData`.
}
\description{
This function adds polygon data to a `SpatialExperiment` object.
}
\examples{
# Assuming `spe` is a SpatialExperiment object and `polygons` is an sf object:
# spe <- addPolygonsToSPE(spe, polygons)
}
