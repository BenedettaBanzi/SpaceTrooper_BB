% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polygons.R
\name{readPolygonsMerfish}
\alias{readPolygonsMerfish}
\title{readPolygonsMerfish}
\usage{
readPolygonsMerfish(
  polygonsFolder,
  type = c("HDF5", "parquet"),
  keepMultiPol = TRUE,
  hdf5pattern = "hdf5",
  z_lev = 3L,
  zcolumn = "ZIndex",
  geometry = "Geometry",
  verbose = FALSE
)
}
\arguments{
\item{polygonsFolder}{A character string specifying the folder containing the
polygon data files.}

\item{type}{A character string specifying the file type ("HDF5" or "parquet").}

\item{keepMultiPol}{A logical value indicating whether to keep multipolygons.}

\item{hdf5pattern}{A character string specifying the pattern to match HDF5
files.}

\item{z_lev}{An integer specifying the Z level to filter the data. Default is
`3L`.}

\item{zcolumn}{A character string specifying the column name for the Z index.}

\item{geometry}{A character string specifying the geometry column name.}

\item{verbose}{A logical value indicating whether to print additional
information.}
}
\value{
An `sf` object containing the MERFISH polygon data.
}
\description{
This function reads polygon data specific to MERFISH technology.
}
\examples{
# Read MERFISH polygon data from a Parquet file:
# polygons <- readPolygonsMerfish("path/to/polygon_folder", type="parquet")
}
