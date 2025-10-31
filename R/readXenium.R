#' readXeniumSPE
#' @name readXeniumSPE
#' @rdname readXeniumSPE
#' @title Load data from a 10x Genomics Xenium experiment
#'
#' @description
#' Creates a [`SpatialExperiment`] from an unzipped Xenium Output Bundle
#' directory containing spatial gene expression data.
#'
#' @param dirName `character(1)`
#'   Path to a Xenium Output Bundle directory.
#' @param sampleName `character(1)`
#'   Sample identifier to assign to `sample_id`. Default: `"sample01"`.
#' @param type `character(1)`
#'   One of `"HDF5"` or `"sparse"`; method to read the feature matrix.
#' @param coordNames `character(2)`
#'   Names of X/Y spatial coordinate columns. Default:
#'   `c("x_centroid", "y_centroid")`.
#' @param boundariesType `character(1)`
#'   One of `"parquet"` or `"csv"`; format of the polygon file.
#' @param computeMissingMetrics `logical(1)`
#'   If `TRUE`, compute area and aspect‐ratio from boundary polygons.
#' @param keepPolygons `logical(1)`
#'   If `TRUE`, append raw polygon geometries to `colData`.
#' @param countsFilePattern `character(1)`
#'   Pattern to locate the feature matrix file. Default:
#'   `"cell_feature_matrix"`.
#' @param metadataFPattern `character(1)`
#'   Pattern to locate the cell metadata file. Default: `"cells"`.
#' @param polygonsFPattern `character(1)`
#'   Pattern to locate the cell boundaries file. Default:
#'   `"cell_boundaries"`.
#' @param polygonsCol `character(1)`
#'   Name of the polygons column to add to `colData`. Default:
#'   `"polygons"`.
#' @param txPattern `character(1)`
#'   Pattern (base filename, without extension) to locate the transcript file
#'   (usually a `.parquet` file) from which to extract Field-Of-View (FOV)
#'   information for each cell. Default: `"transcripts"`.
#' @param addFOVs `logical(1)`
#'   If `TRUE`, extract Field-Of-View (FOV) information from the transcript file
#'   (as located by `txPattern`) and append it to cell metadata (`colData`).
#'   Default: `FALSE`.

#'
#' @details
#' Expects the unzipped bundle to contain an `outs/` folder with:
#' - `cell_feature_matrix.h5` or `cell_feature_matrix/`
#' - `cells.csv.gz`
#'
#' @return A [`SpatialExperiment`] object with assays, `colData`, spatial
#'   coordinates, and `metadata$polygons` & `metadata$technology`.
#'
#' @author Dario Righelli, Benedetta Banzi
#'
#' @importFrom DropletUtils read10xCounts
#' @importFrom data.table fread
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SummarizedExperiment assays rowData
#' @importFrom SpatialExperimentIO readXeniumSXE
#' @export
#' @examples
#' xepath <- system.file(
#'   "extdata", "Xenium_small", package = "SpaceTrooper"
#' )
#' (spe <- readXeniumSPE(
#'   dirName = xepath,
#'   keepPolygons = TRUE
#' ))
readXeniumSPE <- function(dirName, sampleName="sample01",
    type=c("HDF5", "sparse"), coordNames=c("x_centroid", "y_centroid"),
    boundariesType=c("parquet", "csv"), computeMissingMetrics=TRUE,
    keepPolygons=FALSE, countsFilePattern="cell_feature_matrix",
    metadataFPattern="cells.csv.gz", polygonsFPattern="cell_boundaries",
    polygonsCol="polygons", txPattern="transcripts", addFOVs=FALSE) {

    stopifnot(file.exists(dirName))
    type <- match.arg(type)
    boundariesType <- match.arg(boundariesType)
    if(basename(dirName) != "outs") { # add "outs/" dir if not already included
        dirbkup <- dirName
        dirName <- file.path(dirName, "outs")
        if (!file.exists(dirName)) {
            dirName <- dirbkup
        } else {
            warning("automatically detected/added outs dir in the 10x filepath")
        }
    }

    cfm <- paste0(countsFilePattern, switch(type, HDF5=".h5", ""))
    spe <- SpatialExperimentIO::readXeniumSXE(dirName=dirName,
            countMatPattern=cfm, metaDataPattern=metadataFPattern,
            coordNames=coordNames, returnType="SPE", addExperimentXenium=FALSE,
            altExps=NULL, addParquetPaths=FALSE)
    spe$sample_id <- sampleName
    rownames(colData(spe)) <- spe$cell_id
    pex <- paste0(polygonsFPattern, switch(boundariesType, parquet=".parquet",
                                                            csv=".csv.gz"))
    polfile <- list.files(dirName, pex, full.names=TRUE)

    if (computeMissingMetrics) {
        message("Computing missing metrics, this could take some time...")
        cd <- computeMissingMetricsXenium(polfile, colData(spe), keepPolygons,
                                        polygonsCol)
    }
    if (addFOVs) {
        cd <- .addFovFromTx(file.path(dirName, paste0(txPattern, ".parquet")),
                    cd)
    }

    colData(spe) <- cd
    metadata(spe) <- list(polygons=polfile, technology="10X_Xenium")
    return(spe)
}

#' computeMissingMetricsXenium
#' @name computeMissingMetricsXenium
#' @rdname computeMissingMetricsXenium
#'
#' @description
#' Compute Missing Metrics for Xenium Data
#'
#' This function computes missing metrics, such as the aspect ratio, from
#' polygon data in a Xenium dataset and optionally appends the polygon data to
#' the resulting `colData`.
#'
#' @param polFile A character string specifying the file path to the polygon
#' data.
#' @param colData A `DataFrame` containing the `colData` for the Xenium dataset.
#' @param keepPolygons A logical value indicating whether to keep the polygon
#' data in the resulting `colData`. Default is `FALSE`.
#' @param polygonsCol character indicating the name of the polygons column to
#' add into the colData (default is `polygons`).
#'
#' @return A `DataFrame` containing the updated `colData` with computed metrics.
#' If `keepPolygons` is `TRUE`, the polygon data is also included.
#'
#' @details The function reads the polygon data from the specified file,
#' computes the aspect ratio for each polygon, and merges these metrics with
#' the provided `colData`. Optionally, the polygon data can be kept in the
#' returned `colData`.
#'
#' @importFrom S4Vectors cbind.DataFrame
#' @export
#'
#' @examples
#' example(readXeniumSPE)
#' colData(spe) <- computeMissingMetricsXenium(metadata(spe)$polygons,
#'     colData(spe), keepPolygons=TRUE)
computeMissingMetricsXenium <- function(polFile, colData, keepPolygons=FALSE,
                                polygonsCol="polygons")
{
    stopifnot(file.exists(polFile))
    polygons <- readPolygonsXenium(polFile, keepMultiPol=TRUE)
    cd <- colData
    cd$AspectRatio <- computeAspectRatioFromPolygons(polygons)
    if(keepPolygons) cd <- .addPolygonsToCD(cd, polygons, polygonsCol)
    return(cd)
}

#' .addFovFromTx
#' @name dot-addFovFromTx
#' @rdname dot-addFovFromTx
#'
#' @description
#' Add FOV information from transcript file to cell metadata.
#'
#' This function retrieves FOV information from transcript file and appends
#' the data to the resulting `colData`.
#'
#' @param txFile `character(1)` path to a Xenium Output tx file.
#' @param colData A `DataFrame` containing the `colData` for the Xenium dataset.
#'
#' @return A `DataFrame` containing the updated `colData` with FOV information.
#'
#' @details The function reads the transcript file then groups it by cell_id
#' and merges the FOV information to the cell metadata in `colData`. Only
#' parquet file is supported for this operation
#'
#' @importFrom S4Vectors cbind.DataFrame
#' @importFrom arrow read_parquet
#' @importFrom dplyr group_by select distinct left_join
#' @keywords internal
#'
.addFovFromTx <- function(txFile, colData) {
    stopifnot(file.exists(txFile))
    df <- data.frame(colData)
    tx <- arrow::read_parquet(txFile)
    if (!"fov_name" %in% colnames(tx)) {
        stop("No fov_name column was found in tx file. \r\n",
            "Rerun readXeniumSPE without adding FOV information.")}
    g_tx <- group_by(tx, cell_id) |> select(cell_id, fov=fov_name) |>
        distinct(cell_id, .keep_all = TRUE)
    df <- left_join(df, g_tx, by="cell_id")
    colData$fov <- df$fov
    return(colData)
}

