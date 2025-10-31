#' readCosmxSPE
#' @name readCosmxSPE
#' @rdname readCosmxSPE
#' @aliases readCosmxSPE readCosmxProteinSPE
#' @description
#' Read and Construct a SpatialExperiment Object from CosMx Data
#'
#' This function reads in data from Nanostring CosMx files and constructs a
#' `SpatialExperiment` object, optionally including polygon data.
#'
#' @param dirName A character string specifying the directory containing the
#' CosMx data files.
#' @param sampleName A character string specifying the sample name. Default is
#' `"sample01"`.
#' @param coordNames A character vector specifying the names of the spatial
#' coordinate columns in the data. Default is `c("CenterX_global_px",
#' "CenterY_global_px")`.
#' @param countMatFPattern A character string specifying the pattern to match
#' the count matrix file. Default is `"exprMat_file.csv"`.
#' @param metadataFPattern A character string specifying the pattern to match
#' the metadata file. Default is `"metadata_file.csv"`.
#' @param polygonsFPattern A character string specifying the pattern to match
#' the polygons file. Default is `"polygons.csv"`.
#' @param fovPosFPattern A character string specifying the pattern to match the
#' FOV positions file. Default is `"fov_positions_file.csv"`.
#' @param fovdims A named numeric vector specifying the dimensions of the FOV
#' in pixels. Default is `c(xdim=4256, ydim=4256)`.
#'
#' @return A `SpatialExperiment` object containing the read CosMx data,
#' including count matrices, metadata, and optionally polygons.
#'
#' @details The function reads in the specified files for count matrices,
#' metadata, and FOV positions, and constructs a `SpatialExperiment` object.
#' Optionally, polygon data can be read and added to the object.
#'
#' readCosmxProteinSPE is a wrapper of readCosmxSPE, it only changes the
#' technology metadata in Nanostring_CosMx_Protein.
#'
#' @author Dario Righelli, Benedetta Banzi
#'
#' @importFrom data.table fread merge.data.table
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom dplyr mutate
#' @importFrom SpatialExperimentIO readCosmxSXE
#' @export
#'
#' @examples
#' cospath <- system.file(file.path("extdata", "CosMx_DBKero_Tiny"),
#'    package="SpaceTrooper")
#' spe <- readCosmxSPE(cospath, sampleName="DBKero_Tiny")
readCosmxSPE <- function(dirName, sampleName="sample01",
    coordNames=c("CenterX_global_px", "CenterY_global_px"),
    countMatFPattern="exprMat_file.csv", metadataFPattern="metadata_file.csv",
    polygonsFPattern="polygons.csv", fovPosFPattern="fov_positions_file.csv",
    fovdims=c(xdim=4256, ydim=4256)) {

    stopifnot(all(names(fovdims) == c("xdim", "ydim"), file.exists(dirName)))

    spe <- SpatialExperimentIO::readCosmxSXE(dirName=dirName, returnType="SPE",
        countMatPattern=countMatFPattern, metaDataPattern=metadataFPattern,
        coordNames=coordNames, addFovPos=TRUE, fovPosPattern=fovPosFPattern,
        altExps=NULL,addParquetPaths=FALSE)

    pol_file <- list.files(dirName, polygonsFPattern, full.names=TRUE)#parquet?
    cn <- paste0("f", spe$fov, "_c", spe$cell_ID)
    colnames(spe) <- cn
    rownames(colData(spe)) <- cn
    if(length(grep("cell_id", colnames(colData)))!=0)
        warning("Overwriting existing cell_id column in colData")
    spe$cell_id <- cn
    spe <- .checkFovPositionVersion(spe)
    metadata(spe) <- list(fov_positions=metadata(spe)$fov_positions,
        fov_dim=fovdims, polygons=pol_file, technology="Nanostring_CosMx")

    names(colData(spe))[names(colData(spe)) == "cell_ID"] <- "cellID"
    spe$sample_id <- sampleName
    return(spe)
}

#' @export
readCosmxProteinSPE <- function(dirName, sampleName="sample01",
    coordNames=c("CenterX_global_px", "CenterY_global_px"),
    countMatFPattern="exprMat_file.csv", metadataFPattern="metadata_file.csv",
    polygonsFPattern="polygons.csv", fovPosFPattern="fov_positions_file.csv",
    fovdims=c(xdim=4256, ydim=4256)) {

    spe <- readCosmxSPE(dirName, sampleName, coordNames, countMatFPattern,
        metadataFPattern, polygonsFPattern, fovPosFPattern, fovdims)

    metadata(spe)$technology <- "Nanostring_CosMx_Protein"
    return(spe)
}

#' Check and Standardize FOV Position Column Names
#'
#' This internal utility function standardizes column names of a data frame
#' containing Field of View (FOV) positional information.
#' It modifies column names to ensure compatibility with expected naming
#' conventions, including support for older formats.
#'
#' Specifically, it:
#' - Renames any column containing "FOV" to "fov"
#' - Converts columns with coordinates matching "X", "Y", or "Z" to lowercase
#' - Replaces suffix "_px" with "_global_px" for coordinate pixel columns
#' - If the input contains `x_mm` and `y_mm` columns, the function computes
#' corresponding `x_global_px` and `y_global_px` values by converting from
#' millimeters to pixels using a fixed resolution factor (0.12028 mm/pixel).
#' @param spe A `SpatialExperiment` containing FOV position information
#' in the metadata to be standardized.
#'
#' @return A `SpatialExperiment` with updated and standardized column names
#' for the metadata `fov_position` `data.frame`.
#'
#' @keywords internal
#' @noRd

.checkFovPositionVersion <- function(spe)
{
    fovpos <- metadata(spe)$fov_positions
    fovcidx <- grep("FOV", colnames(fovpos)) # works also with older vers
    if(length(fovcidx)!=0) colnames(fovpos)[fovcidx] <- "fov"
    fovcrdx <- grep("[X|Y|Z]", colnames(fovpos))
    if(length(fovcrdx)!=0) colnames(fovpos)[fovcrdx] <-
        tolower(colnames(fovpos)[fovcrdx])
    fovccdx <- grep("[x|y]_px", colnames(fovpos))
    if(length(fovccdx)!=0) colnames(fovpos)[fovccdx] <-
        gsub("_px", "_global_px", colnames(fovpos)[fovccdx])

    if(length(grep("x_mm", colnames(fovpos))!=0)) {
        fovpos <- fovpos |>
            dplyr::mutate(x_global_px = x_mm/0.12028*10^3,
                        y_global_px = (y_mm/0.12028*10^3) - 4256)
    }
    idx <- fovpos$fov %in% unique(spe$fov)
    fovpos <- fovpos[idx, ]

    fovpos <- fovpos[order(fovpos$fov), ]
    metadata(spe)$fov_positions <- fovpos
    return(spe)
}
