#' spatialPerCellQC
#' @name spatialPerCellQC
#' @rdname spatialPerCellQC
#' @description
#' Computes quality‐control metrics for each cell and adds them to `colData`.
#'
#' @param spe A `SpatialExperiment` object containing spatial data.
#' @param micronConvFact Numeric factor to convert pixels to microns. Default
#'   `0.12`.
#' @param rmZeros logical for removing zero counts cells (default is TRUE).
#' @param negProbList Character vector of patterns to identify negative probes.
#'   Defaults include:
#'   - Nanostring CosMx: `"NegPrb"`, `"Negative"`, `"SystemControl"`
#'   - Xenium: `"NegControlProbe"`, `"NegControlCodeword"`,
#'     `"UnassignedCodeword"`
#'   - MERFISH: `"Blank"`
#' @param use_altexps logical for `use_altexps` in `scuttle` package.
#' If TRUE uses the altexps for computing some metrics on it.
#' Useful for interoperability with `SpatialExperimentIO`.
#' (See \link[scuttle]{addPerCellQC} for additional details).
#'
#' @return A `SpatialExperiment` object with added QC metrics in `colData`.
#'
#' @details
#' Calculates sums and detected counts for control and target probes,
#' computes ratio and count‐area metrics, converts coords to microns for
#' CosMx, and drops zero‐count cells.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom scater addPerCellQC
#' @importFrom S4Vectors cbind.DataFrame
#' @export
#' @examples
#' example(readCosmxSPE)
#' spe <- spatialPerCellQC(spe)
spatialPerCellQC <- function(spe, micronConvFact=0.12, rmZeros=TRUE,
    negProbList=c("NegPrb", "Negative", "SystemControl", "Ms IgG1", "Rb IgG",
    "BLANK_", "NegControlProbe", "NegControlCodeword", "UnassignedCodeword",
    "Blank"),
    use_altexps=NULL) {
    stopifnot(is(object=spe, "SpatialExperiment"))
    idxlist <- lapply(negProbList, function(ng) {
        grep(paste0("^", ng), rownames(spe))
    })
    names(idxlist) <- negProbList
    idxlist <- idxlist[which(lengths(idxlist)!=0)]
    spe <- addPerCellQC(spe, subsets=idxlist, use_altexps=use_altexps)
    idx <- grep("^subsets_.*_sum$", colnames(colData(spe)))
    npc <- npd <- 0
    if ( length(idx) !=0 ) {
        npc <- rowSums(as.matrix(colData(spe)[ , idx, drop=FALSE]))
        # TODO: not robust at all! the +1 is not a really good choice
        npd <- rowSums(as.matrix(colData(spe)[ , idx+1, drop=FALSE]))
    }
    spe$control_sum <- npc
    spe$control_detected <- npd
    spe$target_sum <- spe$sum - npc
    spe$target_detected <- spe$detected - npd
    if(!all(spatialCoordsNames(spe) %in% names(colData(spe)))) {
        # TODO: CHANGE SPE constructor WITH COORDINATES IN COLDATA
        colData(spe) <- cbind.DataFrame(colData(spe), spatialCoords(spe))
    }

    spe$ctrl_total_ratio <- spe$control_sum/spe$total
    spe$ctrl_total_ratio[which(is.na(spe$ctrl_total_ratio))] <- 0
    if(metadata(spe)$technology == "Nanostring_CosMx_Protein") {
        # Only for proteins will be included in QScore
        spe$log2Ctrl_total_ratio <- log2(spe$ctrl_total_ratio)
        idx <- which(names(colData(spe)) == "Area.um2")
        if(length(idx)!=0) { names(colData(spe))[idx] <- "Area_um" }
    }

    if(any(metadata(spe)$technology %in%
            c("Nanostring_CosMx", "Nanostring_CosMx_Protein"))) {
        spnc <- spatialCoords(spe) * micronConvFact
        colnames(spnc) <- gsub("px", "um", spatialCoordsNames(spe))
        colData(spe) <- cbind.DataFrame(colData(spe), spnc)
        spe$Area_um <- spe$Area * (micronConvFact^2)
        spe <- .computeBorderDistanceCosMx(spe)
    }

    if (metadata(spe)$technology == "10X_Xenium") {
        spe$Area_um <- spe$cell_area # standardized across other techs
    }
    if ("AspectRatio" %in% colnames(colData(spe))) {
        spe$log2AspectRatio <- log2(spe$AspectRatio) # not cosmx
    } else { warning("Missing aspect ratio in colData") }

    spe$CountArea <- spe$sum/spe$Area_um
    spe$log2CountArea <- log2(spe$CountArea)
    if (rmZeros) {
        if (sum(spe$sum==0) > 0) {
            message("Removing ", dim(spe[,spe$sum==0])[2],
                    " cells with 0 counts!")
            spe <- spe[,!spe$sum==0]
        }
    }
    return(spe)
}

#' .computeBorderDistanceCosMx
#' @name .computeBorderDistanceCosMx
#' @rdname dot-computeBorderDistanceCosMx
#' @description
#' Calculates the minimum distance of each cell to the field‐of‐view border
#' and adds it to `colData`.
#'
#' @param spe A `SpatialExperiment` object with CosMx data.
#' @param xwindim Width of FOV in x (default from `metadata(spe)$fov_dim`).
#' @param ywindim Height of FOV in y (default from `metadata(spe)$fov_dim`).
#'
#' @return A `SpatialExperiment` object with `dist_border` columns in
#' `colData`.
#'
#' @importFrom dplyr left_join
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @keywords internal
.computeBorderDistanceCosMx <- function(spe,
                                    xwindim=metadata(spe)$fov_dim[["xdim"]],
                                    ywindim=metadata(spe)$fov_dim[["ydim"]]) {
    stopifnot(is(spe, "SpatialExperiment"))
    cd <- colData(spe)
    cdf <- left_join(as.data.frame(cd),metadata(spe)$fov_positions,by="fov")
    spcn <- spatialCoordsNames(spe)
    fovpn <- colnames(metadata(spe)$fov_positions)[colnames(
        metadata(spe)$fov_positions) %in% c("x_global_px", "y_global_px")]
    cd$dist_border_x <- pmin(cdf[,spcn[1]] - cdf[,fovpn[1]],
                            (cdf[,fovpn[1]] + xwindim) - cdf[,spcn[1]])
    cd$dist_border_y <- pmin(cdf[,spcn[2]] - cdf[,fovpn[2]],
                            (cdf[,fovpn[2]] + ywindim) - cdf[,spcn[2]])
    cd$dist_border <- pmin(cd$dist_border_x, cd$dist_border_y)
    colData(spe) <- cd
    return(spe)
}

#' computeSpatialOutlier
#' @name computeSpatialOutlier
#' @rdname computeSpatialOutlier
#' @description
#' Computes outliers based on the Area (in micron) of the experiment.
#' It gives the possibility to choose between the medcouple (mc method argument)
#' and the MADs (scuttle method argument).
#'
#' @details
#' The medcouple method is a measure for the skeweness of univariate
#' distribution as described in Hubert M. et al. (2008).
#' In particular, the computed medcouple value must be in a range between -0.6
#' and 0.6 to computed adjusted boxplots and perform the outlier detection.
#' For median absolute deviations (MADs) method we just wrap the isOutlier
#' function in the scuttle package. Please see McCarthy DJ et al (2017)
#' for further details.
#'
#' @param spe a SpatialExperiment object with target_counts, area in micron
#' and log2 of the aspect ratio in the `colData`.
#' @param computeBy character indicating a `colData` column name on which
#' compute the outlier.
#' @param method one of `mc`, `scuttle`, `both`.
#' Use `mc` for medcouple, `scuttle` for median absolute deviations as computed
#' in `scuttle`, `both` for computing both of them.
#' @param mcDoScale logical indicating if the values to compute the medcouple
#' for the outlier detection should be scaled (default is FALSE, as suggested
#' by the original Medcouple authors.). See \link[robustbase]{mc} for further
#' readings.
#' @param scuttleType One of `"both"`, `"lower"`, `"higher"` for scuttle method.
#'
#' @return a SpatialExperiment object with additional column(s) (named as
#' the column name indicated in `column_by` followed by the outlier_sc/mc
#' nomenclature) with the outlier detection as `outlier.filter` logical class
#' object. This allows to store the thresholds as attributes of the column.
#' use attr(,"thresholds") to retrieve them.
#'
#' @export
#' @importFrom robustbase mc adjbox
#' @importFrom e1071 skewness
#' @importFrom scuttle isOutlier outlier.filter
#'
#' @examples
#' example(spatialPerCellQC)
#' spe <- computeSpatialOutlier(spe, computeBy="log2CountArea", method="both")
#' table(spe$log2CountArea_outlier_mc)
#' table(spe$log2CountArea_outlier_sc)
computeSpatialOutlier <- function(spe, computeBy=NULL,
    method=c("mc", "scuttle", "both"), mcDoScale=FALSE,
    scuttleType=c("both", "lower", "higher")) {
    stopifnot(all(is(spe, "SpatialExperiment"), !is.null(computeBy)))
    stopifnot(computeBy %in% names(colData(spe)))
    options(mc_doScale_quiet=TRUE)
    method <- match.arg(method)
    scuttleType <- match.arg(scuttleType)
    cd <- colData(spe)
    cdcol <- cd[[computeBy]]
    mcfl <- scuttlefl <- FALSE
    switch(method, both={ mcfl <- scuttlefl <- TRUE },
            mc={ mcfl <- TRUE }, scuttle={ scuttlefl <- TRUE },
            {stop("Method is not one of allowed methods")} )
    if (mcfl) {
        skw <- e1071::skewness(cdcol, na.rm = TRUE) # NAs arise problems
        if (skw >- 1 & skw < 1) warning("Distribution is symmetric: ",
                "mc is for asymmetric distributions. Use scuttle instead.")
        mcval <- robustbase::mc(cdcol, doScale=mcDoScale, na.rm=TRUE)
        if ( any( (mcval <= -0.6), (mcval >= 0.6) ) )
            stop("mc is: ",round(mcval, digits=4),"outliers reqs not satisfied")
        names(cdcol) <- colnames(spe)
        outl <- robustbase::adjbox(cdcol, plot=FALSE)
        outsmc <- rep("NO", dim(cd)[1])
        outsmc[rownames(cd) %in% names(outl$out)] <-
            ifelse(outl$out <= outl$fence[1], "LOW", "HIGH")
        outlier_mc <- scuttle::outlier.filter(outsmc) #using scuttle class
        thrs <- as.numeric(outl$fence)
        names(thrs) <- c("lower", "higher")
        attr(outlier_mc, "thresholds") <- thrs
        cd$outlier_mc <- outlier_mc
        names(cd)[names(cd) =="outlier_mc"] <- paste0(computeBy, "_outlier_mc")
        # TODO: compute distributions in the adjusted boxplots to store in cd
    }
    if (scuttlefl) {
        outssc <- scuttle::isOutlier(cdcol, type=scuttleType)
        sctri <- rep("NO", dim(cd)[1])
        sctri <- ifelse(outssc == TRUE & cdcol <= attr(outssc, "thresholds")[1],
                        "LOW", sctri)
        outlier_sc <- ifelse(outssc == TRUE &
                            cdcol >= attr(outssc, "thresholds")[2],
                            "HIGH", sctri)
        outlier_sc <- scuttle::outlier.filter(outlier_sc)
        attr(outlier_sc, "thresholds") <- attr(outssc, "thresholds")
        cd$outlier_sc <- outlier_sc
        names(cd)[names(cd)=="outlier_sc"] <- paste0(computeBy, "_outlier_sc")
    }
    colData(spe) <- cd
    return(spe)
}

#' computeThresholdFlags
#' @name computeThresholdFlags
#' @rdname computeThresholdFlags
#' @description
#' Compute Flagged cells using fixed thresholds for SpatialExperiment.
#'
#' This function calculates flagged cells only for total counts and control on
#' total probe counts ratio using fixed thresholds for a `SpatialExperiment`
#' object.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#' @param totalThreshold A numeric value for the threshold of total counts to
#' identify cells with low counts. Default is `0`.
#' @param ctrlTotRatioThreshold A numeric value for the threshold of
#' control-to-total ratio to flag cells over a certain threshold. Default is
#' `0.1`.
#'
#' @return The `SpatialExperiment` object with added filter flags in `colData`.
#'
#' @details The function flags cells basing on zero counts and control-to-total
#' ratio to identify junk cells.
#' It also combines these flags into a single filter flag.
#'
#' @importFrom SummarizedExperiment colData
#' @export
#' @examples
#' example(readCosmxSPE)
#' spe <- spatialPerCellQC(spe)
#' spe <- computeThresholdFlags(spe)
#' table(spe$threshold_flags)
computeThresholdFlags <- function(spe, totalThreshold=0,
                            ctrlTotRatioThreshold=0.1)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("total" %in% names(colData(spe)))
    stopifnot("ctrl_total_ratio" %in% names(colData(spe)))

    spe$is_zero_counts <- ifelse(spe$total == totalThreshold, TRUE, FALSE)
    #flagging cells with probe counts on total counts ratio > 0.1
    spe$is_ctrl_tot_outlier <- ifelse(spe$ctrl_total_ratio >
                                        ctrlTotRatioThreshold, TRUE, FALSE)

    spe$threshold_flags <- (spe$is_ctrl_tot_outlier &
                                spe$is_zero_counts)
    return(spe)
}


#' computeLambda
#' @description
#' Compute Optimal Ridge Regularization Parameter \eqn{\lambda} via
#' Cross-Validation
#'
#' \code{computeLambda} performs ridge (L2) logistic regression with
#' cross-validation to identify the optimal regularization parameter
#' \eqn{\lambda} for a binary response.
#'
#' @param technology  \[character\]
#'   The name of the experimental technology. Passed to
#'   \code{getModelFormula()} to retrieve the corresponding model formula.
#'
#' @param trainDF  \[data.frame\]
#'   A data frame for training that must include:
#'   \describe{
#'     \item{Predictor columns}{All columns referenced in the formula returned
#'     by \code{getModelFormula()}.}
#'     \item{\code{qscore_train}}{A binary (0/1) response vector to be modeled.}
#'   }
#'
#' @return
#' \[numeric\]
#'   The value of \eqn{\lambda} (i.e., \code{lambda.min}) from
#'   \code{\link[glmnet]{cv.glmnet}} that minimizes the cross-validation error.
#'
#' @details
#' Internally, the function:
#' \enumerate{
#'   \item Calls \code{getModelFormula(technology)} to obtain a model formula
#'   as text,
#'   \item Constructs the design matrix via \code{model.matrix()},
#'   \item Runs ridge logistic regression cross-validation using
#'         \code{\link[glmnet]{cv.glmnet}} with \code{alpha = 0},
#'   \item Extracts and returns \code{ridge_cv$lambda.min}.
#' }
#'
#' @examples
#' example(spatialPerCellQC)
#' withr::with_seed(1998, trainDF <- computeTrainDF(spe))
#' best_lambda <- computeLambda(metadata(spe)$technology, trainDF)
#' print(best_lambda)
#'
#' @seealso
#' \code{\link[glmnet]{cv.glmnet}}
#'
#' @export
computeLambda <- function(technology, trainDF) {
    model_formula <- getModelFormula(technology)
    model_matrix <- model.matrix(as.formula(model_formula), data=trainDF)
    ridge_cv <- cv.glmnet(model_matrix, trainDF$qscore_train,
                        family="binomial", alpha=0, lambda=NULL)
    bestLambda <- ridge_cv$lambda.min
    return(bestLambda)
}

#' computeQScore
#' @name computeQScore
#' @rdname computeQScore
#' @description
#' Compute quality score and automatically define weights for quality score
#' through glm training. This function computes quality score with a formula
#' that depends on the technology.
#'
#' @details
#' For CosMx datasets, the Quality Score formula is defined as follows:
#'
#' quality score ~ count density - aspect ratio - interaction term
#'
#' count density is total counts-to-area ratio, aspect ratio represents
#' border effect typical of CosMx datasets and the last one is the
#' interaction term of the previous two terms.
#'
#' For Xenium and Merscope datasets, quality score depends solely on count
#' density, as no border effect has been observed for these two technologies.
#'
#' To automatically define the formula coefficient weights, model training
#' is performed through ridge regression.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#' @param verbose logical for having a verbose output. Default is FALSE.
#' @param bestLambda the best lambda typically computed using `computeLambda`.
#'
#' @return The `SpatialExperiment` object with added quality score in `colData`.
#' @export
#' @importFrom dplyr case_when filter mutate distinct
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats as.formula model.matrix predict quantile
#' @examples
#' example(spatialPerCellQC)
#' set.seed(1998)
#' spe <- computeQScore(spe)
#' summary(spe$training_status)
#' summary(spe$quality_score)
computeQScore <- function(spe, bestLambda=NULL, verbose=FALSE) {
    stopifnot(is(spe, "SpatialExperiment"))

    trainDF <- computeTrainDF(spe, verbose)
    model_formula <- getModelFormula(metadata(spe)$technology)
    model_matrix <- model.matrix(as.formula(model_formula), data=trainDF)
    model <- trainModel(model_matrix, trainDF)
    if(is.null(bestLambda)) {
        bestLambda <- computeLambda(metadata(spe)$technology,
                                trainDF)
    }
    cd <- data.frame(colData(spe))
    full_matrix <- model.matrix(as.formula(model_formula), data = cd)
    cd$quality_score <- as.vector(predict(model, s=bestLambda,
                                        newx = full_matrix,
                                        type = "response"))
    spe$quality_score <- cd$quality_score
    train_identity <- rep("TEST", dim(spe)[2])
    train_bad <- trainDF$cell_id[trainDF$qscore_train==0]
    train_good <- trainDF$cell_id[trainDF$qscore_train==1]
    spe$training_status <- dplyr::case_when(
        spe$cell_id %in% train_bad ~ "BAD",
        spe$cell_id %in% train_good ~ "GOOD",
        TRUE ~ train_identity)
    return(spe)
}

#' trainModel
#' @name trainModel
#' @rdname trainModel
#' @description
#' Fit a Ridge Logistic Regression Model
#'
#' \code{trainModel} fits an L2-regularized (ridge) logistic regression
#' using \pkg{glmnet}, given a design matrix and a training data frame.
#'
#' @param trainDF \[data.frame\]
#'   A data frame containing at least the response column
#'   \code{qscore_train}, coded as 0/1.
#' @param modelMatrix a matrix describing the model variables, tipically created
#' with `getModelFormula` and `model.matrix` functions.
#'
#' @return
#' A \code{\link[glmnet]{glmnet}} model object fitted with
#' \code{family="binomial"}, \code{alpha=0} (ridge), and a sequence of
#' \eqn{\lambda} values.
#'
#' @examples
#' example(computeTrainDF)
#' model_formula <- getModelFormula(metadata(spe)$technology)
#' model_matrix <- model.matrix(as.formula(model_formula), data=df_train)
#' fit <- trainModel(model_matrix, df_train)
#' coef(fit, s = 0.01)
#'
#' @export
trainModel <- function(modelMatrix, trainDF)
{
    model <- glmnet(x=modelMatrix, y=trainDF$qscore_train,
                    family="binomial", lambda=NULL, alpha=0)
    return(model)
}

#' computeTrainDF
#' @name computeTrainDF
#' @rdname computeTrainDF
#' @description
#' Build a Balanced Training Data Frame from a SpatialExperiment
#'
#' \code{computeTrainDF} takes a \linkS4class{SpatialExperiment} object,
#' flags spatial outliers on “log2CountArea”, then assembles a
#' balanced training set of “good” vs “bad” cells for subsequent model fitting.
#'
#' @param spe \linkS4class{SpatialExperiment}
#'   A SpatialExperiment containing at least:
#'   \itemize{
#'     \item assay(s) with nonzero \code{total} counts,
#'     \item \code{colData(spe)} columns including \code{log2CountArea},
#'     \code{dist_border}, etc.,
#'     \item \code{metadata(spe)$technology} indicating the platform.
#'   }
#'
#' @param verbose \[logical(1)\] (default \code{FALSE})
#'   If \code{TRUE}, prints the number of “bad” and “good” cells selected.
#'
#' @return
#' A \code{data.frame} with one row per cell, including:
#' \itemize{
#'   \item \code{qscore_train} (0/1) indicating “bad” vs “good”,
#'   \item relevant \code{colData} columns used for modeling.
#' }
#'
#' @details
#' Internally the function:
#' \enumerate{
#'    \item Filters out zero-count cells,
#'    \item Calls \code{computeSpatialOutlier()} on “log2CountArea” to get
#'    fences,
#'    \item Labels cells as “LOW”/“HIGH” outliers or “NO”,
#'    \item Delegates to either \code{.computeCosmxTrainSet()} or
#'    \code{.computeXenMerTrainSet()} based on \code{metadata(spe)$technology},
#'    \item Deduplicates and down-samples “good” cells to match the number of
#'    “bad” cells.
#' }
#'
#' @examples
#' example(spatialPerCellQC)
#' df_train <- computeTrainDF(spe, verbose = TRUE)
#' table(df_train$qscore_train)
#'
#' @export
computeTrainDF <- function(spe, verbose=FALSE)
{
    spe_temp <- computeSpatialOutlier(spe[,spe$total>0],
        computeBy="log2CountArea", method="both")

    if(getFencesOutlier(spe_temp, "log2CountArea_outlier_mc", "lower") <
        min(spe_temp$log2CountArea)) {
            low_thr <- quantile(spe$log2CountArea, probs = 0.01)
    } else {
        low_thr <- getFencesOutlier(spe_temp, "log2CountArea_outlier_mc",
                                    "lower")
    }

    high_thr <- getFencesOutlier(spe_temp, "log2CountArea_outlier_mc", "higher")
    spe$log2CountArea_outlier_train <- case_when(spe$total==0 ~ "NO",
        spe$log2CountArea<low_thr ~ "LOW",spe$log2CountArea>high_thr ~ "HIGH",
        TRUE ~ "NO")

    attr(spe$log2CountArea_outlier_train, "thresholds") <-
        getFencesOutlier(spe_temp, "log2CountArea_outlier_mc")
    attr(spe$log2CountArea_outlier_train, "thresholds")[1] <- low_thr

    if(metadata(spe)$technology == "Nanostring_CosMx") {
        ts <- .computeCosmxTrainSet(spe)
    }
    if(metadata(spe)$technology == "Nanostring_CosMx_Protein") {
        ts <- .computeCosmxProteinTrainSet(spe)
    }
    if(any(metadata(spe)$technology %in% c("10X_Xenium", "Vizgen_MERFISH"))) {
        ts <- .computeXenMerTrainSet(spe)
    }
    train_bad <- ts$bad
    train_good <- ts$good
    train_bad <- train_bad |> distinct(cell_id, .keep_all = TRUE)

    if(verbose) message("Chosen low quality examples: ", dim(train_bad)[1])

    train_good <- train_good |> distinct(cell_id, .keep_all = TRUE)
    train_good <- train_good[!train_good$is_a_bad_boy,]
    # set.seed(1998) # not needed if run is encapsulated in with_seed funct
    train_good <- train_good[sample(rownames(train_good), dim(train_bad)[1],
                                    replace=FALSE),]
    if(verbose) message("Chosen good quality examples: ", dim(train_good)[1])

    train_good$is_a_bad_boy <- NULL
    trainDF <- rbind(train_bad, train_good)
    trainDF <- trainDF |> distinct(cell_id, .keep_all = TRUE)
    return(trainDF)
}

#' getModelFormula
#' @name getModelFormula
#' @rdname getModelFormula
#' @description
#' Returns the right‐hand side of a model formula string based on technology.
#' @param technology \[character\]
#'   Technology name to decide which predictors to include.
#' @return \[character\]
#'   A one‐sided formula as a string (e.g. "~ log2CountArea + ...").
#' @export
#' @examples
#' example(spatialPerCellQC)
#' getModelFormula(metadata(spe)$technology)
getModelFormula <- function(technology)
{
    model_formula <- "~log2CountArea" # xen and merf
    if(technology == "Nanostring_CosMx") {
        model_formula <- paste0("~ log2CountArea + I(abs(log2AspectRatio) ",
                                "* as.numeric(dist_border<50)) + ",
                                " log2CountArea:I(abs(log2AspectRatio)",
                                "* as.numeric(dist_border<50))") #for cosmx
    }
    if(technology == "Nanostring_CosMx_Protein") {
    model_formula <- paste0("~ log2CountArea + I(abs(log2AspectRatio) ",
        "*as.numeric(dist_border<50)) + log2Ctrl_total_ratio + ",
        " log2CountArea:I(abs(log2AspectRatio) ",
        "* as.numeric(dist_border<50)) + log2CountArea:log2Ctrl_total_ratio ",
        " + log2Ctrl_total_ratio:I(abs(log2AspectRatio) ",
        " * as.numeric(dist_border<50))")
    }
    return(model_formula)
}

#' .computeXenMerTrainSet
#' @name dot-computeXenMerTrainSet
#' @rdname dot-computeXenMerTrainSet
#' @description
#' Internal: Build Training Set for Xenium & MERFISH
#' Splits a SpatialExperiment into “bad” vs “good” cells based on
#' pre-computed outlier labels on log2CountArea.
#' @param spe \linkS4class{SpatialExperiment}
#' @return
#' A list with elements \code{bad} and \code{good}, each a data.frame
#' with \code{qscore_train} and (for “good”) an \code{is_a_bad_boy} flag.
#' @keywords internal
.computeXenMerTrainSet <- function(spe)
{
    train_bad <- data.frame(colData(spe)) |>
        filter(log2CountArea_outlier_train=="LOW") |> mutate(qscore_train=0)
    train_good <- data.frame(colData(spe)) |>
        filter((log2CountArea > quantile(log2CountArea, probs = 0.90) &
                    log2CountArea < quantile(log2CountArea, probs = 0.99))) |>
        mutate(qscore_train=1, is_a_bad_boy=cell_id %in% train_bad$cell_id)
    return(list(bad=train_bad, good=train_good))
}

#' .computeCosmxTrainSet
#' @name dot-computeCosmxTrainSet
#' @rdname dot-computeCosmxTrainSet
#' @description
#' Internal: Build Training Set for CosMx
#' Splits a SpatialExperiment into “bad” vs “good” cells based on
#' outliers in aspect ratio near tissue border or low count area.
#' @param spe \linkS4class{SpatialExperiment}
#' @return
#' A list with elements \code{bad} and \code{good}, each a data.frame
#' with \code{qscore_train} and (for “good”) an \code{is_a_bad_boy} flag.
#' @keywords internal
.computeCosmxTrainSet <- function(spe)
{
    spe <- computeSpatialOutlier(spe, "log2AspectRatio", "scuttle")
    train_bad <- data.frame(colData(spe)) |>
        filter((log2AspectRatio_outlier_sc == "HIGH" & dist_border < 50) |
                (log2AspectRatio_outlier_sc == "LOW" & dist_border < 50) |
                log2CountArea_outlier_train == "LOW") |>
        mutate(qscore_train = 0)

    train_good <- data.frame(colData(spe)) |>
        filter((log2AspectRatio > quantile(log2AspectRatio, probs = 0.25) &
                    log2AspectRatio < quantile(log2AspectRatio, probs = 0.75) &
                    dist_border > 50) |
                    (log2CountArea > quantile(log2CountArea, probs = 0.90) &
                    log2CountArea < quantile(log2CountArea, probs = 0.99))) |>
        mutate(qscore_train=1, is_a_bad_boy=cell_id %in% train_bad$cell_id)
    return(list(bad=train_bad, good=train_good))
}

#' .computeCosmxProteinTrainSet
#' @name dot-computeCosmxProteinTrainSet
#' @rdname dot-computeCosmxProteinTrainSet
#' @description
#' Internal: Build Training Set for CosMx-Protein
#' Splits a SpatialExperiment into “bad” vs “good” cells based on
#' outliers in aspect ratio near tissue border or low count area.
#' @param spe \linkS4class{SpatialExperiment}
#' @return
#' A list with elements \code{bad} and \code{good}, each a data.frame
#' with \code{qscore_train} and (for “good”) an \code{is_a_bad_boy} flag.
#' @keywords internal
.computeCosmxProteinTrainSet <- function(spe)
{
    spe <- computeSpatialOutlier(spe, "log2AspectRatio", method="both")
    spe <- computeSpatialOutlier(spe, "log2Ctrl_total_ratio", method="both")
    train_bad <- data.frame(colData(spe)) |>
        filter((log2AspectRatio_outlier_mc == "HIGH" & dist_border < 50) |
                (log2AspectRatio_outlier_mc == "LOW" & dist_border < 50) |
                log2CountArea_outlier_train == "LOW" |
                log2Ctrl_total_ratio_outlier_sc == "HIGH") |>
        mutate(qscore_train = 0)

    train_good <- data.frame(colData(spe)) |>
        filter((log2AspectRatio > quantile(log2AspectRatio, probs = 0.25) &
                log2AspectRatio < quantile(log2AspectRatio, probs = 0.75) &
                dist_border > 50) |
                (log2CountArea > quantile(log2CountArea, probs = 0.90) &
                log2CountArea < quantile(log2CountArea, probs = 0.99))) |>
        mutate(qscore_train=1, is_a_bad_boy=cell_id %in% train_bad$cell_id)
    return(list(bad=train_bad, good=train_good))
}


#' computeQScoreFlags
#' @name computeQScoreFlags
#' @rdname computeQScoreFlags
#' @description
#' Compute flagged cells based on a manually chosen threshold on quality score
#'
#' This function Compute flagged cells based on a manually chosen threshold on
#' quality score stored in `SpatialExperiment` object.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#' @param qsThreshold Numeric threshold or quantile for quality score. Default
#'   `0.5`.
#' @param useQSQuantiles Logical; if `TRUE`, treat `qsThreshold` as a
#'   percentile.
#'
#' @return The `SpatialExperiment` object with added filter flags in `colData`.
#'
#'
#' @importFrom SummarizedExperiment colData
#' @export
#' @examples
#' example(computeQScore)
#' spe <- computeQScoreFlags(spe)
#' table(spe$low_qscore)
#' # if fixed filters are defined we have an additional column
#' spe <- computeThresholdFlags(spe)
#' spe <- computeQScoreFlags(spe)
#' table(spe$low_threshold_qscore)
computeQScoreFlags <- function(spe, qsThreshold=0.5, useQSQuantiles=FALSE) {
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("quality_score" %in% names(colData(spe)))

    if(useQSQuantiles) {
        spe$low_qscore <- ifelse(
            spe$quality_score < quantile(spe$quality_score, probs=qsThreshold),
            TRUE, FALSE)
    } else {
        spe$low_qscore <- spe$quality_score < qsThreshold
    }

    if("threshold_flags" %in% names(colData(spe))) {
        spe$low_threshold_qscore <- (spe$low_qscore &
                                        spe$threshold_flags)
    }
    return(spe)
}
