#' spatialPerCellQC
#' @description
#' Perform Per-Cell Quality Control on a SpatialExperiment Object
#'
#' This function calculates quality control metrics for each cell in a
#' `SpatialExperiment` object and adds them to `colData`.
#'
#' @param spe A `SpatialExperiment` object containing spatial transcriptomics
#' data.
#' @param micronConvFact A numeric value for converting pixel dimensions to
#' microns. Default is `0.12`.
#' @param negProbList A character vector of patterns used to identify negative
#' probes.
#' Default values are:
#' - For Nanostring CosMx: `"NegPrb"`, `"Negative"`, `"SystemControl"`
#' - For Xenium: `"NegControlProbe"`, `"NegControlCodeWord"`,
#' `"UnassignedCodeWord"`
#' - For MERFISH: `"Blank"`
#'
#'
#' @return The `SpatialExperiment` object with added quality control metrics in
#' `colData`.
#'
#' @details The function computes several QC metrics, including control probe
#' sums, target probe sums, and the ratio of control probes to the total.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom scater addPerCellQC
#' @importFrom S4Vectors cbind.DataFrame
#' @export
#' @example
#' # TBD
spatialPerCellQC <- function(spe, micronConvFact=0.12,
                             negProbList=c("NegPrb", "Negative", "SystemControl", # CosMx
                                           "Ms IgG1", "Rb IgG", # CosMx Protein
                                           "NegControlProbe", "NegControlCodeWord", "UnassignedCodeWord", # Xenium
                                           "Blank" # MERFISH
                             ))
{
    ## CHECK EXISTENCE OF POLYGONS/AREAS ETC -> create function for metrics creation
    stopifnot(is(object=spe, "SpatialExperiment"))
    idxlist <- lapply(negProbList, function(ng){
        grep(paste0("^", ng), rownames(spe))
    })
    names(idxlist) <- negProbList
    idxlist <- idxlist[which(lengths(idxlist)!=0)] ## lengths instead of length

    spe <- addPerCellQC(spe, subsets=idxlist)
    idx <- grep("^subsets_.*_sum$", colnames(colData(spe)))

    if ( length(idx) !=0 )
    {
        ## TO TEST -> BENEDETTA
        npc <- rowSums(as.matrix(colData(spe)[ , idx, drop=FALSE])) # meglio dataframe, perché rowsums non funziona
        ## getting detected probes as the column suddenly after the sum column
        npd <- rowSums(as.matrix(colData(spe)[ , idx+1, drop=FALSE])) # not robust at all!
    }

    spe$control_sum <- npc
    spe$control_detected <- npd
    spe$target_sum <- spe$sum - npc
    spe$target_detected <- spe$detected - npd

    #### CHANGE SPE constructor WITH COORDINATES IN COLDATA #########
    colData(spe) <- cbind.DataFrame(colData(spe), spatialCoords(spe))

    if(metadata(spe)$technology == "Nanostring_CosMx")
    {
        spnc <- spatialCoords(spe) * micronConvFact
        colnames(spnc) <- gsub("px", "um", spatialCoordsNames(spe))
        colData(spe) <- cbind.DataFrame(colData(spe), spnc)
        spe$Area_um <- spe$Area * (micronConvFact^2)
    }

    #### compute AspectRatio for other technologies ####
    if ("AspectRatio" %in% colnames(colData(spe)))
    {
        spe$log2AspectRatio <- log2(spe$AspectRatio)
    } else {
        warning(paste0("Missing aspect ratio in colData...\n",
                       "NB: This could lead to additional warnings or errors.\n",
                       "Missing AspectRatio can be computed by loading polygons."))
    }

    spe$ctrl_total_ratio <- spe$control_sum/spe$total
    spe$ctrl_total_ratio[which(is.na(spe$ctrl_total_ratio))] <- 0
    spe$log2CountArea <- log2(spe$sum/spe$Area_um)

    message("Removing ", dim(spe[,spe$sum==0])[2], " cells, they have 0 counts")
    spe <- spe[,!spe$sum==0]

    return(spe)
}

#' computeBorderDistanceCosMx
#' @description
#' Compute Distance to FoV Border in SpatialExperiment for CosMx technology.
#'
#' Calculates the minimum distance of each coordinate in a `SpatialExperiment`
#' object to the nearest border of the field of view (FOV) and adds it to
#' `colData`.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#' @param xwindim Width of the FOV in the x-dimension. Defaults to
#' `metadata(spe)$fov_dim[["xdim"]]`.
#' @param ywindim Height of the FOV in the y-dimension. Defaults to
#' `metadata(spe)$fov_dim[["ydim"]]`.
#'
#' @return The `SpatialExperiment` object with added border distance data in
#' `colData`.
#'
#' @importFrom dplyr left_join
#' @importFrom SummarizedExperiment colData
#' @export
#' @example
#' # TBD
computeBorderDistanceCosMx <- function(spe,
                                       xwindim=metadata(spe)$fov_dim[["xdim"]],
                                       ywindim=metadata(spe)$fov_dim[["ydim"]])
{
    stopifnot(is(spe, "SpatialExperiment"))

    cd <- colData(spe)
    cdf <- left_join(as.data.frame(cd), metadata(spe)$fov_positions, by="fov")
    spcn <- spatialCoordsNames(spe)
    fovpn <- colnames(metadata(spe)$fov_positions)[colnames(
        metadata(spe)$fov_positions)%in%c("x_global_px", "y_global_px")]

    cd$dist_border_x <- pmin(cdf[,spcn[1]] - cdf[,fovpn[1]],
                             (cdf[,fovpn[1]] + xwindim) - cdf[,spcn[1]])

    cd$dist_border_y <- pmin(cdf[,spcn[2]] - cdf[,fovpn[2]],
                             (cdf[,fovpn[2]] + ywindim) - cdf[,spcn[2]])

    cd$dist_border <- pmin(cd$dist_border_x, cd$dist_border_y)
    colData(spe) <- cd
    return(spe)
}

#' computeQCScore
#' @description
#' Compute for each cell a quality score based on the target counts, area in micron
#' and log2 of the aspect ratio.
#'
#' @param spe a SpatialExperiment object with target_counts, area in micron
#' and log2 of the aspect ratio in the `colData`, tipically computed with
#' \link{spatialPerCellQC}
#' @param custom a boolean to set if you are using custom polygons
#'
#' @return a `SpatialExperiment` object with a `flag_score` stored in
#' the `colData`.
#' @export
#'
#' @examples
#' # TBD
computeQCScore <- function(spe, a=1, b=1, custom = FALSE)#0.3, b=0.8)
{
    stopifnot(is(spe, "SpatialExperiment"))
    cd <- colData(spe)

    if(!all(c("target_sum", "log2CountArea", "log2AspectRatio") %in% colnames(cd)))
    {
        stop("One of target_sum, Area_um, log2AspectRatio is missing, ",
             "did you run spatialPerCellQC?")
    }
    # spe$log2CountArea <- log2(spe$target_sum/spe$Area_um)
    if (metadata(spe)$technology=="Nanostring_CosMx")
    {
        spe <- computeBorderDistanceCosMx(spe)

        if(custom == TRUE)
        {
            if(!all(c("target_sum", "cust_log2CountArea",
                      "cust_log2AspectRatio") %in% colnames(cd)))
            {
                stop("One of target_sum, cust_log2CountArea, cust_log2AspectRatio
                is missing, did you run custom customPolyMetrics?")
            }
            qs <- 1 / (1 + exp(-a * spe$cust_log2CountArea +
                                   b * abs(spe$cust_log2AspectRatio) *
                                   as.numeric(spe$dist_border < 50)))
        } else {

            qs <- 1 / (1 + exp(-a * spe$log2CountArea +
                                   b * abs(spe$log2AspectRatio) *
                                   as.numeric(spe$dist_border < 50)))
        }
    } else {
        qs <- 1 / (1 + exp(-a * spe$log2CountArea +
                               b * abs(spe$log2AspectRatio)))
    }

    spe$quality_score <- qs
    # spe$quality_score <- 1/(1 + exp(-a*spe$target_sum/spe$Area_um +
    #b*abs(spe$log2AspectRatio)))
    return(spe)
}



#' computeSpatialOutlier
#' @description
#' Computes outliers based on the Area (in micron) of the experiment.
#' It gives the possibility to choose between the medcouple (mc method argument)
#' and the MADs (scuttle method argument).
#'
#' @details
#' The medcouple method is a measure for the skeweness of univariate distribution
#' as described in Hubert M. et al. (2008).
#' In particular, the computed medcouple value must be in a range between -0.6
#' and 0.6 to computed adjusted boxplots and perform the outlier detection.
#' For median absolute deviations (MADs) method we just wrap the isOutlier
#' function in the scuttle package. Please see McCarthy DJ et al (2017)
#' for further details.
#'
#' @param spe a SpatialExperiment object with target_counts, area in micron
#' and log2 of the aspect ratio in the `colData`.
#' @param compute_by character indicating a `colData` column name on which compute
#' the outlier.
#' @param method one of `mc`, `scuttle`, `both`.
#' Use `mc` for medcouple, `scuttle` for median absolute deviations as computed
#' in `scuttle`, `both` for computing both of them.
#' @param mcDoScale logical indicating if the values to compute the medcouple
#' for the outlier detection should be scaled (default is FALSE, as suggested
#' by the original Medcouple authors.). See \link[robustbase]{mc} for further
#' readings.
#'
#' @return a SpatialExperiment object with additional column(s) (named as
#' the column name indicated in `column_by` followed by the outlier_sc/mc
#' nomenclature) with the outlier detection as `outlier.filter` logical class
#' object. This allows to store the thresholds as attributes of the column.
#' use attr(,"thresholds") to retrieve them.
#' @export
#' @importFrom robustbase mc adjbox
#' @importFrom e1071 skewness
#' @importFrom scuttle isOutlier outlier.filter
#'
#' @examples
#' # TBD
computeSpatialOutlier <- function(spe, compute_by=NULL,
                                method=c("mc", "scuttle", "both"),
                                mcDoScale=FALSE,
                                scuttleType=c("both", "lower", "higher"))
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot(!is.null(compute_by))
    stopifnot(compute_by %in% names(colData(spe)))

    options(mc_doScale_quiet=TRUE)

    method <- match.arg(method)
    scuttleType <- match.arg(scuttleType)
    cd <- colData(spe)
    cdcol <- cd[[compute_by]]
    mcfl=scuttlefl=FALSE
    switch(method,
           both={ mcfl=scuttlefl=TRUE },
           mc={ mcfl=TRUE },
           scuttle={ scuttlefl=TRUE },
           {stop("Method is not one of allowed methods")}
    )

    if (mcfl)
    {
        skw <- e1071::skewness(cdcol)
        if (skw<=0) warning("The distribution is symmetric. ",
                    "The medcouple is not suited for symmetric distributions. ",
                    "In your case we suggest to use the scater method instead.")

        mcval <- robustbase::mc(cdcol, doScale=mcDoScale)
        if ( any( (mcval <= -0.6), (mcval >= 0.6) ) )
            stop("Obtained medcouple value is: ", round(mcval, digits=4),
                 "\nIt doesn't meet the needed requirements for outlier",
                 " identification with this method.")
        names(cdcol) <- colnames(spe)
        outl <- robustbase::adjbox(cdcol, plot=FALSE)
        outsmc <- rep(FALSE, dim(cd)[1])
        outsmc[rownames(cd) %in% names(outl$out)] <- TRUE
        ## using scuttle outlier.filter class to store fences for each filtering
        ## in the attributes of the class
        outlier_mc <- scuttle::outlier.filter(outsmc)
        thrs <- as.numeric(outl$fence)
        names(thrs) <- c("lower", "higher")
        attr(outlier_mc, "thresholds") <- thrs
        cd$outlier_mc <- outlier_mc
        names(cd)[names(cd)=="outlier_mc"] <- paste0(compute_by, "_outlier_mc")
        # metadata(spe)$outlier_fences[[compute_by]] <- outl$fence ## DEPRECATED

        ## TODO: compute distributions in the adjusted boxplots to store them
        ## in the coldata
    }

    if (scuttlefl)
    {
        outssc <- scuttle::isOutlier(cdcol, type=scuttleType)
        cd$outlier_sc <- scuttle::outlier.filter(outssc)
        names(cd)[names(cd)=="outlier_sc"] <- paste0(compute_by, "_outlier_sc")
    }

    colData(spe) <- cd
    return(spe)
}

#' computeFixedFlags
#' @description
#' Compute Flagged cells using fixed thresholds for SpatialExperiment
#'
#' This function calculates various flags to identify outliers in a
#' `SpatialExperiment` object based on quality control metrics.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#' @param total_threshold A numeric value for the threshold of total counts to
#' identify cells with zero counts. Default is `0`.
#' @param ctrl_tot_ratio_threshold A numeric value for the threshold of
#' control-to-total ratio to flag outliers. Default is `0.1`.
#'
#' @return The `SpatialExperiment` object with added filter flags in `colData`.
#'
#' @details The function flags cells based on zero counts, control-to-total
#' ratio, and `quality_score` to identify potential outliers. It also combines
#' these flags into a single filter flag.
#'
#' @importFrom SummarizedExperiment colData
#' @export computeFixedFlag
#' @examples
#' # TBD
computeFixedFlags <- function(spe,
                              total_threshold=0,
                              ctrl_tot_ratio_threshold=0.1)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("total" %in% names(colData(spe)))
    stopifnot("ctrl_total_ratio" %in% names(colData(spe)))

    spe$is_zero_counts <- ifelse(spe$total == total_threshold, TRUE, FALSE)
    #flagging cells with probe counts on total counts ratio > 0.1
    spe$is_ctrl_tot_outlier <- ifelse(spe$ctrl_total_ratio >
                                          ctrl_tot_ratio_threshold, TRUE, FALSE)

    spe$fixed_filter_out <- (spe$is_ctrl_tot_outlier &
                                 spe$is_zero_counts)
    return(spe)
}

#' qscoreOpt
#' @description
#' Compute automatically weights for quality score through glm training
#'
#' This function optimizes quality score by computing a glm to determine weights
#' for the quality score formula after selecting a training set from the
#' `SpatialExperiment` object.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#' @param plot a boolean to set to TRUE if you want to see the glm plot
#' @param custom a boolean value to set to TRUE if you are using custom polygons
#'
#' @return The `SpatialExperiment` object with added filter flags in `colData`.
#'
#' @details
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr filter mutate sample_n full_join
#' @importFrom ggplot2 ggplot aes geom_point stat_smooth geom_hline theme_bw
#' @importFrom ggplot2 xlab ylab
#' @export qscoreOpt
#' @examples
#' # TBD

qscoreOpt <- function(spe = spe, plot = FALSE, custom = FALSE){

    if(custom==TRUE){
        train_df1 <- data.frame(colData(spe)) |> dplyr::filter (is_ctrl_tot_outlier |
                     cust_Area_um_outlier_mc | Mean.DAPI_outlier_mc) |>
            dplyr::mutate(qscore_train = 0)

        message(paste0("Chosen low quality examples: ", dim(train_df1)[1]))

        train_df2 <- data.frame(colData(spe)) |> dplyr::filter(
            (cust_log2AspectRatio > summary(cust_log2AspectRatio)[2] &
                 cust_log2AspectRatio < summary(cust_log2AspectRatio)[5]) |
                (cust_log2CountArea > summary(cust_log2CountArea)[2] &
                 cust_log2CountArea >summary(cust_log2CountArea)[5]) |
                dist_border > 50) |>
            dplyr::mutate(qscore_train = 1) |> dplyr::sample_n(dim(train_df1)[1])

        train_df <- train_df1 %>%
            dplyr::mutate(rn = data.table::rowid(cell_id)) %>%
            dplyr::full_join(train_df2 %>%
                          dplyr::mutate(rn = data.table::rowid(cell_id))) %>%
            dplyr::select(-rn)

        model <- glm(qscore_train ~ cust_log2CountArea +
                         I(abs(cust_log2AspectRatio) * as.numeric(dist_border<50)),
                     family = binomial(link = "logit"), data = train_df)

        quality_score_opt <- 1 / (1 + exp(-model$coef[1] - model$coef[2] *
                                              spe$cust_log2CountArea - model$coef[3] *
                                              abs(spe$cust_log2AspectRatio) *
                                              as.numeric(spe$dist_border<50)))

        message(paste0("Computing optimized quality score with coefficients: ",
                       round(model$coef[1],2), ", ", round(model$coef[2],2),
                       ", ", round(model$coef[3],2)))

        if(any(is.na(quality_score_opt))){
            message(paste0(
                table(
                    is.na(quality_score_opt))[["TRUE"]]," NA quality scores found.
      Subsetting spe to keep only cells not NA for quality score"))
            spe$quality_score_opt <- quality_score_opt
            spe <- spe[,!is.na(spe$quality_score_opt)]
        }else{
            spe$quality_score_opt <- quality_score_opt
        }

        if(plot == TRUE){
            print(ggplot2::ggplot(train_df, ggplot2::aes(x = cust_log2CountArea +
                                           abs(cust_log2AspectRatio) *
                                           as.numeric(dist_border<50),
                                       y = qscore_train)) +
                      ggplot2::geom_point(col = "blue4") +
                      ggplot2::stat_smooth(method = "glm", se = FALSE, method.args =
                                      list(family=binomial), color = "pink2") +
                      ggplot2::geom_hline(yintercept = c(0, 1), color = "lightblue") +
                      ggplot2::theme_bw() + ggplot2::ylab("Training quality score") +
                      ggplot2::xlab("log2(count/area) + |aspectRatio| * distance from border index"))
            return(spe)
        }
    }else{
        train_df1 <- data.frame(colData(spe)) |> dplyr::filter (is_ctrl_tot_outlier |
                     Area_um_outlier_mc | Mean.DAPI_outlier_mc) |>
            dplyr::mutate(qscore_train = 0)

        message(paste0("Chosen low quality examples: ", dim(train_df1)[1]))

        train_df2 <- data.frame(colData(spe)) |> dplyr::filter(
            (log2AspectRatio > summary(log2AspectRatio)[2] & log2AspectRatio <
                 summary(log2AspectRatio)[5]) |
                (log2CountArea > summary(log2CountArea)[2] & log2CountArea >
                     summary(log2CountArea)[5]) |
                dist_border > 50) |>
            dplyr::mutate(qscore_train = 1) |> dplyr::sample_n(dim(train_df1)[1])

        train_df <- train_df1 %>%
            dplyr::mutate(rn = data.table::rowid(cell_id)) %>%
            dplyr::full_join(train_df2 %>%
                          dplyr::mutate(rn = data.table::rowid(cell_id))) %>%
            dplyr::select(-rn)

        model <- glm(qscore_train ~ log2CountArea + I(abs(log2AspectRatio) *
                     as.numeric(dist_border<50)),
                     family = binomial(link = "logit"), data = train_df)

        quality_score_opt <- 1 / (1 + exp(-model$coef[1] - model$coef[2] *
                                              spe$log2CountArea - model$coef[3] *
                                              abs(spe$log2AspectRatio) *
                                              as.numeric(spe$dist_border<50)))

        message(paste0("Computing optimized quality score with coefficients: ",
                       round(model$coef[1],2), ", ", round(model$coef[2],2),
                       ", ", round(model$coef[3],2)))

        if(any(is.na(quality_score_opt))){
            message(paste0(
                table(
                    is.na(quality_score_opt))[["TRUE"]]," NA quality scores found.
      Subsetting spe to keep only cells not NA for quality score"))
            spe$quality_score_opt <- quality_score_opt
            spe <- spe[,!is.na(spe$quality_score_opt)]
        }else{
            spe$quality_score_opt <- quality_score_opt
        }

        if(plot == TRUE){
            print(ggplot2::ggplot(train_df, ggplot2::aes(x = log2CountArea + abs(log2AspectRatio) *
                  as.numeric(dist_border<50), y = qscore_train)) +
                      ggplot2::geom_point(col = "blue4") +
                      ggplot2::stat_smooth(method = "glm", se = FALSE,
                                  method.args = list(family=binomial),
                                  color = "pink2") +
                      ggplot2::geom_hline(yintercept = c(0, 1), color = "lightblue") +
                      ggplot2::theme_bw() + ggplot2::ylab("Training quality score") +
                      ggplot2::xlab("log2(count/area) + |aspectRatio| * distance from border index"))
            return(spe)
        }
    }
    return(spe)
}

#' computeQscoreFlags
#' @description
#' Compute flagged cells based on a manually chosen threshold on quality score
#'
#' This function Compute flagged cells based on a manually chosen threshold on
#' quality score stored in `SpatialExperiment` object.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#' @param use_qs_quantiles a boolean. If TRUE uses the value specified in
#' qs_threshold as a percentile
#' @param opt a boolean value to set to TRUE if you want to compute flagged cells
#' for optimized quality score
#'
#' @return The `SpatialExperiment` object with added filter flags in `colData`.
#'
#' @details
#'
#' @importFrom SummarizedExperiment colData
#' @export computeQscoreFlags
#' @examples
#' # TBD
#'
computeQscoreFlags <- function(spe, qs_threshold=0.5,
                               use_qs_quantiles=FALSE,
                               opt = FALSE)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("quality_score" %in% names(colData(spe)))

    if(use_qs_quantiles)
    {
        if(opt==TRUE){
            spe$is_qscore_opt_outlier <- ifelse(spe$quality_score_opt <
                                                    quantile(spe$quality_score_opt,
                                                             probs=qs_threshold),
                                                TRUE, FALSE)
        }else{
            spe$is_qscore_outlier <- ifelse(spe$quality_score <
                                                quantile(spe$quality_score,
                                                         probs=qs_threshold),
                                            TRUE, FALSE)
        }

    } else {
        if(opt==TRUE){
            spe$is_qscore_opt_outlier <- spe$quality_score_opt < qs_threshold
        } else {
            spe$is_qscore_outlier <- spe$quality_score < qs_threshold
        }
    }
    if(opt==TRUE){
        spe$fixed_qscore_opt_out <- (spe$is_qscore_opt_outlier & spe$fixed_filter_out)
    } else {
        spe$fixed_qscore_out <- (spe$is_qscore_outlier & spe$fixed_filter_out)
    }
    return(spe)
}
