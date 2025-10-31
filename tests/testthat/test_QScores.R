library(testthat)
library(SpaceTrooper)

# load the example SpatialExperiment
spe0 <- example(readCosmxSPE)$value

test_that("QC functions are exported", {
    expect_true(exists("spatialPerCellQC",   mode = "function"))
    expect_true(exists("computeQScore",      mode = "function"))
    expect_true(exists("computeSpatialOutlier", mode = "function"))
    expect_true(exists("computeQScoreFlags",  mode = "function"))
    expect_true(exists("computeThresholdFlags",  mode = "function"))
})


test_that("spatialPerCellQC adds per‐cell metrics to colData", {
    spe <- spatialPerCellQC(spe0, micronConvFact = 0.15)
    expect_s4_class(spe, "SpatialExperiment")
    cd <- colData(spe)
    # at least these should now be in colData()
    required <- c("sum", "detected", "total", "control_sum", "target_sum",
                  "Area_um", "log2AspectRatio", "log2CountArea")
    expect_true(all(required %in% colnames(cd)))
})

test_that("computeQScore adds a flag_score between 0 and 1", {
    spe <- spatialPerCellQC(spe0)
    spe2 <- computeQScore(spe)
    cd2 <- colData(spe2)
    expect_true("quality_score" %in% colnames(cd2))
    fs <- cd2$quality_score
    expect_true(is.numeric(fs))
    expect_true(all(fs >= 0 & fs <= 1))
})


test_that("computeSpatialOutlier flags outliers for a chosen metric", {
    spe <- spatialPerCellQC(spe0)
    # use 'Area_um' and both methods
    out <- computeSpatialOutlier(spe, computeBy="Area_um", method="scuttle")
    cd_out <- colData(out)
    expect_true(is.logical(cd_out$Area_um_outlier_sc) ||
                    is(cd_out$Area_um_outlier_sc, "outlier.filter"))
})


test_that("computeQScoreFlags combines filters and returns filter_out", {
    spe <- spatialPerCellQC(spe0)
    spe <- computeQScore(spe)
    ff <- computeThresholdFlags(spe,
                            totalThreshold = 10,
                            ctrlTotRatioThreshold = 0.2)
    cd_ff <- colData(ff)
    # expect logical columns and combined filter_out
    expect_true(all(c("is_zero_counts", "is_ctrl_tot_outlier",
                      "threshold_flags") %in% colnames(cd_ff)))
    expect_true(is.logical(cd_ff$is_zero_counts))
    expect_true(is.logical(cd_ff$is_ctrl_tot_outlier))
    expect_true(is.logical(cd_ff$threshold_flags))
    # filter_out should only be TRUE where all three flags are TRUE
    combined <- (ff$is_zero_counts & ff$is_ctrl_tot_outlier)
    expect_identical(combined, cd_ff$threshold_flags)
})


# test_that("computeQScoreMods identical to computeQScore on CosMx", {
#     sample_dir <- system.file(file.path("extdata", "CosMx_DBKero_Tiny"),
#                                package="SpaceTrooper")
#     skip_if(sample_dir == "" || !dir.exists(sample_dir),
#             "no CosMx example data found in inst/extdata; salto il test")
#
#     spe <- readCosmxSPE(sample_dir)
#
#     spe <- spatialPerCellQC(spe)
#     library(withr)
#     with_seed(1998,
#         qs_orig <- computeQScore(spe)
#     )
#     with_seed(1998,
#         qs_mod  <- computeQScoreMods(spe)
#     )
#
#     expect_equal(qs_mod$quality_score, qs_orig$quality_score, tolerance = 1e-8)
# })
