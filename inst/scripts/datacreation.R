outfolder <- "personal_folder"

outfolder <- "~/Library/CloudStorage/GoogleDrive-drighelli@gmail.com/My\ Drive/works/coding"
## creating data for cosmx starting from DBKero
cospath <- "~/Downloads/CosMx_data/DBKero/CosMx_Breast/CosMx_data_Case2"
# debug(readCosmxSPE)
spe<- readCosmxSPE(cospath)
f=c(16)
# c=c(1:50)
# gc<- 1:12

countmat_file <- list.files(cospath, "exprMat_file.csv", full.names=TRUE)
countmat <- data.table::fread(countmat_file, showProgress=FALSE) # cell count matrix
c1 <- countmat[countmat$fov%in%f,]
# c1 <- c1[c1$cell_ID%in%c,1:12]
write.csv(x=c1, file=file.path(outfolder, "/SpaceTrooper/inst/extdata/CosMx_DBKero_Tiny/DBKero_exprMat_file16.csv"), row.names=FALSE)

# spe32 <- spe[,spe$fov==32]
# spe32 <- spe32[1:22, 1:100]
# spe32 <- spatialPerCellQC(spe32)
# spe32 <- computeQScore(spe32)
# spe32$quality_score
# spe32 <- computeFixedFlags(spe32)
# FirstFilterPlot(spe32)
metadata_file <- list.files(cospath, "metadata_file.csv", full.names=TRUE)
metadata <- data.table::fread(metadata_file, showProgress=FALSE) # cell metadata
m1 <- metadata[metadata$fov==f,]
# m1 <- m1[m1$cell_ID%in%c,]
write.csv(x=m1, file=file.path(outfolder,"/SpaceTrooper/inst/extdata/CosMx_DBKero_Tiny/DBKero_metadata_file16.csv"), row.names=FALSE)

fovpos_file <- list.files(cospath, "fov_positions_file.csv", full.names=TRUE)
fov_positions <- as.data.frame(data.table::fread(fovpos_file, header=TRUE))
f1 <- fov_positions[fov_positions$fov==f,]
write.csv(x=f1, file=file.path(outfolder,"/SpaceTrooper/inst/extdata/CosMx_DBKero_Tiny/DBKero_fov_positions_file16.csv"), row.names=FALSE)

pol_file <- metadata(spe)$polygons

spat_obj <- data.table::fread(pol_file)
s1 <- spat_obj[spat_obj$fov==f, ]
write.csv(x=s1, file=file.path(outfolder,"/SpaceTrooper/inst/extdata/CosMx_DBKero_Tiny/DBKero-polygons16.csv"), row.names=FALSE)

## data creation for merscope

library(devtools)
load_all()
# folder <- "~/Downloads/MERSCOPE_data/Human_brain"
folder <- "~/Downloads/Merfish_data/human_uterine_cancer_patient2"

debug(readMerfishSPE)
spe <- readMerfishSPE(folder, compute_missing_metrics=FALSE, keep_polygons=FALSE)
countmat_file <- list.files(dirname, countmatfpattern, full.names=TRUE)
metadata_file <- list.files(dirname, metadatafpattern, full.names=TRUE)
# Read in
countmat <- data.table::fread(countmat_file)
c1 <- countmat[1:50,1:11]
write.csv(x=c1, file=file.path(outfolder, "/SpaceTrooper/inst/extdata/Merfish_Tiny/HumanUterineCancerPatient2_cell_by_gene.csv"), row.names=FALSE)

metadata <- data.table::fread(metadata_file) # cell metadata
m1 <- metadata[1:50,]
write.csv(x=m1, file=file.path(outfolder, "/SpaceTrooper/inst/extdata/Merfish_Tiny/HumanUterineCancerPatient2_cell_metadata.csv"), row.names=FALSE)

debug(readPolygonsMerfish)
pols <- readPolygonsMerfish(metadata(spe)$polygons, type="parquet")
mpols <- pols$cell_id[pols$is_multi]
polfile <- metadata(spe)$polygons
polygons <- arrow::read_parquet(polfile, as_data_frame=FALSE)
polygons$metadata
m1 <- data.table::fread(file.path(outfolder, "/SpaceTrooper/inst/extdata/Merfish_Tiny/HumanUterineCancerPatient2_cell_metadata.csv")) # cell metadata
library(dplyr)
polygons
p1 <- polygons %>%
    filter(EntityID %in% m1$EntityID)
p2 <- polygons[polygons$EntityID %in% m1$EntityID,]
p1$Geometry
arrow::write_parquet(p1, sink=file.path(outfolder, "/SpaceTrooper/inst/extdata/Merfish_Tiny/HumanUterineCancerPatient2_cell_boundaries.parquet"))
p3 <- arrow::read_parquet(file.path(outfolder, "/SpaceTrooper/inst/extdata/Merfish_Tiny/HumanUterineCancerPatient2_cell_boundaries.parquet"), as_data_frame=TRUE)
polfiles <- file.path(outfolder, "/SpaceTrooper/inst/extdata/Merfish_Tiny/HumanUterineCancerPatient2_cell_boundaries.parquet")
readPolygonsMerfish(polfiles, "parquet")


#### CosMx Protein
spe <- readCosmxProteinSPE("~/Downloads/CosMx_data/S0_prot")
# spe <- spe1
spe1 <- spe
# x_shift_px <- 1
y_shift_px <- 180
fov_px <- 4239.291

metadata(spe)$fov_positions <- metadata(spe)$fov_positions |>
    dplyr::mutate(
        # x_global_px = x_global_mm / 0.12028 * 1000 + x_shift_px,
        y_global_px = y_global_mm / 0.12028 * 1000 - fov_px + y_shift_px
    )
fovs <- c(16,17,28,29)
fovs <- c(30:100)
fovs <- c(59,60,71,72)
spe10<-spe[, spe$fov%in%fovs]
metadata(spe10)$fov_positions <- metadata(spe10)$fov_positions[metadata(spe10)$fov_positions$fov%in%fovs,]
plotCellsFovs(spe10)

fovs <- c(60)
cospath <- "~/Downloads/CosMx_data/S0_prot"
countmat_file <- list.files(cospath, "exprMat_file.csv", full.names=TRUE)
countmat <- data.table::fread(countmat_file, showProgress=FALSE) # cell count matrix
c1 <- countmat[countmat$fov%in%fovs,]
# c1 <- c1[c1$cell_ID%in%c,1:12]
write.csv(x=c1, file=file.path(outfolder, "/SpaceTrooper/inst/extdata/S0_prot/S0_exprMat_file.csv.gz"), row.names=FALSE)

metadata_file <- list.files(cospath, "metadata_file.csv", full.names=TRUE)
metadata <- data.table::fread(metadata_file, showProgress=FALSE) # cell metadata
m1 <- metadata[metadata$fov==fovs,]
write.csv(x=m1, file=file.path(outfolder,"/SpaceTrooper/inst/extdata/S0_prot/S0_metadata_file.csv"), row.names=FALSE)

fovpos_file <- list.files(cospath, "fov_positions_file.csv", full.names=TRUE)
fov_positions <- as.data.frame(data.table::fread(fovpos_file, header=TRUE))
f1 <- fov_positions[fov_positions$FOV==fovs,]
y_shift_px <- 180
fov_px <- 4239.291

f1 <- f1 |>
    dplyr::mutate(
        # x_global_px = x_global_mm / 0.12028 * 1000 + x_shift_px,
        y_global_px = y_global_mm / 0.12028 * 1000 - fov_px + y_shift_px
    )
write.csv(x=f1, file=file.path(outfolder,"/SpaceTrooper/inst/extdata/S0_prot/S0_fov_positions_file.csv"), row.names=FALSE)

pol_file <- metadata(spe)$polygons

spat_obj <- data.table::fread(pol_file)
s1 <- spat_obj[spat_obj$fov==fovs, ]
write.csv(x=s1, file=file.path(outfolder,"/SpaceTrooper/inst/extdata/S0_prot/S0-polygons.csv"), row.names=FALSE)

spe60 <- readCosmxProteinSPE("~/Library/CloudStorage/GoogleDrive-drighelli@gmail.com/My\ Drive/works/coding/SpaceTrooper/inst/extdata/S0_prot")
