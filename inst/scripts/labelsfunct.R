addLabelsDBKero <- function(spe)
{
    labf <- system.file(file.path("extdata", "CosMx_DBKero_Tiny",
                            "labels_tiny.tsv"), package="SpaceTrooper")
    labs <-  read.table(file=labf, sep="\t", header=TRUE)
    spe$labels <- NA
    spe$labels_colors <- "black"
    spe$labels[match(labs$cell_id, spe$cell_id)] <- labs$label
    spe$labels_colors[match(labs$cell_id, spe$cell_id)] <- labs$lab_color
    spe$labels <- as.factor(spe$labels)
    spe$labels_colors <- as.factor(spe$labels_colors)
    spe
}
