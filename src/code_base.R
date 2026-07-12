# Helper functions that can be use in all the project

.tiff_pdf_plot <- function(
    PLOT,
    HEIGHT,
    WIDTH,
    DIRECTORY,
    FILENAME) {
    #' Generate TIFF and PDF plots
    #'
    #' @param PLOT ggplot2 object to save
    #' @param HEIGHT Plot height in inches
    #' @param WIDTH Plot width in inches
    #' @param DIRECTORY Path to the output directory
    #' @return TIFF and PDF files written to disk

    ggplot2::ggsave(
        filename = file.path(DIRECTORY, paste0(FILENAME, ".tiff")),
        plot = PLOT,
        height = HEIGHT,
        width = WIDTH,
        dpi = 300
    )

    ggplot2::ggsave(
        filename = file.path(DIRECTORY, paste0(FILENAME, ".pdf")),
        plot = PLOT,
        height = HEIGHT,
        width = WIDTH
    )
}

.ridge_plot <- function(
    SEURAT_OBJECT,
    MARKERS = NULL,
    PALETTE = NULL,
    COLUMNS = NULL,
    ASSAY = "RNA") {
    #' Recreate Ridge Plot from Seurat by integrating MetBrewer
    #'
    #' @param SEURAT_OBJECT Seurat object post QC, processing and marker identification
    #' @param MARKERS character vector of gene markers for expression level plotting
    #' @param PALETTE MetBrewer palette name
    #' @param COLUMNS Number of columns for the layout of plots
    #' @param ASSAY Assay to use, default `RNA` but `SCT` available as well
    #' @return patchwork grid of ridge plots

    N_CLUSTERS <- length(levels(Seurat::Idents(SEURAT_OBJECT)))
    PALETTE_MAX <- length(MetBrewer::met.brewer(PALETTE))

    COLORS <- setNames(
        if (N_CLUSTERS > PALETTE_MAX) {
            MetBrewer::met.brewer(PALETTE, n = N_CLUSTERS, type = "continuous")
        } else {
            MetBrewer::met.brewer(PALETTE, n = N_CLUSTERS, type = "discrete")
        },
        levels(Seurat::Idents(SEURAT_OBJECT))
    )

    purrr::map(MARKERS, ~ SCpubr::do_RidgePlot(
        sample = SEURAT_OBJECT,
        feature = .x,
        assay = ASSAY,
        colors.use = COLORS,
        legend.position = "bottom"
    )) %>%
        patchwork::wrap_plots(ncol = COLUMNS) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom")
}

.dim_umap_scpub <- function(
    SEURAT_OBJECT,
    POINT_SIZE = NULL,
    LABEL_SIZE = NULL,
    LABEL = FALSE,
    PALETTE = NULL) {
    #' Generate UMAP publication-ready plot
    #'
    #' @param SEURAT_OBJECT Seurat object post clustering and UMAP reduction
    #' @param POINT_SIZE Point size for cells
    #' @param LABEL_SIZE Cluster label size
    #' @param PALETTE MetBrewer palette name
    #' @return ggplot2 object

    CLUSTER_IDS <- levels(Seurat::Idents(SEURAT_OBJECT))
    N_CLUSTERS <- length(CLUSTER_IDS)
    PALETTE_MAX <- length(MetBrewer::met.brewer(PALETTE))

    COLORS <- setNames(
        if (N_CLUSTERS > PALETTE_MAX) {
            MetBrewer::met.brewer(PALETTE, n = N_CLUSTERS, type = "continuous")
        } else {
            MetBrewer::met.brewer(PALETTE, n = N_CLUSTERS, type = "discrete")
        },
        CLUSTER_IDS
    )

    SCpubr::do_DimPlot(
    sample = SEURAT_OBJECT,
    reduction = "umap",
    pt.size = POINT_SIZE,
    label = LABEL,
    label.box = FALSE,
    label.size = LABEL_SIZE,
    colors.use = COLORS
    )
}

.feature_plot_markers <- function(
    SEURAT_OBJECT,
    MARKERS = NULL,
    PALATTE = NULL,
    COLORS = NULL) {
    #' Feature plots for marker genes arranged in a grid
    #'
    #' @param SEURAT_OBJECT Seurat object post marker identification
    #' @param MARKERS character vector of marker genes to plot
    #' @param PALETTE MetBrewer palette name, defaults to "Redon"
    #' @return patchwork grid of feature plots

    # COLORS <- MetBrewer::met.brewer(PALETTE, type = "continuous")

    purrr::map(MARKERS, ~ SCpubr::do_FeaturePlot(
        sample = SEURAT_OBJECT,
        features = .x,
        sequential.palette = COLORS
    )) %>%
        patchwork::wrap_plots()
}
