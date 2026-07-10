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

