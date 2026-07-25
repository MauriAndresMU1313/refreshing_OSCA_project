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

.get_palette <- function(
    SEURAT_OBJECT, 
    PALETTE) {
    #' Handle long list of colors by adapting from `continuous` to `dicrete` palette type
    #' 
    #' @param SEURAT_OBJECT that is necessary for plotting
    #' @param PALETTE to use from `MetBrewer`
    #' @return allows to extend the `n` for color palettes

    N <- length(levels(Seurat::Idents(SEURAT_OBJECT)))
    MetBrewer::met.brewer(
        PALETTE,
        n = N,
        type = if (N > length(MetBrewer::met.brewer(PALETTE))) "continuous" else "discrete"
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
    REDUCTION = NULL,
    GROUP_BY = NULL,
    SPLIT_BY = NULL,
    PALETTE = NULL) {
    #' Generate UMAP plot
    #'
    #' @param SEURAT_OBJECT Seurat object post clustering and UMAP reduction
    #' @param POINT_SIZE Point size for cells
    #' @param LABEL_SIZE Cluster label size
    #' @param LABEL Whether to show cluster labels
    #' @param REDUCTION Dimensional reduction to use
    #' @param GROUP_BY Metadata column to group cells by
    #' @param PALETTE MetBrewer palette name
    #' @return ggplot2 object

    # Handle metadata
    CLUSTER_IDS <- if (is.null(GROUP_BY)) {
        levels(Seurat::Idents(SEURAT_OBJECT))
    } else {
        unique(SEURAT_OBJECT@meta.data[[GROUP_BY]])
    }

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
        reduction = REDUCTION,
        pt.size = POINT_SIZE,
        split.by = SPLIT_BY,
        label = LABEL,
        label.box = FALSE,
        label.size = LABEL_SIZE,
        group.by = GROUP_BY,
        colors.use = COLORS
    )
}

.feature_plot_markers <- function(
    SEURAT_OBJECT,
    MARKERS = NULL,
    PALETTE = NULL,
    SPLIT_BY = NULL,
    MAX_CUTOFF = NULL,
    REDUCTION = "umap",
    COLUMNS = NULL) {
    #' Feature plots for marker genes arranged in a grid
    #'
    #' @param SEURAT_OBJECT Seurat object post marker identification
    #' @param MARKERS character vector of marker genes to plot
    #' @param PALETTE sequential palette name for continuous color scale
    #' @param SPLIT_BY metadata column to split plots by condition
    #' @param MAX_CUTOFF maximum expression cutoff value
    #' @param REDUCTION dimensional reduction to use, defaults to "umap"
    #' @param COLUMNS number of columns in patchwork layout
    #' @return patchwork grid of feature plots

    if (is.null(SPLIT_BY)) {
        purrr::map(MARKERS, ~ SCpubr::do_FeaturePlot(
            sample = SEURAT_OBJECT,
            features = .x,
            reduction = REDUCTION,
            max.cutoff = MAX_CUTOFF,
            sequential.palette = PALETTE,
            ncol = COLUMNS
        )) %>%
            patchwork::wrap_plots()
    } else {
        SCpubr::do_FeaturePlot(
            sample = SEURAT_OBJECT,
            features = MARKERS,
            split.by = SPLIT_BY,
            reduction = REDUCTION,
            max.cutoff = MAX_CUTOFF,
            sequential.palette = PALETTE,
            ncol = COLUMNS
        )
    }
}

.vln_plot_markers <- function(
    SEURAT_OBJECT,
    MARKERS = NULL,
    PALETTE = NULL,
    COLUMNS = NULL,
    ASSAY = NULL,
    SPLIT_BY = NULL,
    GROUP_BY = NULL) {
    #' Violin plots for marker genes arranged in a grid
    #'
    #' @param SEURAT_OBJECT post marker identification
    #' @param MARKERS character vector of marker
    #' @param PALETTE MetBrewer palette name
    #' @param COLUMNS Number of columns for the layout of plots
    #' @param ASSAY Assay to use
    #' @param SPLIT_BY metadata column to split plots by condition
    #' @param GROUP_BY metadata column to group cells by
    #' @return patchwork grid

    # Handle SPLIT_BY or GROUP_BY
    CLUSTER_IDS <- if (is.null(SPLIT_BY)) {
        if (is.null(GROUP_BY)) {
            levels(Seurat::Idents(SEURAT_OBJECT))
        } else {
            unique(SEURAT_OBJECT@meta.data[[GROUP_BY]])
        }
    } else {
        unique(SEURAT_OBJECT@meta.data[[SPLIT_BY]])
    }

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

    purrr::map(MARKERS, ~ SCpubr::do_ViolinPlot(
        sample = SEURAT_OBJECT,
        features = .x,
        assay = ASSAY,
        colors.use = COLORS,
        split.by = SPLIT_BY,
        group.by = GROUP_BY,
        order = is.null(SPLIT_BY),
        plot_boxplot = is.null(SPLIT_BY)
    )) %>%
        patchwork::wrap_plots(ncol = COLUMNS) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom")
}

.dot_plot_cons_markers <- function(
    SEURAT_OBJECT,
    MARKERS = NULL,
    SPLIT_BY = NULL,
    FLIP = FALSE,
    DOT_SCALE = NULL,
    TITLE = NULL) {
    #' Dot plot for conserved markers across conditions
    #'
    #' @param SEURAT_OBJECT Seurat object post marker identification
    #' @param MARKERS character vector or named list of marker genes
    #' @param SPLIT_BY metadata column to split plot by condition, defaults to NULL
    #' @param FLIP logical, flip axes, defaults to FALSE
    #' @param DOT_SCALE size scaling for dots, defaults to NULL
    #' @param TITLE plot title
    #' @return ggplot2 dot plot

    SCpubr::do_DotPlot(
        sample = SEURAT_OBJECT,
        features = MARKERS,
        split.by = SPLIT_BY,
        flip = FLIP,
        dot.scale = DOT_SCALE,
        plot.title = TITLE,
        axis.text.x.angle = 45,
        use_viridis = TRUE,
        font.size = 16
    )
}

.cell_scatter_comparison <- function(
    AGGREGATE_OBJECT,
    CELL_TYPES = NULL,
    CONDITION_1 = NULL,
    CONDITION_2 = NULL,
    MARKERS = NULL,
    COLUMNS = NULL,
    SIZE = NULL) {
    #' Cell scatter plots comparing two conditions per cell type
    #'
    #' @param AGGREGATE_OBJECT Aggregated Seurat object
    #' @param CELL_TYPES character vector of cell types to compare
    #' @param CONDITION_1
    #' @param CONDITION_2
    #' @param MARKERS character vector of genes to highlight and label
    #' @param COLUMNS number of columns in patchwork layout
    #' @return patchwork grid of labeled cell scatter plots

    purrr::map(CELL_TYPES, function(CELL_TYPE) {
        LABEL_1 <- paste0(CELL_TYPE, "_", CONDITION_1)
        LABEL_2 <- paste0(CELL_TYPE, "_", CONDITION_2)

        Seurat::CellScatter(
            AGGREGATE_OBJECT,
            cell1 = LABEL_1,
            cell2 = LABEL_2,
            highlight = MARKERS
        ) %>%
            Seurat::LabelPoints(
                points = MARKERS,
                repel = TRUE,
                size = SIZE
            )
    }) %>%
        patchwork::wrap_plots(ncol = COLUMNS)
}
