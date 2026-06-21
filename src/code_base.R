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
