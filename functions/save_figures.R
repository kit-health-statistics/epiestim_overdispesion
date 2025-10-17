#' Save a plot to the figures directory
#'
#' @param plot the plot object to save
#' @param filename string, the filename (without path and file extension)
#' @param width a positive real, the width in inches
#' @param height a positive real, the height in inches
#' @param dpi integer, resolution in dots per inch
#' @export
save_plot <- function(plot, filename, width = 14, height = 12.5, dpi = 400) {
  # Save as a PDF
  ggplot2::ggsave(
    filename = paste0(here::here("figure", filename), ".pdf"),
    plot = plot,
    width = width,
    height = height,
    create.dir = TRUE
  )

  # Save as a PNG to make the comparison in a PR on GitHub easier.
  # Can be deleted later.
  ggplot2::ggsave(
    filename = paste0(here::here("figure", filename), ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    create.dir = TRUE
  )
}
