#' Plot multiple events or genes as a heatmap
#'
#' Visualize multiple PSI or cRPKM values in a single plot.
#'
#' @details
#' Input is similar to \code{\link{plot_event}} or \code{\link{plot_expr}}.
#'
#' Uses \code{\link[ggplot2]{geom_tile}} to generate heatmap.
#'
#' @param df A data frame of input values (PSI or cRPKM). If the latter, need to
#' set \code{expr = TRUE}.
#' @param config Optional configuration settings. Can be
#' a path to the \code{.config} file, or 4-column data frame of the \code{.config}
#' file.
#' @param expr Logical - \code{TRUE} if plotting cRPKMs, \code{FALSE} otherwise
#' @param xlab The x-axis label
#' @param ylab The y-axis label
#' @param title Title of the plot
#' @param cluster_rows Logical to cluster rows using heirarchical clustering
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
#' @examples
#' \dontrun{
#' plot_multi(psi, config = config)
#' plot_multi(crpkm, config = config, expr = TRUE)
#' plot_multi(crpkm, config = config, expr = TRUE, cluster_rows)
#' }
plot_multi <- function(df, config = NULL, expr = FALSE, xlab = "", ylab = "",
           title = "Heatmap", cluster_rows = FALSE) {
  # Format input
  df <- format_table(df, expr = expr)
  reordered <- preprocess_sample_colors(df, config = config, expr = expr)

  if (cluster_rows) {
    # Perform heirarchical clustering of events/genes
    hr <- hclust(dist(reordered$data))
    reordered$data <- reordered$data[hr$order,]
  }

  reordered$data$id <- rownames(reordered$data)

  m <- melt(reordered$data)
  m$variable <- factor(m$variable, levels = colnames(reordered$data),
                       ordered = TRUE)
  m$id <- factor(m$id, levels = reordered$data$id, ordered = TRU)

  gp <- ggplot(m, aes(x = variable, y = id)) +
    geom_tile(aes(fill = value)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle= 45, hjust = 0.5, vjust = 0, size = 9),
          axis.text.y = element_text(size = 8)) +
    xlab(xlab) + ylab(ylab) + ggtitle(title) +
    coord_fixed(ratio = ratio)

  if (expr) {
    gp <- gp +
      scale_fill_gradientn(colours = RColorBrewer::brewer.pal(6, "YlOrRd"),
                          name = "Expr", na.value = "white")
  } else {
    gp <- gp +
      scale_fill_gradient2(low = "yellow2", mid = "white", high = "blue2",
                                        na.value = "white", midpoint = 50,
                           name = "PSI")
  }
  return(gp)
}
