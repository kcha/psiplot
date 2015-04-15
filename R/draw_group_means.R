#' Draw group means as horizontal lines
#'
#' When this function is called, the group means are calculated and plotted
#' as horizontal lines. For internal use only.
#'
#' @param gp ggplot2 object
#' @param reordered a list containing input data and plot configuration. Produced
#' by \code{\link{preprocess_sample_colors}}
#' @param offset an integer indicating the offset to print group label from line
#' @return ggplot2 object
#' @seealso \code{\link{preprocess_sample_colors}}, \code{\link{plot_event}},
#' \code{\link{plot_expr}}
draw_group_means <- function(gp, reordered, offset = 2) {
  groups <- names(reordered$group.index)
  for (t in 1:length(groups)) {
    if (groups[t] %in% names(reordered$group.index)) {
      mu <- mean(as.numeric(
        reordered$data[reordered$group.index[[groups[t]]]]),
        na.rm=TRUE)
      gp <- gp + geom_hline(yintercept = mu,
                   colour = reordered$group.col[groups[t]]) +
        annotate("text", 1, min(100, mu + offset),
                 label = paste(groups[t], "Avg"),
                 size = 3,
                 color = reordered$group.col[groups[t]])

    }
  }
  return(gp)
}
