#' Draw group means as horizontal lines
#'
#' When this function is called, the group means are calculated and plotted
#' as horizontal lines. For internal use only.
#'
#' @param reordered a list containing input data and plot configuration. Produced
#' by \code{\link{preprocess_sample_colors}}
#' @seealso \code{\link{preprocess_sample_colors}}, \code{\link{plot_event}},
#' \code{\link{plot_expr}}
draw_group_means <- function(reordered) {
  # Draw horizontal lines for groups
  seen <- vector()
  groups <- names(reordered$group.index)
  for (t in 1:length(groups)) {
    mu <- mean(as.numeric(
      reordered$data[reordered$group.index[[groups[t]]]]
    ), na.rm=TRUE)
    abline(h=mu, col=reordered$group.col[groups[t]], lwd=0.5)
    seen <- append(seen, paste(groups[t]))
  }

  # plot legend for mean group values
  if (length(seen) > 0) {
    legend_position <- ifelse(reordered$data[ncol(reordered$data)] > 50,
                              "bottomright", "topright")
    legend(legend_position, legend = seen, lty = 1, col = reordered$group.col,
           title = "Group Means", cex = 0.6, ncol = 2)
  }
}
