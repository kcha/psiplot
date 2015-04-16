#' Draw group means as horizontal lines
#'
#' When this function is called, the group means are calculated and plotted
#' as horizontal lines. For internal use only.
#'
#' @param gp ggplot2 object
#' @param values a melted data frame of PSI/cRPKM values. Usually comes from \code{\link{plot_event}} or \code{\link{plot_expr}}.
#' @param cfg configuration as data frame
#' @param group.col group colours as computed from the \code{group.col} attribute returned by \code{\link{preprocess_sample_colors}}
#' @param reordered a list containing input data and plot configuration. Produced
#' by \code{\link{preprocess_sample_colors}}
#' @param offset an integer indicating the offset to print group label from line
#' @return ggplot2 object
#' @import plyr
#' @seealso \code{\link{preprocess_sample_colors}}, \code{\link{plot_event}},
#' \code{\link{plot_expr}}
draw_group_means <- function(gp, values, cfg, group.col, offset = 3) {
  m <- suppressMessages(join(values, cfg))
  msum <- ddply(m, .(GroupName), summarize, mu = mean(value, na.rm=TRUE))
  gp <- gp + geom_hline(data = msum, aes(yintercept = mu, colour = GroupName),
                  show_guide = TRUE) +
    scale_colour_manual("", values = group.col)
  return(gp)
}
