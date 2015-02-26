#' Plot cRPKM values of a given gene
#'
#' Generate a plot of cRPKM values for a given gene. The cRPKM values are
#' obtained from the
#' \href{https://www.github.com/vastgroup/vast-tools}{vast-tools} pipeline.
#'
#' @details
#' Like \code{\link{plot_event}}, plots can be customized via the \code{config}
#' option. Either a data frame or
#' the filepath to the config file can be used. Alternatively, plots can be
#' customized using a limited set of graphical parameters as described above.
#'
#' See Details of \code{\link{plot_event}} for more information on costumizing
#' plots.
#'
#' cRPKM values that have \emph{NA} value are omitted and not plotted.
#'
#' @param x A 1-row data frame containing cRPKM values to be plotted
#' @param config Optional configuration settings for \code{plot_expr}. Can be
#' a path to the \code{.config} file, or 4-column data frame of the \code{.config}
#' file. Use the latter option if you are calling \code{plot_expr} multiple times.
#' @param groupmean Logical indicating whether grouped means should be drawn.
#' Requires \code{config}.
#' @param col Vector of colors with length matching the number of samples. If
#' specificed, this will override the color settings specified in \code{config}.
#' @param title Title of the plot
#' @param xlab The x-axis label
#' @param ylab The y-axis label
#' @param ylim Range of y-axis
#' @param cex.main Plot title size
#' @param cex.yaxis Y-axis font size
#' @param cex.xaxis X-axis font size (i.e. the sample names)
#' @param pch Point symbol
#' @param cex.pch Size of datapoints
#' @param lines Draw a connecting line between points for an event
#' @param gridlines Logical indicating whether grid lines should be drawn
#' @seealso
#' \code{\link{format_table}} for performing some initial conversion steps of \code{x}
#'
#' \code{\link{preprocess_sample_colors}} for pre-processing of \code{x} using
#' \code{config}
#'
#' @export
#' @examples
#' \dontrun{
#' plot_expr(crpkm[1,])
#'
#' # Plot with custom configuration
#' config
#' plot_expr(crpkm[1,], config = config, groupmean=TRUE)
#' plot_expr(crpkm[1,], config = "/path/to/config")
#'
#' # Plot using custom configuration, changing point sympol, and y-axis
#' # scale
#' plot_expr(crpkm[1,], config = config, pch = 9, ylim = c(20, 80))
#' }
plot_expr <- function(
  x, config = NULL,
  groupmean = ifelse(is.null(config), FALSE, TRUE), col = NULL,
  title = NULL, xlab = "", ylab = "Expression", ylim = NULL,
  cex.main = 0.9, cex.yaxis = 0.8, cex.xaxis = 0.6,
  pch = 20, cex.pch = 1, lines = FALSE, gridlines = TRUE) {
  if (nrow(x) != 1) {
    stop("Too many rows!")
  }

  N <- ncol(x) - 2
  if (N < 2) {
    stop("Need two or more samples!")
  }

  x <- format_table(x, expr = TRUE)
  reordered <- preprocess_sample_colors(x, config = config, expr = TRUE, col = col)
  crpkm <- reordered$data

  if (all(is.na(crpkm))) {
    warning("Did not find any points to plot")
  }

  # Set plot title
  if (is.null(title)) {
    title <- rownames(x)
  }

  if (is.null(ylim)) {
    ylim <- c(0, round(max(x)) + 1)
  }

  # Set up plot
  plot(NA,
       main=title,
       ylim = ylim, xlim = c(1, N),
       ylab=ylab, xlab=xlab, xaxt="n",
       cex.main=cex.main, cex.axis=cex.yaxis, las = 1)
  axis(1, at=seq(1, N, by=1), labels = FALSE)
  text(seq(1, N, by=1),
       par("usr")[3] - 0.2,
       labels = colnames(crpkm),
       srt = 45, adj=c(1,1), xpd = TRUE, cex=cex.xaxis)


  # Draw horizontal lines for groups
  if (!is.null(config) && groupmean) {
    draw_group_means(reordered)
  }

  # Draw grid lines
  if (gridlines) {
    abline(v=1:N, col="grey", lwd=0.5, lty=2)
    abline(h=seq(0,ylim[2],5), col="grey", lwd=0.5, lty=2)
  }

  # Draw line
  if (lines) {
    lines(1:N, as.numeric(crpkm), col = "black")
  }

  # Draw crpkm
  points(1:N, as.numeric(crpkm), col=as.character(reordered$col),
         pch=pch, cex = cex.pch)
}
