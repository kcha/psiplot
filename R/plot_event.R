#' Plot PSI values of a given alternative splicing event
#'
#' Generate a plot of PSI values for a given exon. The PSI values and
#' corresponding quality scores are typically obtained from the
#' \href{https://www.github.com/vastgroup/vast-tools}{vast-tools}
#' pipeline.
#'
#' @details
#' Plots can be customized via the \code{config} option. Either a data frame or
#' the filepath to the config file can be used. Alternatively, plots can be
#' customized using a limited set of graphical parameters as described above.
#'
#' The order of samples (e.g. columns in \emph{x}) as it appears on the
#' resulting plot can be customized using a config file. If a config file is not
#' used and re-ordering is desired, then it must be done manually before calling
#' \code{\link{plot_event}} by altering the columns of \emph{x}.
#'
#' If groups are defined in \code{config} and \code{groupmean=TRUE}, the mean
#' PSI of the samples within each group are drawn as horizontal lines. The
#' colors of the lines are determined by the color set to the first sample of
#' each group by \code{RColorCode} in \code{config}. A corresponding legend key
#' will also be drawn.
#'
#' PSI values that have \emph{NA} value are omitted and not plotted.
#'
#' Error bars based on the confidence interval of the PSI estimation can be
#' shown by setting \code{errorbar=TRUE}.
#'
#' @param x A 1-row data frame containing PSI values to be plotted
#' @param config Optional configuration settings for \code{plot_event}. Can be
#' a path to the \code{.config} file, or 4-column data frame of the \code{.config}
#' file. Use the latter option if you are calling \code{plot_event} multiple times.
#' @param errorbar Logical indicating whether error bars should be drawn
#' @param groupmean Logical indicating whether grouped means should be drawn.
#' Requires \code{config}.
#' @param col Vector of colors with length matching the number of samples. If
#' specificed, this will override the color settings specified in \code{config}.
#' @param title Title of the plot
#' @param xlab The x-axis label
#' @param ylab The y-axis label
#' @param ylim Range of y-axis
#' @param xlim Range of x-axis
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
#' plot_event(psi[1,])
#'
#' # Plot with custom configuration
#' config
#' plot_event(psi[1,], config = config, groupmean=TRUE)
#' plot_event(psi[1,], config = "/path/to/config")
#'
#' # Plot using custom configuration, changing point sympol, and y-axis
#' # scale
#' plot_event(psi[1,], config = config, pch = 9, ylim = c(20, 80))
#' }
plot_event <- function(
  x, config = NULL, errorbar = TRUE,
  groupmean = ifelse(is.null(config), FALSE, TRUE), col = NULL,
  title = NULL, xlab = "", ylab = "PSI", ylim = c(1,100),
  xlim = c(1, ncol(x)/2), cex.main = 0.9, cex.yaxis = 0.8, cex.xaxis = 0.6,
  pch = 20, cex.pch = 1, lines = FALSE, gridlines = TRUE) {
  if (nrow(x) != 1) {
    stop("Too many rows!")
  }

  # Format input
  x <- format_table(x)
  reordered <- preprocess_sample_colors(x, config, col = col)
  psi <- reordered$data

  if (all(is.na(psi))) {
    warning("Did not find any points to plot")
  }

  # Set plot title
  if (is.null(title)) {
    title <- make_title(rownames(x))
  }

  # Set up plot
  plot(NA,
       main=title,
       ylab=ylab, xlab=xlab, xaxt="n",
       ylim=ylim, xlim=xlim,
       cex.main=cex.main, cex.axis=cex.yaxis)
  axis(1, at=seq(1, ncol(psi), by=1), labels = FALSE)
  text(seq(1, ncol(psi), by=1),
       par("usr")[3] - 3.5,
       labels = colnames(psi),
       srt = 45, adj=c(1,1), xpd = TRUE, cex=cex.xaxis)


  # Draw error bars
  if (errorbar) {
    ci <- get_beta_ci(reordered$qual)
    ci[which(is.na(psi)),] <- NA

    arrows(1:ncol(psi), ci[,1],
           1:ncol(psi), ci[,2],
           length = 0.025,
           angle = 90,
           code = 3, col = as.character(reordered$col))
  }

  # Draw horizontal lines for groups
  if (!is.null(config) && groupmean) {
    seen <- vector()
    groups <- names(reordered$group.index)
    for (t in 1:length(groups)) {
      mu <- mean(as.numeric(
        psi[reordered$group.index[[groups[t]]]]
        ), na.rm=TRUE)
      abline(h=mu, col=reordered$group.col[groups[t]], lwd=0.5)
      seen <- append(seen, paste(groups[t]))
    }

    # plot legend for mean group values
    if (length(seen) > 0) {
      legend_position <- ifelse(psi[ncol(psi)] > 50, "bottomright", "topright")
      legend(legend_position, legend = seen, lty = 1, col = reordered$group.col,
             title = "Group Means", cex = 0.7, ncol = 2)
    }
  }

  # Draw grid lines
  if (gridlines) {
    abline(v=1:ncol(psi), col="grey", lwd=0.3, lty=2)
    abline(h=seq(0,100,10), col="grey", lwd=0.3, lty=2)
  }

  # Draw line
  if (lines) {
    lines(1:ncol(psi), as.numeric(psi), col = "black")
  }

  # Draw psi
  points(1:ncol(psi), as.numeric(psi), col=as.character(reordered$col),
         pch=pch, cex = cex.pch)

}

#' Make plot title
#'
#' Create a plot title using the event ID
#'
#' @param x A character containing the event ID like
#' "\code{S|TSPAN6|chrX:99885756-99885863|108}"
#' @return A character with a human-friendly title
#' @export
#' @examples
#' f <- format_table(psi)
#' print(make_title(rownames(f)[1]))
make_title <- function(x) {
  event <- strsplit(x, split = "\\|")[[1]]
  sprintf("%s (position = %s, length = %s, type = %s)",
                   event[2], event[3], event[4], event[1])
}
