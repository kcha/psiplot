#' Plot PSI values of a given alternative splicing event
#'
#' @param x A 1-row data frame containing PSI values to be plotted
#' @param config Optional configuration settings for \code{plot_event}. Can be
#' a path to the \code{.config} file, or 4-column data frame of the \code{.config} file. Use the latter option if you are calling \code{plot_event} multiple times.
#' @param errorbar Logical indicating whether error bars should be drawn
#' @param groupmean Logical indicating whether grouped means should be drawn.
#' Requires \code{config}.
#' @param col Vector of colors for the points on the plot
#' @param title Title of the plot
#' @param xlab The x-axis label
#' @param ylab The y-axis label
#' @param ylim
#' @param xlim
#' @param cex.main
#' @param cex.axis
#' @param cex.xaxis
#' @param pch
#' @param cex
#' @param gridlines Logical indicating whether grid lines should be drawn
#' @export
#' @examples
#' data(mm.psi)
#' head(mm.psi)
#' plot_event(mm.psi[1,]) 
#'
#' # Plot with custom configuration
#' data(mm.psi.config)
#' plot_event(mm.psi[1,], config = mm.psi.config)
plot_event <- function(
  x, config = NULL, errorbar = TRUE,
  groupmean = ifelse(is.null(config), FALSE, TRUE), col = NULL,
  title = NULL, xlab = "", ylab = "", ylim = c(1,100),
  xlim = c(1, ncol(x)/2), cex.main = 0.9, cex.axis = 0.8, cex.xaxis = 0.6,
  pch = 20, cex = 1, gridlines = TRUE) {
  # Add some checks
#   if (!is.null(config) && is.character(config)) {
#     config <- read.csv(config, sep="\t")
#   }

  # Format input
  x <- format_table(x)
  reordered <- preprocess_sample_colors(x, config)
  psi <- reordered$data

  # Set plot title
  if (is.null(title)) {
    title <- make_title(rownames(x))
  }

  # Set up plot
  plot(NA,
       main=title,
       ylab=ylab, xlab=xlab, xaxt="n",
       ylim=ylim, xlim=xlim,
       cex.main=cex.main, cex.axis=cex.axis)
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

  # Draw psi
  points(1:ncol(psi), as.numeric(psi), col=as.character(reordered$col),
         pch=pch, cex = cex)

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
#' f <- format_table(mm.psi)
#' print(make_title(rownames(f)[1]))
make_title <- function(x) {
  event <- strsplit(x, split = "\\|")[[1]]
  sprintf("%s (position = %s, length = %s, type = %s)",
                   event[2], event[3], event[4], event[1])
}
