#' Plot PSI values of a given alternative splicing event
#'
#' @param x A 1-row data frame containing PSI values to be plotted
#' @param config An optional psiplot config object
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
#' @return Nothing
plot_event <- function(
  x, config = NULL, errorbar = TRUE, groupmean = config, col = NULL,
  title = NULL, xlab = "", ylab = "", ylim = c(1,100),
  xlim = c(1, ncol(x)), cex.main = 0.9, cex.axis = 0.8, cex.xaxis = 0.6,
  pch = 20, gridlines = TRUE) {
  # Add some checks


  # Set plot title
  if (is.null(title)) {
    title <- make_title(rownames(x))
  }

  # Set up plot
  plot(NA,
       main=title,
       ylab=ylab, xlab=xlab, xaxt="n",
       ylim=ylim, xlim=xlim,
       cex.main=cex.main, cex.axis=x.axis)
  axis(1, at=seq(1, ncol(PSIs), by=1), labels = FALSE)
  text(seq(1, ncol(PSIs), by=1),
       par("usr")[3] - 3.5,
       labels = samples,
       srt = 45, adj=c(1,1), xpd = TRUE,cex=cex.xaxis)


  # Draw error bars
  if (errorbar) {
    ci <- get_beta_ci(reordered$qual[i,])
    ci[which(is.na(reordered$data[i,])),] <- NA

    arrows(1:ncol(PSIs), ci[,1],
           1:ncol(PSIs), ci[,2],
           length = 0.025,
           angle = 90,
           code = 3, col = as.character(supercolors))
  }

  # Draw horizontal lines for groups
  if (!is.null(config_file) && groupmean) {
    seen <- vector()
    groups <- names(reordered$group.index)
    for (t in 1:length(groups)) {
      abline(h=mean(PSIs[i, reordered$group.index[[groups[t]]] ],
                    na.rm=TRUE),
             col=reordered$group.col[groups[t]], lwd=0.5)
      seen <- append(seen, paste(groups[t]))
    }

    # plot legend for mean group values
    if (length(seen) > 0) {
      legend_position <- ifelse(PSIs[i,ncol(PSIs)] > 50, "bottomright", "topright")
      legend(legend_position, legend = seen, lty = 1, col = reordered$group.col,
             title = "Group Means", cex = 0.7, ncol = 2)
    }
  }

  # Draw grid lines
  if (gridlines) {
    abline(v=1:ncol(PSIs), col="grey", lwd=0.3, lty=2)
    abline(h=seq(0,100,10), col="grey", lwd=0.3, lty=2)
  }

  # Draw PSIs
  points(1:ncol(PSIs), as.numeric(PSIs[i,]), col=as.character(supercolors),
         pch=pch, cex = cex)

}

#' Make plot title
#'
#' Create a plot title using the event ID
#'
#' @param x A character containing the event ID like
#' \code{S|TSPAN6|chrX:99885756-99885863|108}
#' @return A character with a human-friendly title
make_title <- function(x) {
  event <- strsplit(x, split = "\\|")[[1]]
  sprintf("%s (position = %s, length = %s, type = %s)",
                   event[2], event[3], event[4], event[1])
}
