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
#' @param cex.main Plot title size (pts)
#' @param cex.yaxis Y-axis font size (pts)
#' @param cex.xaxis X-axis font size (i.e. the sample names) (pts)
#' @param pch Point symbol
#' @param cex.pch Size of datapoints
#' @param gridlines Logical indicating whether grid lines should be drawn
#' @param plot (deprecated) prints the plot
#' @return ggplot2 object
#' @seealso
#' \code{\link{format_table}} for performing some initial conversion steps of \code{x}
#'
#' \code{\link{preprocess_sample_colors}} for pre-processing of \code{x} using
#' \code{config}
#'
#' @export
#' @import methods
#' @import ggplot2
#' @importFrom reshape2 melt
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
  groupmean = ifelse(is.null(config), FALSE, TRUE), counts= FALSE, col = NULL,
  title = NULL, xlab = "", ylab = "Expression", ylim = NULL,
  cex.main = 14, cex.yaxis = 12, cex.xaxis = 12,
  pch = 20, cex.pch = 3, plot = NULL, gridlines = TRUE) {
  if (!missing(plot)) {
    warning("The option 'plot' has been deprecated")
  }

  if (nrow(x) != 1) {
    stop("Too many rows!")
  }

  N <- ncol(x) - 2
  if (N < 2) {
    stop("Need two or more samples!")
  }

  x <- format_table(x, expr = TRUE, counts=counts)
  reordered <- preprocess_sample_colors(x, config = config, expr = TRUE, col = col)
  crpkm <- reordered$data

  subg <- "SubgroupName" %in% colnames(reordered$config)

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

  mdata <- suppressMessages(melt(crpkm,
                                 measure.vars = names(crpkm),
                                 variable.name = "SampleName"))

  if(subg){
    sm <- suppressMessages(join(mdata,reordered$config))
    smsum <- ddply(sm,.(SubgroupName),summarize, value=mean(value,na.rm=T))
    smsum2 <- smsum[sapply(names(reordered$subgroup.order),
                           function(x) which(smsum$SubgroupName==x)),]
    smsum2$SubgroupName <- factor(smsum2$SubgroupName,smsum2$SubgroupName)
  }

  if(subg){
    gp <- ggplot(smsum2,aes(x=SubgroupName, y=value)) +
      geom_point(colour = reordered$subgroup.col, size = cex.pch, shape=pch)
  } else{
    gp <- ggplot(mdata, aes(x = SampleName, y = value)) +
      geom_point(colour = reordered$col, size=cex.pch, shape = pch)
  }

  gp <- gp +
    ylab("cRPKM") +
    xlab("") +
    ylim(ylim) +
    ggtitle(title) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = cex.xaxis),
          axis.text.y = element_text(size = cex.yaxis),
          axis.title.y = element_text(size = cex.yaxis),
          title = element_text(size = cex.main))
  


  # Draw horizontal lines for groups
  if (!is.null(config) && groupmean) {
    gp <- draw_group_means(gp, mdata, reordered$config, reordered$group.col)
  }
  
  if (!gridlines) {
    gp <- gp + theme(panel.grid = element_blank())
  }

  return(gp)
}
