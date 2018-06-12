#' Plot cRPKM values of a given gene
#'
#' Generate a plot of cRPKM values for a given gene. The cRPKM values are
#' obtained from the
#' \href{https://www.github.com/vastgroup/vast-tools}{vast-tools} pipeline.
#'
#' @details
#' Like \code{\link{plot_event}} and \code{\link{plot_multievent}}, plots can be
#' customized via the \code{config} option. Either a data frame or
#' the filepath to the config file can be used. Alternatively, plots can be
#' customized using a limited set of graphical parameters as described above.
#'
#' See Details of \code{\link{plot_event}} and \code{\link{preprocess_sample_colors}}
#' for more information on customizing plots. Note that the \code{errorbar}
#' argument is not available for cRPKMs.
#'
#' cRPKM values that have \emph{NA} value are omitted and not plotted.
#'
#' @param x A 1-row data frame containing cRPKM values to be plotted
#' @param trim_colnames String that must be searched for and trimmed at the end
#' of every sample column in x. Useful to trim the "-cRPKM" suffix from expression
#' tables. If no string must be trimmed, leave as \code{NULL}.
#' @param config Optional configuration settings for \code{plot_expr}. Can be
#' a path to the \code{.config} file, or 4/5-column data frame of the \code{.config}
#' file. Use the latter option if you are calling \code{plot_expr} multiple times.
#' @param subg Logical indicating whether samples should be subgrouped for plotting.
#' If \code{TRUE}, the average of all samples in a subgroup is plotted as a single
#' data point. See \code{\link{plot_event}} and \code{\link{preprocess_sample_colors}}
#' for more details on subgrouping.
#' @param counts Logical indicating whether the data frame contains read counts.
#' Set to \code{TRUE} if the data frame contains two columns per sample (cRPKM and
#' counts), otherwise leave as \code{FALSE} (default).
#' @param groupmean Logical indicating whether grouped means should be drawn.
#' Requires \code{config}.
#' @param col Vector of colors with length matching the number of samples. If
#' specified, this will override the color settings specified in \code{config}.
#' @param title Title of the plot. If \code{NULL} (default), the title will be
#' the content of the \code{ID} column in \code{x}.
#' @param xlab The x-axis label.
#' @param ylab The y-axis label.
#' @param ylim Range of y-axis.
#' @param cex.main Plot title size (pts).
#' @param cex.yaxis Y-axis font size (pts).
#' @param cex.xaxis X-axis font size (i.e. the sample names) (pts).
#' @param pch Point symbol.
#' @param cex.pch Size of datapoints.
#' @param gridlines Logical indicating whether grid lines should be drawn.
#' @param plot (deprecated) prints the plot.
#' @param show_group_legend Set to FALSE to avoid showing a legend with the sample
#' groups and their colors.
#' @return ggplot2 object.
#' @seealso
#' \code{\link{format_table}} for performing some initial conversion steps of \code{x}
#'
#' \code{\link{preprocess_sample_colors}} for pre-processing of \code{x} using
#' \code{config}
#'
#' @export
#' @import methods
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @importFrom magrittr "%>%"
#' @examples
#' \dontrun{
#' plot_expr(crpkm[1,])
#'
#' # Plot with custom configuration
#' config
#' plot_expr(crpkm[1,], config = config, groupmean=TRUE)
#' plot_expr(crpkm[1,], config = "/path/to/config")
#'
#' # Plot using custom configuration, changing point symbol, and y-axis
#' # scale
#' plot_expr(crpkm[1,], config = config, pch = 9, ylim = c(20, 80))
#' }
#'
#' #' # Plot with subgrouped samples
#' plot_expr(crpkm[1,], config=config, subg=TRUE)
#'
#' # Plot directly from a table with suffixes and read counts
#' crpkm_counts
#' plot_expr(crpkm_counts[1,], config = config, trim_colnames = "-cRPKM", counts = TRUE)
#'
plot_expr <- function(
  x, config = NULL, subg = FALSE, trim_colnames = NULL, counts= FALSE,
  groupmean = ifelse(is.null(config), FALSE, TRUE), col = NULL,
  title = NULL, xlab = "", ylab = "Expression (cRPKM)", ylim = NULL,
  cex.main = 14, cex.yaxis = 12, cex.xaxis = 12,
  pch = 20, cex.pch = 3, plot = NULL, gridlines = TRUE, show_group_legend = TRUE) {

  if (!missing(plot)) {
    warning("The option 'plot' has been deprecated")
  }

  if (nrow(x) != 1) {
    stop("Too many rows!")
  }

  x <- format_table(x, expr = TRUE, counts=counts, trim_colnames=trim_colnames,
                    short_ids= FALSE)

  reordered <- preprocess_sample_colors(x,
                                        expr=T,
                                        config = config,
                                        col = col,
                                        subg = subg,
                                        multi_col = NULL)
  crpkm <- reordered$data

  subg <- all(subg==TRUE,
              "SubgroupName" %in% colnames(reordered$original_config))

  if (all(is.na(crpkm))) {
    warning("Did not find any points to plot")
  }

  N <- ifelse(is.null(ncol(crpkm)),1,ncol(crpkm))
  if (N < 2) {
    stop("Need two or more samples!")
  }

  # Set plot title
  if (is.null(title)) {
    title <- x$ID
  }

  if (is.null(ylim)) {
    ylim <- c(0, round(max(x[,-1])) + 1)
  }

  mdata <- suppressMessages(gather(crpkm,
                                 key = "SampleName",
                                 value = "value"))

  #Here the only difference between subg==T and subg==F is that samples need to
  #be pooled for cRPKM calculations. The rest is just the same (error bars are)
  #not an option for expression at the moment

  sm <- left_join(mdata,reordered$subgroup,by="SampleName")

  if(subg){

    smsum <- sm %>%
      dplyr::group_by(SubgroupName) %>%
      dplyr::summarise(value=mean(value,na.rm=T)) %>%
      mutate(value=replace(value,is.na(value),NA))

    } else {

      smsum <- sm

    }

  smsum <- smsum %>%
    dplyr::select(SubgroupName,value)

  smsum <- left_join(reordered$subgroup_order,smsum,by="SubgroupName") %>%
    dplyr::arrange(SubgroupOrder)

  smsum <- left_join(smsum,reordered$group,by="SubgroupName") %>%
    left_join(reordered$group_order,by="GroupName") %>%
    arrange(SubgroupOrder) %>%
    mutate(SubgroupName=factor(SubgroupName,levels=unique(SubgroupName))) %>%
    arrange(GroupOrder) %>%
    mutate(GroupName=factor(GroupName,levels=unique(GroupName))) %>%
    dplyr::select(Order=SubgroupOrder,
                  Sample=SubgroupName,
                  value,
                  GroupName,
                  RColorCode)

  gp <- ggplot(data=smsum) +
    geom_blank(data=smsum,
               aes(x=Sample,
                   y=value))

  gp <- gp + geom_point(aes(x=Sample,
                            y=value,
                            colour=GroupName),
                        size=cex.pch,
                        shape=pch,
                        show.legend=show_group_legend) +
    scale_colour_manual("Sample Group", values=reordered$group_order$RColorCode)

  if(!is.null(config) && groupmean) {
    gp <- draw_group_means(gp,
                           mdata,
                           reordered)
  }

  gp <- gp + ylab(ylab) +
    ylim(ylim) +
    ggtitle(title) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = cex.xaxis),
          axis.text.y = element_text(size = cex.yaxis),
          axis.title.y = element_text(size = cex.yaxis),
          title = element_text(size = cex.main))

  if (!gridlines) {
    gp <- gp + theme(panel.grid = element_blank())
  }

  return(gp)
}
