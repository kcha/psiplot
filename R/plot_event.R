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
#' The order of samples (e.g. columns in \code{x}) as it appears on the
#' resulting plot can be customized using a config file. If a config file is not
#' used and re-ordering is desired, then it must be done manually before calling
#' \code{\link{plot_event}} by altering the columns of \code{x}.
#'
#' If subgroups are defined in \code{config} and \code{subg=TRUE}:
#' \itemize{
#'   \item{The mean PSI of the samples within each subgroup are drawn as a single
#'   data point.}
#'   \item{If \code{subg=TRUE} and \code{errorbar=TRUE}, confidence intervals for
#'   the subgroup are derived from a joint beta distribution, fitted to a set of
#'   points derived from individual beta distributions corresponding to each
#'   sample (see \code{\link{get_beta_ci_subg}}). This is an experimental feature,
#'   computationally expensive, and CIs may not be shown for some subgroups,due
#'   to failure in fitting the joint distribution, especially for events with low
#'   coverage and PSI values near 0 or 1.}
#'   \item{Subgroups will be ordered by the minimum \emph{Order} value of their
#'   samples, and assigned to the first group to which they are matched.}
#'   }
#'
#' If no subgroups are defined in \code{config}, or \code{subg=FALSE}, a subgroup
#' will be defined for each sample, preserving their name and order.
#'
#' If groups and colors are defined in \code{config}, all the samples in a group
#' will be colored with the first color assigned to the group. Else, all
#' samples will be plotted in black.
#'
#' If groups are defined in \code{config} and \code{groupmean=TRUE}, the mean
#' PSI of the samples within each group are drawn as horizontal lines. The
#' colors of the lines are determined by the color set to the first sample of
#' each group by \code{RColorCode} in \code{config}. A corresponding legend key
#' will also be drawn.
#'
#' (See also \code{\link{preprocess_sample_colors}} to see details on how configs
#' are used).
#'
#' PSI values that have \emph{NA} value are omitted and not plotted.
#'
#' Error bars based on the confidence interval of the PSI estimation can be
#' shown by setting \code{errorbar=TRUE}.
#'
#' @param x A 1-row data frame containing PSI values to be plotted.
#' @param trim_colnames String that must be searched for and trimmed at the end
#' of every sample column in x. If no string must be trimmed, leave as \code{FALSE}.
#' @param config Optional configuration settings for \code{plot_event}. Can be
#' a path to the \code{.config} file, or 4/5-column data frame of the \code{.config}
#' file. Use the latter option if you are calling \code{plot_event} multiple times.
#' @param subg Logical indicating whether samples should be subgrouped for plotting.
#' @param qual String indicating the minimun \emph{vast-tools} quality score
#' for the PSI to be accepted. Defaults to \code{'VLOW'}. See the
#' \href{https://github.com/vastgroup/vast-tools/blob/master/README.md}{vast-tools
#' documentation} for details.
#' @param errorbar Logical indicating whether error bars should be drawn.
#' @param groupmean Logical indicating whether grouped means should be drawn.
#' Requires \code{config}.
#' @param col Vector of colors with length matching the number of samples. If
#' specified, this will override the color settings specified in \code{config}.
#' @param title Title of the plot. If \code{NULL} (default), the title generated
#' by \code{\link{make_title}}.
#' @param xlab The x-axis label.
#' @param ylab The y-axis label.
#' @param ylim Range of y-axis.
#' @param cex.main Plot title size (pts).
#' @param cex.yaxis Y-axis font size (pts).
#' @param cex.xaxis X-axis font size (i.e. the sample names) (pts).
#' @param pch Point symbol.
#' @param cex.pch Size of datapoints.
#' @param plot (deprecated) Prints the plot.
#' @param gridlines Logical indicating whether grid lines should be drawn.
#' @param show_group_legend Set to FALSE to avoid showing a legend with the sample
#' groups and their colors.
#' @return ggplot2 object
#' @seealso
#' \code{\link{format_table}} for performing some initial conversion steps of
#' \code{x}.
#'
#' \code{\link{preprocess_sample_colors}} for pre-processing of \code{x}
#' using\code{config}.
#'
#' \code{\link{plot_multievent}} for plotting more than one
#' event in the same plot.
#'
#' @export
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @importFrom magrittr "%>%"
#' @import ggplot2
#' @examples
#' \dontrun{
#' plot_event(psi[1,])
#
#' # Plot with custom configuration
#' config
#' plot_event(psi[1,], config = config, groupmean=TRUE)
#' plot_event(psi[1,], config = "/path/to/config")
#'
#' # Plot with subgrouped samples
#' config
#' plot_event(psi[1,], config = config, subg = TRUE, errorbar = FALSE)
#'
#' # Plot with subroups and error bars (slow, experimental feature)
#' \dontrun{
#' plot_event(psi[1,], config = config, subg = TRUE, errorbar = TRUE)
#' }
#'
#' # Plot using custom configuration, changing point symbol, and y-axis
#' # scale
#' plot_event(psi[1,], config = config, pch = 9, ylim = c(20, 80))
#' }
plot_event <- function(
  x, config = NULL, subg = TRUE, trim_colnames = NULL,
  qual = c("VLOW","N","LOW","OK","SOK"), errorbar = TRUE,
  groupmean = ifelse(is.null(config), FALSE, TRUE), col = NULL,
  title = NULL, xlab = "", ylab = "PSI", ylim = c(0,100),
  cex.main = 14, cex.yaxis = 12, cex.xaxis = 12,
  pch = 20, cex.pch = 3, plot = NULL, gridlines = TRUE,
  show_group_legend = TRUE) {

  if (!missing(plot)) {
    warning("The option 'plot' has been deprecated")
  }

  if (nrow(x) != 1) {
    stop("Too many rows!")
  }

  # Format input
  x <- format_table(x,
                    qual = qual,
                    trim_colnames=trim_colnames,
                    short_ids = FALSE)
  reordered <- preprocess_sample_colors(x,
                                        expr=F,
                                        config,
                                        col = col,
                                        subg= subg,
                                        multi_col=NULL)

  psi <- reordered$data
  qual <- reordered$qual

  #Subgroups will only be used if the original config had subgroups AND
  #if the subg argument was set to TRUE. Otherwise subgroups will be overridden
  #with the sample names (1 sample = 1 subgroup)

  subg <- all(c(subg==TRUE,
                "SubgroupName" %in% colnames(reordered$original_config)))

  if (all(is.na(psi))) {
    warning("Did not find any points to plot")
  }

  N <- ifelse(is.null(ncol(psi)), 1, ncol(psi))
  if (N < 2) {
    stop("Need two or more samples!")
  }

  # Set plot title
  if (is.null(title)) {
    title <- make_title(as.character(x$ID))
  }

  mdata <- suppressMessages(gather(psi,
                                 key = "SampleName",
                                 value = "value"))
  mqual <- suppressMessages(gather(qual,
                                 key = "SampleName",
                                 value = "qual"))


  # Pool samples according to subgroups using the info in reordered
  # (Subgroups equal to samples if subg == F anyway)
  sm <- left_join(mdata,mqual,by="SampleName") %>%
    left_join(reordered$subgroup,by="SampleName")


  smsum <- sm %>%
    dplyr::group_by(SubgroupName)

  if(errorbar){
    smsum <- smsum %>%
      do(m_psi = mean(.$value,na.rm=T),
         ci = get_beta_ci_subg(.$value,.$qual)) %>%
      ungroup() %>%
      mutate(m_psi = map_dbl(m_psi,1),
             ci = map(ci,1),
             lo = map_dbl(ci,1),
             hi = map_dbl(ci,2)) %>%
      mutate(m_psi = replace(m_psi,is.na(m_psi),NA)) %>%
      select(SubgroupName,
             value = m_psi,
             lo,hi)
  } else {
    smsum <- smsum %>%
      do(m_psi = mean(.$value,na.rm=T)) %>%
      ungroup() %>%
      mutate(m_psi = map_dbl(m_psi,1)) %>%
      mutate(m_psi = replace(m_psi,is.na(m_psi),NA)) %>%
      select(SubgroupName,
             value = m_psi)
    }

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
                  colnames(smsum),
                  GroupName,
                  RColorCode)

  gp <- ggplot(data=smsum) +
    geom_blank(data=smsum,
               aes(x=Sample,
                   y=value))

  gp <- gp + geom_point(aes(x=Sample,
                            y=value,
                            colour=GroupName),
                        size = cex.pch,
                        shape=pch,
                        show.legend = show_group_legend) +
    scale_colour_manual("Sample Group", values=reordered$group_order$RColorCode)

  if(errorbar){
    gp <- gp + geom_errorbar(inherit.aes = F,
                             aes(x=Sample,
                                 ymin=lo,
                                 ymax=hi,
                                 colour=GroupName),
                             width=0.05,
                             show.legend = F)
  }

  # Draw horizontal lines for groups
  if (!is.null(config) && groupmean) {
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

  browser()
  return(gp)
}

#' Make plot title
#'
#' Create a plot title using the event ID
#'
#' @param x A character containing the event ID like
#' "\code{S|TSPAN6|chrX:99885756-99885863|108}"
#' @return A character with a human-friendly title with
#' event type, gene symbol, event coordinates, and length
#' @examples
#' \dontrun{
#' f <- format_table(psi)
#' print(make_title(f$ID[1]))
#' }
make_title <- function(x) {
  event <- strsplit(as.character(x), split = "\\|")[[1]]
  sprintf("%s\n(%s, %s bp, type %s)",
                   event[2], event[3], event[4], event[1])
}

#' Make plot title (version 2)
#'
#' Create a plot title using the event ID. Returns symbol and event ID.
#'
#' @param gene A vector or character string
#' @param event A vector or character string with same size as \code{gene}
#' @return A character with a human-friendly title with format:
#' \code{GENE (EVENT ID)}
#' @examples
#' \dontrun{
#' print(make_title.2(psi$GENE, psi$EVENT))
#' }
make_title.2 <- function(gene, event) {
  sprintf("%s (%s)", gene, event)
}


