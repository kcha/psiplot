#' Plot PSI values of multiple alternative splicing events
#'
#' Generate a plot with PSI values for several exons. The PSI values and
#' corresponding quality scores are typically obtained from the
#' \href{https://www.github.com/vastgroup/vast-tools}{vast-tools} pipeline.
#'
#' @details
#' Like in \code{\link{plot_event}} and \code{\link{plot_expr}}, plots can be
#' customized via the \code{config} option. Either a data frame or
#' the filepath to the config file can be used. Alternatively, plots can be
#' customized using a limited set of graphical parameters as described above.
#'
#' See Details of \code{\link{plot_event}} and
#' \code{\link{preprocess_sample_colors}} for more information on the usage of
#' the \code{config}, \code{subg}, and \code{errorbar} arguments.
#'
#' Unlike in \code{\link{plot_event}} and \code{\link{plot_expr}},
#' sample colors in \code{config} or \code{col} are now shown only in the
#' background, as the point colors are now used to differentiate PSI values from
#' different events. In addition, a line connects the points from each event, to
#' increase visibility. The color of each event can be set using the
#' \code{event_col} argument.
#'
#' Also, note that using \code{subg=TRUE} and \code{errorbar=TRUE} together is
#' an experimental feature, computationally expensive, and CIs may not be shown
#' for some subgroups, especially for events with low coverage and PSI values near
#' 0 or 1 (see \code{\link{get_beta_ci_subg}} for details on error bar estimation
#' for subgrouped samples).
#'
#' @param x A data frame containing PSI values to be plotted.
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
#' @param errorbar Logical indicating whether error bars should be drawn
#' @param col Vector of colors with length matching the number of samples. If
#' specified, this will override the color settings specified in \code{config}.
#' @param event_col Vector of colors, with length matching the number of events
#' (rows of \code{x}). If left as \code{NULL}, event colors will be set with
#' \code{\link[ggplot2]{scale_colour_hue}}
#' @param title Title of the plot.
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
#' @param lwd Line width for errorbars and the line connecting PSIs from each event.
#' @param show_event_legend Set to FALSE to avoid showing a legend with the event
#' IDs and their colors.
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
#' \code{\link{plot_event}} for plotting single events.
#'
#' @export
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import ggplot2
#' @examples
#' plot_multievent(psi, config = config)
#'
#' # Example with subgrouped samples, custom title and no error bars
#' plot_multievent(psi, config = config, subg = TRUE, errorbar = FALSE, title = "Highlighted events")
#'
#' # Legends can be hidden separately
#' plot_multievent(psi, config = config, show_event_legend = FALSE)
#' plot_multievent(psi, config = config, show_group_legend = FALSE)
#'
#' # Custom event colors
#' plot_multievent(psi[1:3,], config = config, event_col=c("red","black","orange"))
#'
#' # Use of errobar = TRUE and subg = TRUE is experimental and computationally expensive
#' \dontrun{
#' plot_multievent(psi, config = config, subg = TRUE, errorbar = TRUE)
#' }
#'
plot_multievent <- function(
  x, config = NULL,  subg = FALSE, trim_colnames = NULL,
  qual = c("VLOW","N","LOW","OK","SOK"), errorbar = TRUE,
  col = NULL,  event_col = NULL, title = "MULTI EVENT PLOT", xlab = "",
  ylab = "PSI", ylim = c(0,100), cex.main = 14, cex.yaxis = 12, cex.xaxis = 12,
  pch = 20, cex.pch = 3, plot = NULL, gridlines = TRUE, lwd=0.5,
  show_event_legend=T, show_group_legend=T) {

  if (!missing(plot)) {
    warning("The option 'plot' has been deprecated")
  }

  qual = match.arg(qual)

  # Format input
  x <- format_table(x,
                    qual = qual,
                    trim_colnames = trim_colnames,
                    short_ids = T)

  reordered <- preprocess_sample_colors(x,
                                        config,
                                        subg = subg,
                                        col = col,
                                        multi_col = event_col)

  psi <- reordered$data %>%
    mutate(ID=x$ID) %>%
    dplyr::select(ID,colnames(reordered$data))

  qual <- reordered$qual %>%
    mutate(ID=x$ID) %>%
    dplyr::select(ID,colnames(reordered$qual))

  #Check two or more samples
  if(any(rowSums(!is.na(psi[,-1]))/(ncol(psi)-1) == 0)){
    warning("Could not find any points to plot for some events")
  }

  #Subgroups will only be used if the original config had subgroups AND
  #if the subg argument was set to TRUE. Otherwise subgroups will be overridden
  #with the sample names (1 sample = 1 subgroup)

  subg <- all(c(subg==TRUE,
                "SubgroupName" %in% colnames(reordered$original_config)))

  # Set plot title
  if (is.null(title)) {
    title <- "MULTI_EVENT_PLOT"
  }

  mdata <- psi %>%
    group_by(ID) %>%
    gather(key="SampleName",
           -ID,
           value="value")

  mqual <- qual %>%
    group_by(ID) %>%
    gather(key="SampleName",
           -ID,
           value="qual")


  #Pool samples according to subgroups in reordered
  sm <- left_join(mdata,mqual,by=c("ID","SampleName")) %>%
    left_join(reordered$subgroup,by="SampleName")

  smsum <- sm %>%
    dplyr::group_by(ID,SubgroupName)

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
      select(ID,
             SubgroupName,
             value = m_psi,
             lo,hi)
  } else {
    smsum <- smsum %>%
      do(m_psi = mean(.$value,na.rm=T)) %>%
      ungroup() %>%
      mutate(m_psi = map_dbl(m_psi,1)) %>%
      mutate(m_psi = replace(m_psi,is.na(m_psi),NA)) %>%
      select(ID,
             SubgroupName,
             value = m_psi)
  }


  smsum <- left_join(reordered$subgroup_order,smsum,by="SubgroupName") %>%
    dplyr::arrange(ID,SubgroupOrder)

  smsum <- left_join(smsum,reordered$group,by="SubgroupName") %>%
    left_join(reordered$group_order,by="GroupName") %>%
    arrange(ID,SubgroupOrder) %>%
    group_by(ID) %>%
    mutate(SubgroupName=factor(SubgroupName,levels=unique(SubgroupName))) %>%
    arrange(ID,GroupOrder) %>%
    mutate(GroupName=factor(GroupName,levels=unique(GroupName))) %>%
    dplyr::select(ID,
                  Order=SubgroupOrder,
                  Sample=SubgroupName,
                  colnames(smsum),
                  GroupName,
                  RColorCode)

    #Now on to the plotting

  gp <- ggplot() +
    geom_blank(data=smsum,
               aes(x=Sample,
                   y=value))

  if(!is.null(reordered$original_config)){

    #Background with groups

    groupbg <- smsum %>% dplyr::filter(ID==smsum$ID[1])

    gp <- gp + geom_rect(data=groupbg,
                         aes(xmin=as.numeric(Sample)-0.5,
                             xmax=as.numeric(Sample)+0.5,
                             ymin=-Inf,
                             ymax=Inf,
                             fill=GroupName),
                         alpha=0.1)

    gp <- gp + scale_fill_manual("Sample Group",
                                 values=reordered$group_order$RColorCode)

  }

  #Core plotting of events
  gp <- gp + geom_point(data=smsum,
                        aes(x=Sample,
                            y=value,
                            colour=ID),
                        size = cex.pch,
                        shape=pch)

  if(!is.null(event_col)){
    gp <- gp + scale_colour_manual("ID",
                                   values=reordered$multi_col$EventRColorCode)
    } else{
      gp <- gp + scale_colour_hue("ID")

      }

  #Lines joining points from the same event
  gp <- gp + geom_line(data=smsum,
                       linetype="solid",
                       lwd=lwd,
                       aes(x=Sample,
                           y=value,
                           colour=ID,
                           group=ID),
                       show.legend = F)

  #Draw error bars
  if(errorbar==TRUE){
    gp <- gp + geom_errorbar(data=smsum,
                             inherit.aes = F,
                             mapping = aes(x=Sample,
                                           ymin=lo,
                                           ymax=hi,
                                           colour=ID),
                             width=0.1,
                             show.legend = F)
  }

  #Final touches
  gp <- gp +
    ylab("PSI") +
    ylim(ylim)

  gp <- gp + ggtitle(title)

  gp <- gp +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = cex.xaxis),
          axis.text.y = element_text(size = cex.yaxis),
          axis.title.y = element_text(size = cex.yaxis),
          title = element_text(size = cex.main))

  if (!gridlines) {
    gp <- gp + theme(panel.grid = element_blank())
  }

  #Hide legends according to arguments
  if(show_group_legend==FALSE){
    gp <- gp + guides(fill=show_group_legend)
  }
  if(show_event_legend==FALSE){
    gp <- gp + guides(colour=show_event_legend)
  }

  return(gp)

}

