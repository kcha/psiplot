#' Draw group means as horizontal lines
#'
#' When this function is called, the group means are calculated and plotted
#' as horizontal lines. For internal use only.
#'
#' @param gp ggplot2 object
#' @param mdata a gathered data frame of PSI/cRPKM values. Usually comes from
#' \code{\link{plot_event}} or \code{\link{plot_expr}}.
#' @param reordered a list with the subgroup and group structure of the samples.
#' Usually comes from \code{\link{preprocess_sample_colors}}.
#' @return ggplot2 object
#' @import dplyr
#' @seealso \code{\link{preprocess_sample_colors}}, \code{\link{plot_event}},
#' \code{\link{plot_expr}}
#'
draw_group_means <- function(gp,mdata,reordered) {

  m <- mdata %>%
    left_join(reordered$subgroup,by="SampleName") %>%
    left_join(reordered$group,by="SubgroupName") %>%
    group_by(GroupName) %>%
    dplyr::summarise(mu=mean(value,na.rm = T)) %>%
    mutate(mu=replace(mu,is.na(mu),NA)) %>%
    left_join(reordered$group_order,by="GroupName") %>%
    dplyr::select(GroupOrder,GroupName,mu,RColorCode) %>%
    arrange(GroupOrder) %>%
    mutate(GroupName=factor(GroupName,levels=GroupName))


  gp2 <- gp + geom_hline(data=m,
                         aes(yintercept = mu,
                             colour= GroupName),
                         lwd=0.5,
                         show.legend = F)

  return(gp2)
}
