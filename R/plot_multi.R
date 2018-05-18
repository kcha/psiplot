#' Plot multiple events or genes as a heatmap (experimental)
#'
#' Visualize multiple PSI or cRPKM values in a single plot.
#'
#' @details
#' By default, \code{\link[ggplot2]{geom_tile}} is used to generate a heatmap
#' and a \code{ggplot2} object is returned. Alternatively, other heatmap
#' packages are supported including:  \code{\link[pheatmap]{pheatmap}}
#' (use  \code{usepkgs = "pheatmap"}) and
#' \code{\link[gplots]{heatmap.2}} (use \code{usepkgs = "gplots"}).
#'
#' Input is similar to \code{\link{plot_event}}, \code{\link{plot_expr}}, or
#' \code{\link{plot_multievent}}. Subgrouping of samples is also supported (see
#' \code{\link{plot_event}} and \code{\link{preprocess_sample_colors}} for details
#' on subgrouping). The PSI for each subgroup is taken as the average of all the
#' samples in the subgroup.
#'
#' If \code{cluster_rows = TRUE}, then a hierarchical clustering using \code{hclust}
#' will be performed on a distance matrix computed by \code{dist}. If \code{config}
#' is not specified, then the samples will also be clustered.
#'
#' To set the colours, use \code{fill} option to specify vector of colours.
#' By default (\code{fill=NULL}), PSI uses a yellow/blue gradient, while cRPKMs
#' uses the "YlOrRd" color brewer palette.
#'
#' @param df A data frame of input values (PSI or cRPKM). If the latter, need to
#' set \code{expr = TRUE}.
#' @param config Optional configuration settings. Can be
#' a path to the \code{.config} file, or 4/5-column data frame of the \code{.config}
#' file.
#' @param subg Logical indicating whether samples should be subgrouped for plotting.
#' @param expr Logical - \code{TRUE} if plotting cRPKMs, \code{FALSE} otherwise
#' @param counts Logical indicating whether the data frame contains read counts.
#' Set to \code{TRUE} if the data frame contains two rows per sample (cRPKM and
#' counts), otherwise leave as \code{FALSE} (default).
#' @param trim_colnames String that must be searched for and trimmed at the end
#' of every sample column in x. Useful to trim the "-cRPKM" suffix from expression
#' tables. If no string must be trimmed, leave as \code{FALSE}.
#' @param xlab The x-axis label
#' @param ylab The y-axis label
#' @param title Title of the plot
#' @param cluster_rows Logical to cluster rows using hierarchical clustering
#' @param cluster_cols Logical to cluster columns using hierarchical clustering
#' @param fill A vector of colours. e.g. from \code{colorRampPalette}.
#' Default is \code{NULL}, which will choose the palette automatically.
#' @param usepkg Default is \code{ggplot2}, which creates the heatmap using
#' \code{\link[ggplot2]{geom_tile}}. Otherwise, use \code{gplots}, which calls
#' \code{\link[gplots]{heatmap.2}}.
#' @param ... Additional parameters passed to \code{\link[gplots]{heatmap.2}} or
#' \code{\link[pheatmap]{pheatmap}}.
#' @export
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom grDevices colorRampPalette
#' @examples
#' plot_multi(psi)
#' plot_multi(psi, config = config)
#'
#' # Use expr = TRUE for cRPKMs
#' plot_multi(crpkm, expr = TRUE)
#' plot_multi(crpkm, config = config, expr = TRUE)
#' plot_multi(crpkm, config = config, expr = TRUE, cluster_rows = TRUE)
#'
#' # To use ggplot2 (or if gplots is not installed)
#' plot_multi(psi, config = config, usepkg = "ggplot2")
plot_multi <- function(df, config = NULL, subg = TRUE, expr = FALSE, counts = FALSE,
                       trim_colnames = FALSE, xlab = "", ylab = "", title = "",
                       cluster_rows = TRUE,
                       cluster_cols = ifelse(is.null(config), TRUE, FALSE),
                       fill = NULL, usepkg = c("ggplot2","gplots","pheatmap"),
                       ... ) {

  message(paste("plot_multi() is under active development.",
                "Please report bugs or feedback to https://github.com/kcha/psiplot/issues."))
  match.arg(usepkg)

  # Format input
  formatted_df <- format_table(df,
                               expr = expr,
                               counts = counts,
                               trim_colnames = trim_colnames,
                               short_ids = FALSE)

  if (expr == FALSE) {
    formatted_df$ID <- make_title.2(df$GENE, df$EVENT)
  }

  reordered <- preprocess_sample_colors(formatted_df,
                                        config = config,
                                        subg = subg,
                                        expr = expr)

  psi <- reordered$data %>%
    mutate(ID=formatted_df$ID) %>%
    dplyr::select(ID,colnames(reordered$data))

  mdata <- psi %>%
    group_by(ID) %>%
    gather(key="SampleName",
           -ID,
           value="value")


  #Determine if samples must be subgrouped
    subg <- all(c(subg==TRUE,
                "SubgroupName" %in% colnames(reordered$original_config)))

  #Do subgroups if subg==T, use the samples as subgroups if subg==FALSE.
  #This part could be a separate function, called also from plot_multievent.

  if(subg){
    sm <- left_join(mdata,reordered$subgroup,by="SampleName")
    smsum <- sm %>%
      dplyr::group_by(ID,SubgroupName) %>%
      dplyr::summarise(value=mean(value,na.rm=T)) %>%
      mutate(value=replace(value,is.na(value),NA)) %>%
      dplyr::select(ID,SubgroupName,value)

  } else{

    sm <- left_join(mdata,reordered$subgroup,by="SampleName")

    smsum <- sm %>%
      dplyr::select(ID,SubgroupName,value)

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

  heatmap_data  <- smsum %>%
    dplyr::select(ID,Sample,value) %>%
    spread(key=Sample,value=value) %>%
    as.data.frame()

  rownames(heatmap_data) <- heatmap_data$ID
  heatmap_data <- subset(heatmap_data, select=(-ID))

  if (is.null(fill)) {
    if (expr) {
      # Use RColorBrewer::brewer.pal(6, "YlOrRd")
      fill <- c("#FFFFB2", "#FED976", "#FEB24C", "#FD8D3C", "#F03B20", "#BD0026")
    } else {
      fill <- colorRampPalette(c("yellow", "blue"))(20)
    }
  }


  if (usepkg == "gplots" && requireNamespace("gplots", quietly = TRUE)) {
    #Make heatmap with gplots

    #Preformat data for gplots::heatmap.2
    col_colors <- smsum %>%
      ungroup() %>%
      dplyr::filter(ID==smsum$ID[1]) %>%
      pull(RColorCode)

    dendro <- "none"
    if (cluster_cols && cluster_rows) {
      dendro <- "both"
    } else if (cluster_cols) {
      dendro <- "column"
    } else if (cluster_rows) {
      dendro <- "row"
    }
    # Not ideal, but need two function calls if we want option to have
    # ColSideColors or not. Seems like heatmap.2 can't take an empty
    # ColSideColors -- you have to specify something or don't use it at all.
    if (!is.null(reordered$original_config)) {
      gplots::heatmap.2(as.matrix(heatmap_data),
              Colv = cluster_cols, Rowv = cluster_rows,
              dendrogram = dendro,
              ColSideColors = col_colors,
              col = fill,
              margins = c(10, 25),
              trace = "none",
              key.xlab = ifelse(expr, "cRPKM", "PSI"),
              xlab = xlab, ylab = ylab,
              main = title,
              ...
              )

    } else {
      gplots::heatmap.2(as.matrix(reordered$data),
                Colv = cluster_cols, Rowv = cluster_rows,
                dendrogram = dendro,
                col = fill,
                margins = c(10, 25),
                trace = "none",
                key.xlab = ifelse(expr, "cRPKM", "PSI"),
                xlab = xlab, ylab = ylab,
                main = title,
                ...
      )
    }
  } else if (usepkg == "pheatmap" && requireNamespace("pheatmap", quietly = TRUE)) {
    if(!is.null(reordered$original_config)){
      smsum1 <- smsum %>%
        ungroup() %>%
        dplyr::filter(ID==smsum$ID[1])

      col_colors <- smsum1 %>% pull(RColorCode)
      names(col_colors) <- smsum1 %>% pull(Sample)

      anno_colors <- reordered$group_order %>% pull(RColorCode)
      names(anno_colors) <- reordered$group_order %>% pull(GroupName)
      anno_colors <- list(Group=anno_colors)

      anno_col <- smsum1 %>% select(Group=GroupName) %>% as.data.frame()
      rownames(anno_col) <- smsum1 %>% pull(Sample)

    } else{
      col_colors <- NA
      anno_colors <- NA
      anno_col <- NA
    }

    pheatmap::pheatmap(as.matrix(heatmap_data),
                       cluster_rows = cluster_rows,
                       cluster_cols = cluster_cols,
                       main = title,
                       col = fill,
                       annotation_colors = anno_colors,
                       annotation_col = anno_col,
                       ...
    )
  } else {

    #We go with ggplot2

    if (cluster_rows) {
      # Perform hierarchical clustering of events/genes
      hr <- hclust(dist(heatmap_data))
      heatmap_data <- heatmap_data[hr$order,]
    }

    if (is.null(reordered$original_config)) {
      hc <- hclust(dist(t(heatmap_data)))
      heatmap_data <- heatmap_data[,hc$order]
    }

    heatmap_data$id <- rownames(heatmap_data)
    m <- gather(heatmap_data,key="variable",value="value",-id)

    m$variable <- factor(m$variable,
                         levels = colnames(subset(heatmap_data, select = -id)),
                         ordered = TRUE)

    m$id <- factor(m$id, levels = unique(m$id), ordered = TRUE)

    gp <- ggplot(m, aes(x = variable, y = id)) +
      geom_tile(aes(fill = value)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle= 45, hjust = 1, size = 9),
            axis.text.y = element_text(size = 8),
            panel.grid = element_blank(),
            panel.border = element_blank()) +
      xlab(xlab) + ylab(ylab) + ggtitle(title) +
      coord_fixed(ratio = 1)

    gp <- gp + scale_fill_gradientn(colours = fill, na.value = "white",
                                    name = ifelse(expr == FALSE, "PSI", "cRPKM"))


    return(gp)
  }
}
