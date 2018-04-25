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
#' Input is similar to \code{\link{plot_event}} or \code{\link{plot_expr}}.
#'
#' If \code{cluster_rows = TRUE}, then a heirarchical clustering using \code{hclust}
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
#' a path to the \code{.config} file, or 4-column data frame of the \code{.config}
#' file.
#' @param expr Logical - \code{TRUE} if plotting cRPKMs, \code{FALSE} otherwise
#' @param xlab The x-axis label
#' @param ylab The y-axis label
#' @param title Title of the plot
#' @param cluster_rows Logical to cluster rows using heirarchical clustering
#' @param cluster_cols Logical to cluster columns using heirarchical clustering
#' @param fill A vector of colours. e.g. from \code{colorRampPalette}.
#' Default is \code{NULL}, which will choose the palette automatically.
#' @param usepkg Default is \code{pheatmap}, which calls \code{\link[pheatmap]{pheatmap}}.
#' Other methods are also supported: \code{\link[gplots]{heatmap.2}} and
#' \code{\link[ggplot2]{geom_tile}}.
#' @param ... Additional parameters passed to \code{\link[pheatmap]{pheatmap}}
#' or \code{\link[gplots]{heatmap.2}}
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
#' @examples
#' # Uses ggplot2 by default
#' plot_multi(psi)
#' plot_multi(psi, config = config)
#'
#' # Use expr = TRUE for cRPKMs
#' plot_multi(crpkm, expr = TRUE)
#' plot_multi(crpkm, config = config, expr = TRUE)
#' plot_multi(crpkm, config = config, expr = TRUE, cluster_rows = TRUE)
#'
#' # To use pheatmap or gplots for plotting, set usepgk option
#' plot_multi(psi, config = config, usepkg = "pheatmap")
plot_multi <- function(df, config = NULL, expr = FALSE, xlab = "", ylab = "",
           title = "", cluster_rows = TRUE, cluster_cols = ifelse(is.null(config), TRUE, FALSE),
           fill = NULL, usepkg = "ggplot2", ...
           ) {
  match.arg(usepkg, c("gplots", "ggplot2", "pheatmap"))
  # Format input
  formatted_df <- format_table(df, expr = expr)
  if (expr == FALSE) {
    rownames(formatted_df) <- make_title.2(df$GENE, df$EVENT)
  }
  reordered <- preprocess_sample_colors(formatted_df, config = config, expr = expr)

  if (is.null(fill)) {
    if (expr) {
      # Use RColorBrewer::brewer.pal(6, "YlOrRd")
      fill <- c("#FFFFB2", "#FED976", "#FEB24C", "#FD8D3C", "#F03B20", "#BD0026")
    } else {
      fill <- colorRampPalette(c("yellow", "blue"))(20)
    }
  }

  if (usepkg == "gplots" && requireNamespace("gplots", quietly = TRUE)) {
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
    if (!is.null(config)) {
      col_colors <- reordered$col
      gplots::heatmap.2(as.matrix(reordered$data),
              Colv = cluster_cols, Rowv = cluster_rows,
              dendrogram = dendro,
              ColSideColors = col_colors,
              col = fill,
              margins = c(10, 35),
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
                margins = c(10, 35),
                trace = "none",
                key.xlab = ifelse(expr, "cRPKM", "PSI"),
                xlab = xlab, ylab = ylab,
                main = title,
                ...
      )
    }
  } else if (usepkg == "pheatmap" && requireNamespace("pheatmap", quietly = TRUE)) {
    if (!is.null(config)) {
      col_colors <- reordered$col
      anno_colors <- list(Group = reordered$group.col)
      anno_col <- data.frame(Group = reordered$config$GroupName)
      rownames(anno_col) <- reordered$config$SampleName
    } else {
      col_colors <- NA
      anno_colors <- NA
      anno_col <- NA
    }
    pheatmap::pheatmap(as.matrix(reordered$data),
                       cluster_row = cluster_rows,
                       clsuter_col = cluster_cols,
                       main = title,
                       col = fill,
                       annotation_colors = anno_colors,
                       annotation_col = anno_col,
                       ...
    )
  } else {
    if (cluster_rows) {
      # Perform heirarchical clustering of events/genes
      hr <- hclust(dist(reordered$data))
      reordered$data <- reordered$data[hr$order,]
    }

    if (is.null(config)) {
      hc <- hclust(dist(t(reordered$data)))
      reordered$data <- reordered$data[,hc$order]
    }

    reordered$data$id <- rownames(reordered$data)

    m <- suppressMessages(melt(reordered$data))
    m$variable <- factor(m$variable, levels = colnames(reordered$data),
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


    return(gp + scale_fill_gradientn(colours = fill, na.value = "white",
                                     name = ifelse(expr == FALSE, "PSI", "cRPKM")))
  }
}
