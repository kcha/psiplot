#' Draw group means as horizontal lines
draw_group_means <- function(reordered, config) {
  # Draw horizontal lines for groups
  seen <- vector()
  groups <- names(reordered$group.index)
  for (t in 1:length(groups)) {
    mu <- mean(as.numeric(
      reordered$data[reordered$group.index[[groups[t]]]]
    ), na.rm=TRUE)
    abline(h=mu, col=reordered$group.col[groups[t]], lwd=0.5)
    seen <- append(seen, paste(groups[t]))
  }

  # plot legend for mean group values
  if (length(seen) > 0) {
    legend_position <- ifelse(reordered$data[ncol(reordered$data)] > 50,
                              "bottomright", "topright")
    legend(legend_position, legend = seen, lty = 1, col = reordered$group.col,
           title = "Group Means", cex = 0.7, ncol = 2)
  }
}
