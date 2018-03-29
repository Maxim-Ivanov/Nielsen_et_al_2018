drawMetagenePlot <- function(mat_list, drawCI=TRUE, x.axis=FALSE, title = "Metagene_plot", xlabel = "Intervals of interest", ylabel = "Average signal", vline = FALSE, width=NA, height=NA, units="cm", plotPDF=TRUE) {
  require(ggplot2)
  if (identical(x.axis, FALSE)) {
    x.axis <- seq_len(ncol(mat_list[[1]]))
  }
  long_df <- data.frame()
  for (i in seq(length(mat_list))) {
    curr_mat <- mat_list[[i]]
    curr_name <- names(mat_list)[[i]]
    avg <- apply(curr_mat, 2, mean, na.rm = TRUE)
    df <- data.frame("pos" = x.axis, "avg" = avg, "group" = curr_name)
    if (isTRUE(drawCI)) {
      sem <- apply(curr_mat, 2, sd, na.rm = TRUE) / sqrt(nrow(curr_mat))
      lower <- avg - 1.96 * sem; upper <- avg + 1.96 * sem
      df <- cbind(df, data.frame("lower" = lower, "upper" = upper))
    }
    long_df <- rbind(long_df, df)
  }
  p <- ggplot(long_df, aes(x = pos, y = avg, color = group)) + geom_line(size=1) + ggtitle(title) + xlab(xlabel) + ylab(ylabel) + theme_bw()
  if (isTRUE(drawCI)) {
    p <- p + geom_ribbon(aes(x = pos, ymin = lower, ymax = upper, fill = group), alpha = 0.1)
  }
  if (is.numeric(vline)) {
    p <- p + geom_vline(xintercept = vline)
  }
  filename <- gsub(" ", "_", title)
  ggsave(filename = paste0(filename, ".png"), plot = p, width=width, height=height, units=units)
  if (isTRUE(plotPDF)) {
    ggsave(filename = paste0(filename, ".pdf"), plot = p, width=width, height=height, units=units)
  }
  return(long_df)
}

