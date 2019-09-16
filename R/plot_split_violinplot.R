GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

#' Plot split violin plot (convenient to compare distribution before and after)
#' 
#' @inheritParams proBatch
#'
#' @param df data.frame with \code{y_col}, \code{col_for_color}, \code{col_for_box}
#' @param y_col value to explore the distribution of
#' @param col_for_color column to use to map to two colors
#' @param col_for_box column to use to do group comparison
#' @param colors_for_plot colors to map to col_for_color
#' @param hlineintercept NULL: no intercept line; non-null: intercept value
#'
#' @return ggplot object
#' @export
#'
#' @examples
plot_split_violin_with_boxplot <- function(df, y_col = 'y', 
                                           col_for_color = 'm', col_for_box = 'x', 
                                           colors_for_plot = c('#8f1811', '#F8C333'),
                              hlineintercept = NULL, plot_title = NULL,
                              theme = 'classic'){
  
  p <- ggplot(df, aes(x = !!sym(col_for_box), y = !!sym(y_col), fill = !!sym(col_for_color))) + 
    geom_split_violin(alpha = .75) + 
    geom_boxplot(aes(group = interaction(!!sym(col_for_color), !!sym(col_for_box))), 
                 position = 'dodge', width =.2, outlier.shape = NA, coef = 0)
  if (!is.null(colors_for_plot)){
    p = p + scale_fill_manual(values = colors_for_plot)
  }
  
  if (!is.null(hlineintercept)){
    p = p + geom_hline(yintercept = 0, linetype = 'dashed')
  }
  
  if(theme == 'classic'){
    p = p + theme_classic()
  }
  
  if(!is.null(plot_title)){
    p = p +ggtitle(plot_title)
  }
  return(p)  
}

