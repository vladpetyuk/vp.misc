#' Volcano Plot
#'
#' A convenience function for creating a volcano plot.
#'
#' @param logFC a numeric vector of the log2 fold-change for each feature.
#' @param significance a numeric vector of significance values
#'        (p-value, q-value or adjusted p-value).
#' @param sig_threshold a numeric value specifying the significance threshold.
#'        A dashed horizontal line will be drawn at this value. The default
#'        is \code{NULL}.
#' @param threshold_line_color the color of the \code{sig_threshold} line.
#'        Default is \code{"red"}.
#' @param features a character vector of feature names or \code{NULL} (default).
#' @param top_n_features the number of top significant features to label on
#'        the plot if \code{features} is not \code{NULL}. Default is \code{5}.
#' @param feature_labels a character vector with the names of the features to
#'        be labeled. Any features that should not be labeled must be \code{NA}.
#' @param metadata a data frame containing information that should be
#'        accessible by additional layers. Useful for faceting.
#' @param point_color a color or a vector used to color the points.
#'        Default is \code{"grey"}.
#' @param point_size the size of the points. Default is \code{3}.
#' @param use_theme the plot theme. Default is \code{theme_bw()}.
#' @param label_outside_sig_range If \code{sig_threshold} is provided, this
#'        determines whether labels can be plotted between \code{sig_threshold}
#'        and 1. Default is \code{FALSE}, so labels will only be plotted in
#'        the range of significant features.
#' @param ... additional arguments passed to
#'          \code{\link[ggrepel]{geom_text_repel}}
#'
#' @return A ggplot object.
#'
#' @importFrom scales trans_new log_breaks pretty_breaks
#' @importFrom ggplot2 ggplot geom_point aes scale_y_continuous
#'             theme theme_bw xlim geom_hline
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr arrange
#'
#' @export plot_volcano
#'
#' @examples
#' library(MSnSet.utils)
#' data("cptac_oca")
#'
#' res <- limma_a_b(oca.set,
#'                  model.str = "~ PLATINUM.STATUS + AGE",
#'                  coef.str = "PLATINUM.STATUS")
#'
#' # Base plot
#' plot_volcano(logFC = res$logFC,
#'              significance = res$P.Value)
#'
#' # Label top 10 most significant features
#' plot_volcano(logFC = res$logFC,
#'              significance = res$P.Value,
#'              sig_threshold = 0.05,
#'              features = rownames(res),
#'              top_n_features = 10) +
#' ylab("p-value") +
#' xlab(expression(log[2](paste(fold, "-", change))))
#'
#' # Label top 5 most significant features and
#' # color by the sign of the average expression.
#' plot_volcano(logFC = res$logFC,
#'              significance = res$P.Value,
#'              sig_threshold = 0.05,
#'              features = rownames(res),
#'              point_color = factor(sign(res$AveExpr))) +
#' ylab("p-value") +
#' xlab(expression(log[2](paste(fold, "-", change))))


plot_volcano <- function(logFC,
                         significance,
                         sig_threshold = NULL,
                         threshold_line_color = "red",
                         features = NULL,
                         top_n_features = 5,
                         feature_labels = NULL,
                         metadata = NULL,
                         point_color = "grey",
                         point_size = 3,
                         use_theme = theme_bw(),
                         label_outside_sig_range = FALSE,
                         ...) {

  # Check input
  if (!is.null(features) & !is.null(feature_labels)) {
    stop("Can not provide both features and feature_labels.")
  }

  res <- data.frame(logFC, significance)

  # Add extra columns if they are provided
  if (!is.null(feature_labels)) {
    res <- cbind(res, feature_labels)
  }
  if (!is.null(features)) {
    res <- cbind(res, features)
  }
  if (!is.null(metadata)) {
    res <- cbind(res, metadata)
  }
  if (!is.null(point_color)) {
    res <- cbind(res, point_color)
  }

  # Once ggplot() is called, the data can not be changed
  # as new layers are added. This sorting is required by the
  # top n features block.
  res <- arrange(res, significance)

  y_min = min(significance, na.rm = T)
  breaks <- signif(10^(pretty_breaks()(0:floor(log10(y_min) * 2)) / 2), 1)
  x_extreme <- max(abs(range(logFC))) * 1.08

  # Base layer
  p <- ggplot(data = res)

  # Use horizontal lines instead of grid lines so that
  # the cutoff can be included in the y-axis
  if (!is.null(use_theme$panel.grid$colour)) {
    for (i in 1:length(breaks)) {
      p <- p + geom_hline(yintercept = breaks[i],
                          color = use_theme$panel.grid$colour,
                          size = use_theme$line$size)
    }
  }

  # Add sig_threshold to breaks, if provided
  if (!is.null(sig_threshold)) {
    breaks <- c(signif(sig_threshold, 1), breaks)
  }

  # Scatterplot layer
  if (length(point_color) == 1) {
    p <- p + geom_point(aes(x = logFC, y = significance),
                        size = point_size, color = point_color, alpha = 0.5)
  } else {
    p <- p + geom_point(aes(x = logFC, y = significance, color = point_color),
                        size = point_size, alpha = 0.5)
  }

  # Format axes
  log10_rev_trans <- trans_new("log10_rev",
                               function(x) -log10(x),
                               function(x) 10^(-x),
                               breaks = function(x) {
                                 y <- log_breaks(10)(x)
                                 rev(y)
                               },
                               domain = c(1e-100, Inf))

  p <- p +
    scale_y_continuous(trans = log10_rev_trans,
                       breaks = breaks,
                       limits = c(1, min(breaks))) +
    xlim(x_extreme*c(-1, +1)) +
    use_theme +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())

  # Add dashed line for significance cutoff
  if(!is.null(sig_threshold)) {
    p <- p +
      geom_hline(yintercept = sig_threshold,
                 color = threshold_line_color,
                 lty = "longdash")
  }

  # Label top n features or features provided by user
  if (!is.null(features)) {
    if(top_n_features < 1) {
      stop("If features is provided, top_n_features must be at least 1.")
    } else {
      # Get top n features
      res$features[which(!is.na(res$features))[
        top_n_features + 1:nrow(res)]] <- NA
      feature_labels <- res$features
    }
  }

  # Label top n features or user-provided features
  if (!is.null(feature_labels)) {

    if (!(label_outside_sig_range & is.null(sig_threshold))) {
      # Only plot labels within range of significant values
      ylim <- c(-log10(min(sig_threshold, 1)), NA)
    } else {
      ylim <- c(NA, NA)
    }

    p <- p +
      # Labels with positive logFC plotted to right of x = 0
      geom_label_repel(aes(x = logFC, y = significance,
                          label = ifelse(logFC >= 0, feature_labels, NA)),
                      min.segment.length = 0,
                      na.rm = TRUE, # silently remove missing values
                      xlim = c(0, NA),
                      ylim = ylim,
                      # Make label box partially transparent
                      fill = alpha("white", 0.5),
                      # Remove label box outline
                      label.size = NA,
                      ...) +
      # Labels with negative logFC plotted to left of x = 0
      geom_label_repel(aes(x = logFC, y = significance,
                           label = ifelse(logFC < 0, feature_labels, NA)),
                       min.segment.length = 0,
                       na.rm = TRUE, # silently remove missing values
                       xlim = c(NA, 0),
                       ylim = ylim,
                       fill = alpha("white", 0.5),
                       label.size = NA,
                       ...)
  }

  p
}

