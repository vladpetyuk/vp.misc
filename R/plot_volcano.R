#' @title Volcano Plot
#'
#' @description A convenience function for creating a volcano plot.
#'
#' @param df \code{data.frame} or object that can be coerced to a
#'   \code{data.frame}
#' @param logFC character; the name of the column in \code{df} containing log2
#'   fold-changes.
#' @param pvals character; the name of the column in \code{df} containing
#'   p-values (adjusted or not).
#' @param sig_threshold \code{NULL} (default) or numeric between 0 and 1; the
#'   threshold indicating statistical significance. Creates a dashed horizontal
#'   line at that value. Any points above this line in the graph are
#'   statistically significant.
#' @param label \code{NULL} or character; the name of the column in \code{df}
#'   used to label the top most significant features. See details for more.
#' @param num_features numeric; the number of most significant features to
#'   label. Default is 8.
#' @param point_args a list of arguments passed to
#'   \code{\link[ggplot2]{geom_point}}.
#' @param label_args a list of arguments passed to
#'   \code{\link[ggrepel]{geom_label_repel}}.
#'
#'
#' @details \code{sig_threshold} will create a secondary axis if the threshold
#'   is not one of the existing y-axis breaks.
#'
#'   To label specific features (not top \code{num_features}), \code{label}
#'   needs to be the name of a column where all but the features that will be
#'   labeled are \code{NA}. Also, set \code{num_features} to \code{nrow(df)}.
#'
#'
#' @return A ggplot object.
#'
#' @examples
#' library(ggplot2) # additional plot modifications
#' library(MSnSet.utils)
#' data("cptac_oca")
#'
#' # Differential analysis
#' df <- limma_a_b(oca.set,
#'                 model.str = "~ PLATINUM.STATUS + AGE",
#'                 coef.str = "PLATINUM.STATUS")
#' df$features <- rownames(df) # Add feature names column
#'
#' # Base plot
#' p <- plot_volcano(df, logFC = "logFC", pvals = "P.Value")
#' p
#'
#' # Change y-axis title
#' p + labs(y = "Unadjusted p-value")
#'
#' # Add dashed line at y = 0.05
#' plot_volcano(df, logFC = "logFC", pvals = "P.Value",
#'              sig_threshold = 0.05)
#'
#' # Label top 8 most significant features
#' plot_volcano(df, logFC = "logFC", pvals = "P.Value",
#'              label = "features")
#'
#' # Change point opacity, point size, and color of labels
#' plot_volcano(df, logFC = "logFC", pvals = "P.Value",
#'              label = "features",
#'              point_args = list(alpha = 0.3, size = 3),
#'              label_args = list(color = "royalblue"))
#'
#' @importFrom dplyr %>% rename sym arrange select mutate slice_min left_join
#' @importFrom ggplot2 ggplot aes labs theme_bw expansion sec_axis
#'   scale_y_continuous geom_hline
#' @importFrom ggrepel geom_label_repel
#' @importFrom scales pretty_breaks alpha trans_new log_breaks
#' @importFrom utils modifyList
#'
#' @export plot_volcano
#'


plot_volcano <- function(df, logFC, pvals, sig_threshold = NULL,
                         label = NULL, num_features = 8L,
                         point_args = list(), label_args = list()) {

  # Check that logFC and pvals are valid
  if (any(sapply(c(logFC, pvals), is.null))) {
    stop("logFC and pvals must be strings")
  }

  # Check that all columns are present in df
  if (!all(c(logFC, pvals, label) %in% colnames(df))) {
    stop("One or more column names provided is not present in df")
  }

  # Check that num_features is valid
  if (!(num_features %in% 1:nrow(df))) {
    stop("num_features must be an integer between 1 and nrow(df)")
  }

  # Check that sig_threshold is valid --- this could be improved
  if (!is.null(sig_threshold)) {
    if (!is.numeric(sig_threshold) |
        !(0 < sig_threshold & sig_threshold < 1)) {
      stop("sig_threshold must be NULL or a value between 0 and 1")
    }
  }

  # Check that point_args and label_args are valid
  if (!is.list(point_args) | !is.list(label_args)) {
    stop("point_args and label_args must be lists. See ?list for help")
  }

  # Rename columns and sort by significance
  df <- as.data.frame(df) %>%
    dplyr::rename(logFC = !!sym(logFC),
                  pvals = !!sym(pvals)) %>%
    arrange(pvals)

  ## y-axis breaks
  # log10 transform the minimum value, multiply by 2.16,
  # and floor it to extend the lower bound. Divide the
  # sequence from 0 to this value by 2 to partially undo
  # the 2.16 scaling. Create pretty breaks, undo log10
  # transform, and round to 1 significant digit. Use curly
  # braces to prevent piping to first argument.
  breaks <- min(df$pvals, na.rm = TRUE) %>%
    # ifelse(. > 1e-3, 1e-3, .) %>% # set base range
    {signif(10^(pretty_breaks()(0:floor(log10(.) * 2.16)) / 2), 1)}

  # y-axis labels
  labels <- ifelse(breaks > 1e-3,
                   format(breaks, scientific = FALSE, drop0trailing = TRUE),
                   format(breaks))

  # Arguments for geom_point
  point_args <- list(na.rm = TRUE) %>%
    # Allow user-supplied args to overwrite defaults
    modifyList(val = point_args, keep.null = TRUE)

  # Base plot
  p <- ggplot(data = df, mapping = aes(x = logFC, y = pvals)) +
    do.call(what = geom_point, args = point_args) +
    labs(x = expression(paste("log"[2], "(Fold-Change)")),
         y = pvals) +
    theme_bw()

  # Arguments for y-axis/axes
  scale_args <- list(trans = log10_rev_trans,
                     breaks = breaks,
                     labels = labels,
                     expand = expansion(mult = c(0, 0.05)),
                     limits = c(1, min(breaks)))

  # If sig_threshold is outside plotting window, set to NULL
  if (!is.null(sig_threshold)) {
    if (sig_threshold < min(breaks)) {
      message("sig_threshold is outside the y-axis range. Setting to NULL")
      sig_threshold <- NULL
    }
  }

  # If significance threshold is not already in breaks, create a
  # secondary axis. If we were to forcibly include it in the
  # existing breaks, they would be unevenly spaced and the
  # horizontal grid lines would look bad.
  if (!is.null(sig_threshold)) {
    if (!(sig_threshold %in% breaks)) {
      # Create a secondary axis
      scale_args <- c(
        scale_args,
        # Identity transform and break at sig_threshold
        list(sec.axis = sec_axis(trans = ~ .,
                                 breaks = sig_threshold)
        )
      )
    }
  }

  # Modify y axis
  p <- p + do.call(what = scale_y_continuous, args = scale_args) +
    # Dashed line indicating significance cutoff.
    # If NULL, no line is plotted.
    geom_hline(yintercept = sig_threshold,
               lty = "longdash")

  # Label the top n_features
  if (!is.null(label)) {

    # Select top most significant rows to label
    label_df <- df %>%
      # select(pvals, !!sym(label)) %>%
      mutate(feature_labels = !!sym(label)) %>%
      slice_min(order_by = pvals, n = num_features)

    suppressMessages(
      p$data <- left_join(df, label_df) # update plot data
    )

    # Arguments for geom_label_repel
    label_args <- list(
      mapping = aes(label = feature_labels),
      fill = alpha(colour = "white", alpha = 0.5), # translucent background
      min.segment.length = 0,
      max.overlaps = Inf,
      na.rm = TRUE
    ) %>%
      # Allow user-supplied args to overwrite defaults
      modifyList(val = label_args, keep.null = TRUE)

    # Add labels
    p <- p + do.call(what = geom_label_repel, args = label_args)
  }

  return(p)
}


utils::globalVariables(c(".", "feature_labels"))


# More intuitive axis for p-values
log10_rev_trans <- trans_new("log10_rev",
                             transform = function(x) -log10(x),
                             inverse = function(x) 10 ^ (-x),
                             breaks = function(x) rev(log_breaks(10)(x)))

