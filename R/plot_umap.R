#' Create UMAP scatterplots of samples
#'
#' A convenience function for creating UMAP scatter plots of samples in an
#' ExpressionSet/MSnSet object.
#'
#' @param msnset an MSnSet object.
#' @param phenotype character; one of the \code{colnames(pData(msnset))}.
#' @param n_neighbors integer; number of nearest neighbors to use for UMAP.
#'          If \code{NULL} (default) will use the closest integer to the
#'          square root of the number of samples. If a vector of multiple
#'          values are supplied, multiple plots will be generated and arranged
#'          with \code{\link[ggplot2]{facet_wrap}}.
#' @param min_dist numeric; determines how close points appear in the final
#'          layout (default: 0.1).
#' @param n_epochs integer; number of iterations performed during layout
#'          optimization (default: 1000).
#' @param point_size numeric; size of the plotted points (default: 3).
#' @param show_ellipse logical; whether to plot ellipses.
#' @param ellipse_type character; the type of ellipse. "norm" (the default)
#'          assumes a multivariate normal distribution, and "t" assumes a
#'          multivariate t-distribution.
#' @param ellipse_level numeric; if \code{ellipse_type="euclid"}, the radius
#'          of the circle to be drawn; otherwise, the confidence level
#'          (default: 0.5).
#' @param show_na logical; whether points with missing phenotype labels
#'          should be plotted (default: \code{TRUE}).
#' @param ... Other arguments that can be passed to
#'          \code{\link[umap]{umap.defaults}} or
#'          \code{\link[ggplot2]{facet_wrap}}.
#'
#' @return A ggplot object
#'
#' @importFrom umap umap umap.defaults
#' @importFrom MSnbase exprs pData sampleNames
#' @importFrom ggplot2 ggplot geom_point theme_bw stat_ellipse aes
#'             facet_wrap sym vars
#' @importFrom data.table rbindlist
#' @importFrom stats complete.cases
#'
#' @export plot_umap
#'
#' @examples
#' # Load msnset
#' data(srm_msnset)
#' # Do not color by phenotype
#' plot_umap(msnset)
#' # Plots for several values of n_neighbors
#' plot_umap(msnset, n_neighbors = c(6, 8, 10))
#' # Points colored by subject.type
#' plot_umap(msnset, phenotype = "subject.type",
#'           show_ellipse = TRUE)
#' # Plots for several values of n_neighbors with
#' # points colored by subject.type
#' plot_umap(msnset, phenotype = "subject.type",
#'           n_neighbors = c(6, 8, 10),
#'           min_dist = 0.01, show_ellipse = TRUE)


plot_umap <- function(msnset,
                      phenotype = NULL,
                      n_neighbors = NULL,
                      min_dist = 0.1,
                      n_epochs = 1000,
                      point_size = 3,
                      show_ellipse = TRUE,
                      ellipse_type = "norm",
                      ellipse_level = 0.5,
                      show_na = TRUE,
                      ...) {

  ## Check input
  # Subset to complete cases
  stopifnot(sum(complete.cases(exprs(msnset))) > 1)
  if (min_dist == 0) {
    stop("min_dist must be greater than 0.")
  }

  if (!is.null(phenotype)) {
    colorBy <- pData(msnset)[[phenotype]]
    if (!show_na) {
      idx <- !is.na(colorBy)
      msnset <- msnset[, idx]
      colorBy <- colorBy[idx]
    }
  } else {
    colorBy <- ""
    phenotype <- "none"
    show_ellipse <- FALSE
  }

  # If the number of nearest neighbors is not specified,
  # use the integer closest to the square root of the
  # number of samples.
  if (is.null(n_neighbors)) {
    n_neighbors <- round(sqrt(length(sampleNames(msnset))))
  }

  msnset <- msnset[complete.cases(exprs(msnset)), ]

  # Take the transpose so that features are columns
  # and samples are rows.
  z <- t(exprs(msnset))

  # Additional arguments
  user_args <- list(...)
  umap_args <- user_args[names(user_args) %in% names(umap.defaults)]
  facet_args <- user_args[names(user_args) %in% names(formals(facet_wrap))]

  # If n_neighbors is a vector of multiple values, do this
  umap_df <- lapply(n_neighbors, function(k) {
    # UMAP
    u <- do.call(umap, c(list(d = z, n_neighbors = k, n_epochs = n_epochs),
                         umap_args))

    # Data frame with columns UMAP1, UMAP2, n_neighbors, and colorBy
    umap_df <- as.data.frame(u$layout)
    colnames(umap_df) <- c("UMAP1", "UMAP2")
    umap_df$n_neighbors <- k
    umap_df$colorBy <- colorBy
    umap_df
  })
  umap_df <- rbindlist(umap_df)

  # Rename colorBy to the string stored as phenotype
  colnames(umap_df)[4] <- phenotype

  # Base plot
  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) + theme_bw()

  if (show_ellipse) {
    # Plot ellipses beneath points
    p <- p +
      stat_ellipse(aes(fill = !!sym(phenotype)),
                   geom = "polygon", type = ellipse_type,
                   level = ellipse_level, alpha = 0.15)
  }

  # Add points
  if (phenotype != "none") {
    # color by phenotype
    p <- p + geom_point(aes(color = !!sym(phenotype)), size = point_size)
  } else {
    # Do NOT color by phenotype
    p <- p + geom_point(size = point_size)
  }

  if (length(n_neighbors) > 1) {
    p <- p +
      do.call(facet_wrap, c(list(facets = vars(n_neighbors)), facet_args))
  }

  return(p)
}

utils::globalVariables(c("UMAP1", "UMAP2"))

