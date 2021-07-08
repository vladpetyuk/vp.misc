#' UMAP Plot
#'
#' A convenience function for creating a UMAP scatter plot
#' for samples in an ExpressionSet/MSnSet object.
#'
#' @param eset eset (or most likely eset subclass) object
#' @param phenotype character; one of the \code{colnames(pData(eset))}
#' @param standardize logical; should the rows of the expression data be
#'          standardized? (default: \code{FALSE})
#' @param pca logical; whether PCA should be performed prior to UMAP
#'          (default: \code{TRUE}).
#' @param n_perm integer; if pca is \code{TRUE}, this is the number of times
#'          that the expression matrix should be permuted (default: 10).
#'          This is done to determine the optimal number of principal
#'          components to use.
#' @param n_neighbors integer; number of nearest neighbors to use for UMAP
#'          (default: \eqn{\sqrt(n samples)}).
#' @param n_epochs integer; number of iterations performed during layout
#'          optimization.
#' @param min_dist numeric; determines how close points appear in the final
#'          layout (default: 0.01).
#' @param legend_title character; the title for the legend. If NULL, is set to
#'          the wrapped up phenotype.
#' @param legend_title_width integer; wrapping up too long legend titles.
#'          Passed to \code{stringr::str_wrap()} as width argument
#'          (default: 20).
#' @param point_size numeric; size of the plotted points (default: 4).
#' @param show_ellipse logical; whether to plot 50\% confidence ellipses.
#' @param ellipse_type character; the type of ellipse. The default "norm"
#'          assumes a multivariate normal distribution, and "t" assumes a
#'          multivariate t-distribution.
#' @param show_na logical; should the data points for which phenotype is
#'          unknown be shown? (default: \code{TRUE})
#' @param ... Other arguments that can be passed to the UMAP configuration.
#'          See \code{\link[umap]{umap.defaults}} for more details.
#'
#' @return A ggplot object
#'
#' @importFrom MSnbase exprs pData
#' @import     umap
#' @importFrom ggplot2 ggplot geom_point coord_fixed theme_bw
#'                      guides guide_legend stat_ellipse aes
#' @importFrom stringr str_wrap
#'
#' @examples
#' data(srm_msnset)
#' plot_umap(msnset,
#'           phenotype = "subject.type",
#'           show_ellipse = FALSE)
#' plot_umap(msnset,
#'           phenotype = "subject.type",
#'           min_dist = 0.5,
#'           show_ellipse = TRUE)
#' plot_umap(msnset,
#'           phenotype = "subject.type",
#'           pca = FALSE,
#'           show_ellipse = TRUE)
#' plot_umap(msnset)
#'
#' @export plot_umap

plot_umap <- function(eset, phenotype = NULL, standardize = FALSE,
                      pca = TRUE, n_perm = 10,
                      n_neighbors = round(sqrt(length(sampleNames(eset)))),
                      n_epochs = 1000, min_dist = 0.01,
                      legend_title = NULL, legend_title_width = 20,
                      point_size = 4,
                      show_ellipse = TRUE, ellipse_type = "norm",
                      show_na = TRUE,
                      ...) {

  if (!is.null(phenotype)) {
    colorBy <- pData(eset)[[phenotype]]
    if (!show_na) {
      idx <- !is.na(colorBy)
      eset <- eset[, idx]
      colorBy <- colorBy[idx]
    }
  } else {
    colorBy <- ""
    phenotype <- ""
    show_ellipse <- FALSE
  }
  stopifnot(sum(complete.cases(exprs(eset))) > 1)
  eset <- eset[complete.cases(exprs(eset)), ]

  z <- t(exprs(eset))
  if (standardize) {
    z <- sweep(z, 1, rowMeans(z), FUN = "-")
    z <- sweep(z, 1, apply(z, 1, sd), FUN = "/")
  }

  # PCA ---
  if (pca) {
    PC <- prcomp(z, scale. = F)
    expl_var <- PC$sdev^2 / sum(PC$sdev^2)

    # Determine the optimal number of principal components
    expl_var_perm <- matrix(NA, ncol = length(PC$sdev), nrow = n_perm)
    for (k in 1:n_perm) {
      expr_perm <- apply(exprs(eset), 2, sample)
      z_perm <- t(expr_perm)
      if (standardize) {
        z_perm <- sweep(z_perm, 1, rowMeans(z_perm), FUN = "-")
        z_perm <- sweep(z_perm, 1, apply(z_perm, 1, sd), FUN = "/")
      }
      PC_perm <- prcomp(z_perm, scale. = F)
      expl_var_perm[k, ] <- PC_perm$sdev^2 / sum(PC_perm$sdev^2)
    }

    pval <- apply(t(expl_var_perm) >= expl_var, 1, sum) / n_perm
    optPC <- which(pval >= 0.05)[1] - 1

    z <- as.matrix(PC$x[, 1:optPC])
  }

  # UMAP ---
  res_umap <- umap(z,
                   n_neighbors = n_neighbors,
                   n_epochs = n_epochs,
                   min_dist = min_dist,
                   ...)

  umap_df <- data.frame(UMAP1 = res_umap$layout[, 1],
                        UMAP2 = res_umap$layout[, 2],
                        colorBy = colorBy)

  if (is.null(legend_title)) {
    legend_title <- str_wrap(phenotype, legend_title_width)
  }

  if(show_ellipse) {
    p <- ggplot(umap_df) +
      # Plot ellipses behind points
      stat_ellipse(aes(x = UMAP1, y = UMAP2, fill = colorBy),
                   geom = "polygon", type = ellipse_type, level = 0.5,
                   alpha = 0.15, show.legend = TRUE) +
      geom_point(aes(x = UMAP1, y = UMAP2, color = colorBy),
                 size = point_size, shape = 20, show.legend = TRUE) +
      guides(color = guide_legend(legend_title),
             fill = guide_legend(legend_title)) +
      theme_bw()
  } else {
    p <- ggplot(umap_df) +
      geom_point(aes(x = UMAP1, y = UMAP2, color = colorBy),
                 size = point_size, shape = 20, show.legend = TRUE) +
      guides(color = guide_legend(legend_title)) +
      theme_bw()
  }

  return(p)
}

