#' PCA Plot
#'
#' A convenience function for plotting PCA scatter plot for samples in
#' ExpressionSet/MSnSet object. The biplot is essentially the ggplot version of
#' \code{\link[stats]{biplot.prcomp}}.
#'
#' @param eset eset (or most likely eset subclass) object.
#' @param phenotype \code{NULL} or string; one of \code{colnames(pData(eset))}.
#'   This is used to color the points or labels, if \code{label} is not
#'   \code{NULL}. Default is \code{NULL}, which colors all points or labels
#'   black.
#' @param label \code{NULL} or string; one of \code{colnames(pData(eset))}. If a
#'   string is provided, labels will be used instead of points.
#' @param z_score logical; whether to convert values to Z-Scores by sample.
#'   Default is \code{TRUE}.
#' @param standardize logical; if \code{TRUE} (default), the feature loadings
#'   and scores are scaled in opposite directions by the standard deviations of
#'   the principal components. This will produce a biplot similar to
#'   \code{\link[stats]{biplot.prcomp}}. If \code{FALSE}, the result will be
#'   similar to \code{\link[stats]{biplot.default}}.
#' @param show_ellipse logical; whether to show the confidence ellipses if
#'   \code{phenotype} is not \code{NULL}.
#' @param components numeric; a vector of length two specifying the principal
#'   components to plot. Default is \code{c(1, 2)}, which plots PC1 on the
#'   x-axis and PC2 on the y-axis. Order matters.
#' @param biplot logical; whether to display the biplot.
#' @param biplot_labels \code{NULL} or string; the name of a column in
#'   \code{fData(eset)} used to label the biplot features. If \code{NULL}
#'   (default), \code{featureNames(eset)} is used.
#' @param num_features numeric; the number of most influential features from
#'   each principal component to label. Default is \code{6}.
#' @param show_NA logical; whether to include samples with missing phenotype
#'   information. Default is \code{TRUE}.
#' @param legend_title string; title of the plot legend. Defaults to
#'   \code{phenodata}.
#' @param arrow_args a list of arguments passed to
#'   \code{\link[ggplot2]{geom_segment}} to modify the biplot arrows.
#' @param label_args a list of arguments passed to
#'   \code{\link[ggrepel]{geom_label_repel}} to modify the biplot labels.
#' @param ... additional arguments passed to \code{\link[ggplot2]{geom_point}}
#'   or \code{\link[ggplot2]{geom_text}}, such as \code{size} and \code{pch}.
#'
#' @return A ggplot object
#'
#'
#' @import ggplot2
#' @importFrom Biobase exprs pData fData
#' @importFrom ggrepel geom_label_repel
#' @importFrom stats complete.cases prcomp
#' @importFrom utils modifyList
#'
#' @export plot_pca
#'
#' @examples
#' library(MSnSet.utils)
#' data(cptac_oca)
#'
#' # Default plot
#' plot_pca(oca.set)
#'
#' # Color by SUBTYPE with custom legend title
#' plot_pca(oca.set, "SUBTYPE", legend_title = "Subtype")
#'
#' # Color by SUBTYPE and label by Batch
#' plot_pca(oca.set, phenotype = "SUBTYPE", label = "Batch")
#'
#' # Biplot
#' plot_pca(oca.set, biplot = TRUE)


plot_pca <- function(eset, phenotype = NULL, label = NULL, z_score = TRUE,
                     show_ellipse = TRUE, components = 1:2, biplot = FALSE,
                     biplot_labels = NULL, standardize = TRUE,
                     num_features = 6L, show_NA = TRUE,
                     legend_title = phenotype,
                     arrow_args = list(), label_args = list(), ...) {

  # Handling coloring by phenotype. Do this first, in case
  # rows are removed when show_NA = FALSE
  if (!is.null(phenotype)) {
    if (!phenotype %in% colnames(pData(eset))) {
      stop(
        sprintf("'%s' is not the name of a column in pData(eset).", phenotype)
      )
    }

    colorBy <- pData(eset)[, phenotype]
    # If not showing missing values, remove those samples
    if (!show_NA) {
      idx <- !is.na(colorBy)
      eset <- eset[, idx]
      colorBy <- colorBy[idx]
    }
  } else {
    show_ellipse <- FALSE
    colorBy <- NULL
  }

  # Check that components are valid
  if (length(components) != 2) {
    stop(sprintf("components must be a vector of length 2, not %d.",
                 length(components)))
  }
  if (!all(components %in% 1:ncol(eset))) {
    stop(sprintf("The values of components must be between 1 and %d.",
                 ncol(eset)))
  }

  complete_rows <- complete.cases(exprs(eset))

  # Check that there are enough complete rows for PCA
  if (sum(complete_rows) < 2) {
    stop("There are fewer than 2 rows with non-missing data.")
  }

  message(sprintf("Subsetting to %d complete rows for PCA.",
                  sum(complete_rows)))

  # Subset to complete rows
  eset <- eset[complete_rows, ]

  # If z_score, convert to Z-Scores by sample (row when transposed)
  if (z_score) {
    z <- t(scale(exprs(eset), center = TRUE, scale = TRUE))
  } else {
    z <- t(exprs(eset))
  }

  ## PCA
  # By default, center = TRUE, scale. = FALSE
  pca_res <- prcomp(z)

  u <- pca_res$x # Scores
  v <- pca_res$rotation # Eigenvectors

  if (standardize) {
    n <- nrow(u)
    lam <- pca_res$sdev * sqrt(n)

    # Scale u down and v up. Product is still the same
    u <- t(t(u) / lam)
    v <- t(t(v) * lam)
  }

  # Determine ratio between scale of v and u
  u_range <- apply(u[, components], 2, function(x) abs(range(x)))
  v_range <- apply(v[, components], 2, function(x) abs(range(x)))

  ratio <- max(v_range / u_range) # ratio for scaling v and secondary axes
  v <- v / ratio # scale v

  # Data frames for plotting
  df.u <- as.data.frame(u[, components])
  df.v <- as.data.frame(v[, components])

  # Percent of variance explained by each PC
  d <- pca_res$sdev # Standard deviations
  var_expl <- round(100 * d ^ 2 / sum(d ^ 2), digits = 2)[components]
  axis_labs <- sprintf("PC%d (%g%%)", components, var_expl)

  # If colorBy is not NULL, add that column to df
  if (!is.null(colorBy)) {
    df.u$colorBy <- colorBy
  }

  ## Visualization
  # Base plot
  p <- ggplot(data = df.u,
              mapping = aes(x = df.u[, 1], y = df.u[, 2], color = colorBy)) +
    geom_hline(yintercept = 0, lty = "longdash", color = "darkgrey") +
    geom_vline(xintercept = 0, lty = "longdash", color = "darkgrey") +
    labs(x = axis_labs[1], y = axis_labs[2]) +
    theme_bw() +
    theme(aspect.ratio = 1)

  # 50% confidence ellipse layer first so they are
  # beneath the layer of points or labels.
  if (show_ellipse & !is.numeric(colorBy)) {
    p <- p +
      stat_ellipse(mapping = aes(fill = colorBy, color = NULL),
                   geom = "polygon", type = "norm",
                   level = 0.5, alpha = 0.1, show.legend = TRUE)
  }

  # If label is NULL, add points. Otherwise, add labels
  if (is.null(label)) {
    p <- p +
      geom_point(...)
  } else {
    if (!label %in% colnames(pData(eset))) {
      stop(
        sprintf("'%s' is not the name of a column in pData(eset).", label)
      )
    }

    labels <- pData(eset)[, label]
    p <- p +
      geom_text(mapping = aes(label = labels), ...)
  }

  # Set titles for color and fill legend
  p <- p +
    guides(color = guide_legend(title = legend_title),
           fill = guide_legend(title = legend_title))

  # If colorBy is numeric, use a colorbar
  if (is.numeric(colorBy)) {
    p <- p +
      guides(color = guide_colorbar(title = legend_title))
  }

  ## Biplot
  if (biplot) {
    # Get the indices of the top influential features
    # from each principal component. num_features determines how
    # many to select from each component.
    top_features <- lapply(1:2, function(i) {
      order(abs(df.v)[, i], decreasing = TRUE)[1:num_features]
    })
    top_features <- unique(unlist(top_features))

    # Subset loadings to top features and rename columns
    df.v <- df.v[top_features, ]
    colnames(df.v) <- c("xend", "yend")
    df.v$x <- df.v$y <- 0

    # If biplot_labels is not provided, default to row names
    if (is.null(biplot_labels)) {
      df.v$labels <- rownames(df.v)
    } else {
      df.v$labels <- fData(eset)[top_features, biplot_labels]
    }

    scale_args <- list(expand = expansion(mult = rep(0.1, 2)),
                       sec.axis = sec_axis(~ . * ratio))

    # Arguments for geom_segment
    arrow_args <- list(mapping = aes(x = x, y = y, xend = xend, yend = yend),
                       arrow = arrow(length = unit(0.5, "line")),
                       data = df.v, color = "red3") %>%
      # Allow user-supplied args to overwrite defaults
      modifyList(val = arrow_args, keep.null = TRUE)

    # Arguments for geom_label_repel
    label_args <- list(mapping = aes(x = xend, y = yend, label = labels),
                       data = df.v,
                       color = arrow_args[["color"]],
                       max.overlaps = Inf,
                       min.segment.length = 0,
                       fill = alpha("white", 0.5)) %>%
      # Allow user-supplied args to overwrite defaults
      modifyList(val = label_args, keep.null = TRUE)

    # Add segments with arrows and text labels
    p <- p +
      # Add extra padding around plot area and secondary axes for v units
      do.call(scale_x_continuous, scale_args) +
      do.call(scale_y_continuous, scale_args) +
      do.call(geom_segment, arrow_args) +
      do.call(geom_label_repel, label_args) +
      theme(axis.text.y.right = element_text(color = arrow_args[["color"]]),
            axis.text.x.top = element_text(color = arrow_args[["color"]]),
            axis.ticks.y.right = element_line(color = arrow_args[["color"]]),
            axis.ticks.x.top = element_line(color = arrow_args[["color"]]))
  }

  return(p)
}

utils::globalVariables(c("xend", "yend"))







## Deprecated functions (to be removed)
# plot_pca_v1 <- function(eset, phenotype=NULL, show.ellipse=TRUE,
#                         show.NA=TRUE, legend.title.width=20){
#
#   # handling coloring by phenotype
#   if (!is.null(phenotype)) {
#     colorBy <- pData(eset)[[phenotype]]
#     if(!show.NA){
#       idx <- !is.na(colorBy)
#       eset <- eset[,idx]
#       colorBy <- colorBy[idx]
#     }
#   }
#   else {
#     colorBy <- ""
#    phenotype <- ""
#     show.ellipse <- FALSE
#   }
#
#   # get rid of NA values
#   stopifnot(sum(complete.cases(exprs(eset))) > 1)
#   eset <- eset[complete.cases(exprs(eset)),]
#
#   # PCA itself
#   px <- dudi.pca(exprs(eset), scannf = F, scale=T, center=T, nf = 2)
#   px <- px$co # this is plotting loadings rather then scores!!!
#   colnames(px) <- c("PC1", "PC2")
#
#   # visualize
#   ggdata <- data.frame(px, colorBy)
#   p <-
#     ggplot(ggdata) +
#     geom_point(aes(x=PC1, y=PC2, color=colorBy),
#                size=5, shape=20, show.legend = TRUE) +
#     coord_fixed() +
#     theme_bw()
#
#  # Ugly engtanglement of if/else statements. Needs to be improved.
#   phenotype_str <- str_wrap(phenotype, legend.title.width)
#   if(show.ellipse){
#     p <- p +
#       stat_ellipse(aes(x=PC1, y=PC2, fill=colorBy),
#                    geom="polygon", type="norm",
#                    level=0.5, alpha=0.1, show.legend = TRUE) +
#       guides(color=guide_legend(phenotype_str),
#              fill=guide_legend(phenotype_str))
#  }else{
#     if(is.numeric(colorBy)){
#       p <- p + guides(color=guide_colorbar(phenotype_str))
#     } else if (!identical(colorBy,"")){
#       p <- p + guides(color=guide_legend(phenotype_str))
#     }
#   }
#   return(p)
# }
#
# utils::globalVariables(c("PC1", "PC2"))
#
#
# # @describeIn plot_pca_v1 Alternative PCA
# # @importFrom pcaMethods pca
# # @importFrom RColorBrewer "brewer.pal"
# # @importFrom Biobase sampleNames
# # @importFrom graphics legend text grid
# #
# # @export plot_pca_v2
# #
# # @examples
# # data(srm_msnset)
# # plot_pca_v2(msnset)
# # plot_pca_v2(msnset, phenotype = "subject.type")
# # plot_pca_v2(msnset, phenotype = "subject.type", label = "match.group")
#
#
# plot_pca_v2 <- function(eset, phenotype = NULL, label = NULL){
#
#   # library("pcaMethods")
#   # library("RColorBrewer")
#
#   pcaResults <- pca(exprs(eset))
#
#   plot(pcaResults@loadings*1.4, pch='',
#        xlab = sprintf("PC1 (%g%%)", round(pcaResults@R2[1] * 100, digits = 1)),
#        ylab = sprintf("PC2 (%g%%)", round(pcaResults@R2[2] * 100, digits = 1)))
#
#   colz <- 'black'
#   if(!is.null(phenotype)){
#     coloring <- as.factor(pData(eset)[[phenotype]])
#     colScheme <- brewer.pal(max(3, nlevels(coloring)), "Set1")
#     text.col <- colScheme[seq_len(nlevels(coloring))]
#     legend("topleft",
#            legend=levels(coloring),
#            text.col=text.col,
#            bg="#EEEEEE80")
#     colz <- colScheme[coloring]
#     # title(sprintf("PCA colored by %s", "WHIM sample type"))
#   }
#
#   # If sample labels are not provided, use sample names
#   if (is.null(label)) {
#     labels <- sampleNames(eset)
#   } else {
#     labels <- pData(eset)[, label]
#   }
#
#   text(pcaResults@loadings,
#        labels = labels,
#        font=2, col=colz, adj=c(0.25,0))
#   grid()
#}
#
# utils::globalVariables(c("PC1", "PC2"))
#
#
# # @describeIn plot_pca_v1 Alternative PCA
# # @importFrom ggrepel geom_label_repel
# # @export plot_pca_v3
# #
# # @examples
# # data(srm_msnset)
# # plot_pca_v3(msnset, phenotype = "subject.type")
# # plot_pca_v3(msnset, phenotype = "subject.type", label = "sample.id")
# # plot_pca_v3(msnset)
#
# plot_pca_v3 <- function(eset, phenotype=NULL, label=NULL, point_size = 5,
#                         show.ellipse=TRUE,
#                         show.NA=TRUE, legend.title.width=20, ...){
#
#   # handling coloring by phenotype
#   if (!is.null(phenotype)) {
#     colorBy <- pData(eset)[[phenotype]]
#     if(!show.NA){
#       idx <- !is.na(colorBy)
#       eset <- eset[,idx]
#       colorBy <- colorBy[idx]
#     }
#   }
#   else {
#     colorBy <- ""
#     phenotype <- ""
#     show.ellipse <- FALSE
#   }
#
#   # get rid of NA values
#   stopifnot(sum(complete.cases(exprs(eset))) > 1)
#  eset <- eset[complete.cases(exprs(eset)),]
#
#   # PCA itself
#   z <- t(exprs(eset))
#   #  zero-center and scale because we more care about correlation
#   #  rather Eucledian dist
#   z <- sweep(z, 1, rowMeans(z), FUN = "-")
#   z <- sweep(z, 1, apply(z,1,sd), FUN = "/")
#   pca1 <- prcomp(z, scale. = F)
#   # create data frame with scores
#   scores <- as.data.frame(pca1$x)
#
#   exp_var <- 100 * summary(pca1)$importance[2,][c(1,2)] # from ggord
#   axes <- paste0("PC", c(1,2))
#   axes <- paste0(axes, " (", round(exp_var, 2), "%)")
#
#
#   # visualize
#   ggdata <- data.frame(scores[,c(1,2)], colorBy) # first two PCs
#   p <-
#     ggplot(ggdata) +
#     geom_point(aes(x=PC1, y=PC2, color=colorBy),
#                size=point_size, shape=20, show.legend = TRUE) +
#     coord_fixed() +
#     xlab(axes[1]) + ylab(axes[2]) +
#     theme_bw()
#
#
#   # Add labels by 'label'
#   if (!is.null(label)) {
#     custom_args <- list(mapping = aes(x = PC1, y = PC2, fill = colorBy,
#                                       label = pData(eset)[[label]]),
#                         fontface = 'bold', color = 'white',
#                         box.padding = 0.25, point.padding = 0.25,
#                         segment.color='grey50',
#                         label.size = 0.01,
#                         segment.alpha = 0.50,
#                         size = 2.5)
#     user_args <- list(...)
#     custom_args[names(user_args)] <- user_args
#
#     p <- p +
#       do.call(geom_label_repel, custom_args)
#   }
#
#   # Ugly engtanglement of if/else statements. Needs to be improved.
#   phenotype_str <- str_wrap(phenotype, legend.title.width)
#   if(show.ellipse){
#     p <- p +
#       stat_ellipse(aes(x=PC1, y=PC2, fill=colorBy),
#                    geom="polygon", type="norm",
#                    level=0.5, alpha=0.1, show.legend = TRUE) +
#       guides(color=guide_legend(phenotype_str),
#              fill=guide_legend(phenotype_str))
#   }else{
#     if(is.numeric(colorBy)){
#       p <- p + guides(color=guide_colorbar(phenotype_str))
#     } else if (!identical(colorBy,"")){
#       p <- p + guides(color=guide_legend(phenotype_str))
#     }
#   }
#   return(p)
# }
#
#
#
#
# # @describeIn plot_pca_v1 Alternative PCA with option to select pair of principal components
# # @importFrom ggrepel geom_label_repel
# # @export plot_pca_v4
# #
# # @examples
# # data(srm_msnset)
# # plot_pca_v4(msnset, phenotype = "subject.type", pc.index=c(1,3))
# # plot_pca_v4(msnset, phenotype = "subject.type", label = "sample.id", pc.index=c(1,3))
# # plot_pca_v4(msnset, pc.index=c(1,3))
#
# plot_pca_v4 <- function (eset, phenotype = NULL, label = NULL, point_size = 5,
#                          show.ellipse = TRUE, show.NA = TRUE, legend.title.width = 20,
#                          pc.index = c(1, 2),
#                          ...)
# {
#   if (!is.null(phenotype)) {
#     colorBy <- pData(eset)[[phenotype]]
#     if (!show.NA) {
#       idx <- !is.na(colorBy)
#       eset <- eset[, idx]
#       colorBy <- colorBy[idx]
#     }
#   }
#   else {
#     colorBy <- ""
#     phenotype <- ""
#     show.ellipse <- FALSE
#   }
#   stopifnot(sum(complete.cases(exprs(eset))) > 1)
#   eset <- eset[complete.cases(exprs(eset)), ]
#   z <- t(exprs(eset))
#   z <- sweep(z, 1, rowMeans(z), FUN = "-")
#  z <- sweep(z, 1, apply(z, 1, sd), FUN = "/")
#   pca1 <- prcomp(z, scale. = F)
#   scores <- as.data.frame(pca1$x)
#   exp_var <- 100 * summary(pca1)$importance[2, ][pc.index]
#   axes <- paste0("PC", pc.index)
#   axes <- paste0(axes, " (", round(exp_var, 2), "%)")
#   ggdata <- data.frame(scores[, pc.index], colorBy)
#   p <- ggplot(ggdata) + geom_point(aes(x = ggdata[,1], y = ggdata[,2], color = colorBy),
#                                    size = point_size, shape = 20, show.legend = TRUE) +
#    coord_fixed() + xlab(axes[1]) + ylab(axes[2]) + theme_bw()
#   if (!is.null(label)) {
#     custom_args <- list(mapping = aes(x = ggdata[,1], y = ggdata[,2],
#                                       fill = colorBy, label = pData(eset)[[label]]), fontface = "bold",
#                         color = "white", box.padding = 0.25, point.padding = 0.25,
#                         segment.color = "grey50", label.size = 0.01, segment.alpha = 0.5,
#                         size = 2.5)
#     user_args <- list(...)
#     custom_args[names(user_args)] <- user_args
#     p <- p + do.call(geom_label_repel, custom_args)
#   }
#   phenotype_str <- str_wrap(phenotype, legend.title.width)
#   if (show.ellipse) {
#     p <- p + stat_ellipse(aes(x = ggdata[,1], y = ggdata[,2], fill = colorBy),
#                           geom = "polygon", type = "norm", level = 0.5, alpha = 0.1,
#                           show.legend = TRUE) + guides(color = guide_legend(phenotype_str),
#                                                        fill = guide_legend(phenotype_str))
#   }
#   else {
#     if (is.numeric(colorBy)) {
#       p <- p + guides(color = guide_colorbar(phenotype_str))
#     }
#     else if (!identical(colorBy, "")) {
#       p <- p + guides(color = guide_legend(phenotype_str))
#     }
#   }
#   return(p)
# }

