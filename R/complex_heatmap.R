#' @title Annotated Expression or Correlation Heatmaps
#'
#' @description A wrapper around functions from
#'   \code{\link[ComplexHeatmap]{ComplexHeatmap-package}}. Creates expression or
#'   correlation heatmaps from \code{eSet} or \code{MSnSet} objects.
#'
#' @param eset an object of class \code{\link[Biobase]{eSet}} or
#'   \code{\link[MSnbase]{MSnSet-class}}.
#' @param clustering_distance character; the distance measure used to cluster
#'   rows and columns. Passed to \code{\link[stats]{dist}}. One of
#'   (\code{"euclidean"}, \code{"maximum"}, \code{"manhattan"},
#'   \code{"canberra"}, \code{"binary"}, \code{"minkowski"}, \code{"pearson"},
#'   \code{"spearman"}, or \code{"kendall"}). Default is \code{"euclidean"}.
#' @param clustering_method character; the agglomeration method used to cluster
#'   rows and columns. Passed to \code{\link[stats]{hclust}}. One of
#'   (\code{"ward.D"}, \code{"ward.D2"}, \code{"single"}, \code{"average"},
#'   \code{"complete"}, \code{"mcquitty"}, \code{"median"}, or
#'   \code{"centroid"}). Default is \code{"ward.D"}.
#' @param heatmap_type character; the type of heatmap to generate. Must be (an
#'   abbreviation of) either \code{"expression"}, \code{"sample_correlation"},
#'   or \code{"feature_correlation"}. The default (\code{"expression"})
#'   generates a heatmap with features as rows and samples as columns. The other
#'   options calculate the sample or feature correlation matrix and generate a
#'   heatmap with samples or features as both rows and columns.
#' @param cor_method character; the method used to generate the correlation
#'   matrix if \code{heatmap_type} is a correlation heatmap. One of
#'   (\code{"pearson"}, \code{"kendall"}, or \code{"spearman"}). Defaults to
#'   \code{"pearson"}. Passed to \code{\link[stats]{cor}}.
#' @param cluster_columns logical; whether to cluster the columns.
#' @param cluster_rows logical; whether to cluster the rows.
#' @param show_column_dendrogram logical; whether to show the column dendrogram.
#'   Does not affect clustering.
#' @param show_row_dendrogram similar to \code{show_column_dendrogram}, but for
#'   rows.
#' @param show_column_names logical; whether to show the column names.
#' @param show_row_names logical; whether to show the row names.
#' @param heatmap_title character; overall title for the heatmap.
#' @param heatmap_legend_title character; title of the heatmap color legend.
#' @param color_range numeric; vector of length 2 used to restrict the colors of
#'   the heatmap body. Useful when color differences are not easily discernible
#'   due to outliers.
#' @param heatmap_args list of arguments passed to
#'   \code{\link[ComplexHeatmap]{Heatmap}}.
#' @param anno_column character; one or more column names of \code{pData(eset)}
#'   used to annotate the columns of the heatmap. By default, columns are not
#'   annotated.
#' @param anno_row character; one or more column names of \code{fData(eset)}
#'   used to annotate the rows of the heatmap. By default, rows are not
#'   annotated.
#' @param anno_column_titles character; names for the column annotations, if
#'   different from \code{anno_column}. Must be the same length as
#'   \code{anno_column}.
#' @param anno_row_titles similar to \code{anno_column_titles}, but for row
#'   annotations.
#' @param anno_column_colors list of custom colors for column annotations
#'   specified by \code{anno_column}. List names must match one or more names in
#'   \code{anno_column}. If modifying the colors of a continuous column
#'   annotation, **always use \code{\link[circlize]{colorRamp2}}** with desired
#'   breaks and colors. Otherwise, pass a vector of unique colors (order
#'   matters).
#' @param anno_row_colors same as \code{anno_column_colors}, but for
#'   \code{anno_row}.
#' @param anno_args list of arguments passed to
#'   \code{\link[ComplexHeatmap]{HeatmapAnnotation}}.
#' @param filename character; file name used to save the heatmap. Must end in
#'   (".png", ".bmp", ".jpg", ".tif", or ".pdf").
#' @param height numeric; height of heatmap in inches. Default is \code{5}.
#' @param width numeric; width of heatmap in inches. Default is \code{5}.
#' @param draw_args list of arguments passed to
#'   \code{\link[ComplexHeatmap]{draw-HeatmapList-method}}. Unlikely to be used
#'   often.
#'
#' @return An object of class \code{\link[ComplexHeatmap]{HeatmapList-class}} or
#'   nothing (save to file instead).
#'
#' @seealso
#' \href{https://jokergoo.github.io/ComplexHeatmap-reference/book/}{ComplexHeatmap
#' Complete Reference}
#'
#' @note If the `hclust` error `NA/NaN/Inf in foreign function call` occurs, it
#'   means there are cases where two or more features are not present in the
#'   same samples, so their distances will be `NA`, and it is impossible to
#'   perform hierarchical clustering. If this happens, filter to features
#'   present in more than 50% of samples with `m[rowMeans(!is.na(m)) > 0.5, ]`,
#'   where `m` is the `MSnSet`.
#'
#' @references Gu, Z., Eils, R., & Schlesner, M. (2016). Complex heatmaps reveal
#'   patterns and correlations in multidimensional genomic data.
#'   *Bioinformatics*, *32*(18), 2847-2849.
#'   \url{https://doi.org/10.1093/bioinformatics/btw313}
#'
#'   Gu, Z., Gu, L., Eils, R., Schlesner, M., & Brors, B. (2014). Circlize
#'   implements and enhances circular visualization in R. *Bioinformatics*,
#'   *30*(19), 2811-2812. \url{https://doi.org/10.1093/bioinformatics/btu393}
#'
#' @md
#'
#' @import ComplexHeatmap
#' @importFrom dplyr %>% mutate across where
#' @importFrom Biobase exprs pData fData
#' @importFrom stats cor setNames
#' @importFrom grDevices bmp jpeg png tiff pdf dev.off
#' @importFrom purrr map2
#' @importFrom circlize colorRamp2
#' @importFrom utils modifyList
#'
#' @export complex_heatmap
#'
#' @examples
#' # library(MSnSet.utils)
#' data(longitudinal_biomarker_study) # MSnSet
#' ee <- longitudinal_biomarker_study
#' # Missingness filter - clustering may fail otherwise
#' ee <- ee[rowMeans(!is.na(exprs(ee))) >= 0.5, ]
#'
#' # Expression heatmap
#' complex_heatmap(ee)
#'
#' # Limit color range to see differences more easily
#' complex_heatmap(ee, color_range = c(-2, 2))
#'
#' # Sample correlation heatmap
#' complex_heatmap(ee, heatmap_type = "s")
#'
#' # Annotate columns by "Type" (categorical) and "Age" (continuous)
#' # Annotate rows by "isSpike" (logical)
#' complex_heatmap(ee, anno_column = c("Type", "Age"),
#'                 anno_row = "isSpike")
#'


# Generic heatmap function ----
complex_heatmap <- function(
  eset, # MSnSet object
  clustering_distance = c("euclidean", "maximum", "manhattan", "canberra",
                          "binary", "minkowski", "pearson", "spearman",
                          "kendall"),
  clustering_method = c("ward.D", "ward.D2", "single", "average", "complete",
                        "mcquitty", "median", "centroid"),
  heatmap_type = "expression",
  cor_method = c("pearson", "kendall", "spearman"),
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_column_dendrogram = TRUE,
  show_row_dendrogram = TRUE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  heatmap_title = character(0),
  heatmap_legend_title = NULL,
  color_range = NULL,
  heatmap_args = list(),
  anno_column = NULL,
  anno_row = NULL,
  anno_column_titles = anno_column,
  anno_row_titles = anno_row,
  anno_column_colors = list(),
  anno_row_colors = list(),
  anno_args = list(),
  filename = NULL,
  height = 5,
  width = 5,
  draw_args = list()
) {

  ## Check input --------------------------------------------------
  # Match arguments
  clustering_distance <- match.arg(clustering_distance)
  clustering_method <- match.arg(clustering_method)
  cor_method <- match.arg(cor_method)

  # Check for pData
  if (!is.null(anno_column)) {
    if (ncol(pData(eset)) == 0) {
      stop("anno_column is provided, but pData(eset) is empty")
    }

    pData(eset) <- mutate(pData(eset),
                          across(.cols = where(is.character), as.factor))
  }
  # check for fData
  if (!is.null(anno_row)) {
    if (ncol(fData(eset)) == 0) {
      stop("anno_row is provided, but fData(eset) is empty")
    }

    fData(eset) <- mutate(fData(eset),
                          across(.cols = where(is.character), as.factor))
  }

  # Check that args are lists
  if (any(!unlist(lapply(list(heatmap_args, anno_row_colors,
                              anno_column_colors,
                              anno_args, draw_args),
                         is.list)))) {
    stop(paste("heatmap_args, anno_row_colors, anno_col_colors, anno_args,",
               "and draw_args must be lists"))
  }

  # Check color_range
  if (!is.null(color_range)) {
    if(!is.numeric(color_range) | length(color_range) != 2) {
      stop(paste("color_range must be a numeric vector of length 2",
                 "specifying lower and upper bounds for heatmap colors"))
    }
  }

  # Check width and height
  if (!is.numeric(height) | !is.numeric(width)) {
    stop(paste("width and height must be numbers specifying the",
               "width and height of the plot in inches"))
  }

  # Check that all arguments passed as lists are valid
  if (!all(names(heatmap_args) %in% names(formals(Heatmap)))) {
    stop(paste("One or more elements of heatmap_args is not",
               "an argument of ComplexHeatmap::Heatmap"))
  }
  if (!all(names(anno_args) %in% names(formals(HeatmapAnnotation)))) {
    stop(paste("One or more elements of anno_args is not",
               "an argument of ComplexHeatmap::HeatmapAnnotation"))
  }
  # Not sure how to check draw_args because draw is a method

  # Check that the number of annotation titles is the same
  # as the number of annotations
  if(length(anno_column) != length(anno_column_titles)) {
    stop(paste("length of anno_column is not the same as",
               "length of anno_column_titles"))
  }
  if(length(anno_row) != length(anno_row_titles)) {
    stop("length of anno_row is not the same as length of anno_row_titles")
  }


  # Data for heatmap ----------------------------------------------
  pttrn <- paste0("^", heatmap_type)

  if (grepl(pttrn, "expression")) {
    x <- exprs(eset)
    y <- pData(eset)
    z <- fData(eset)
  } else if (grepl(pttrn, "sample_correlation")) {
    x <- cor(exprs(eset), exprs(eset),
             use = "pairwise.complete.obs",
             method = cor_method)
    y <- z <- pData(eset)
    anno_row <- anno_column # row annotation same as column annotation
  } else if (grepl(pttrn, "feature_correlation")) {
    x <- cor(t(exprs(eset)), t(exprs(eset)),
             use = "pairwise.complete.obs",
             method = cor_method)
    y <- z <- fData(eset)
    anno_column <- anno_row # column annotation same as row annotation
  } else {
    stop(paste("heatmap_type must be (an abbreviation of) either 'expression',",
               "'sample_correlation', or 'feature_correlation'"))
  }


  # Row and column annotations ------------------------------------
  if (!is.null(anno_column)) {
    # Use pData for columns
    column_anno <- annotate_heatmap(y, anno = anno_column,
                                    anno_names = anno_column_titles,
                                    anno_colors = anno_column_colors,
                                    choice = "column", anno_args = anno_args)
  } else {
    column_anno <- NULL # No annotation
  }
  if (!is.null(anno_row)) {
    # Use fData for rows
    row_anno <- annotate_heatmap(z, anno = anno_row,
                                 anno_names = anno_row_titles,
                                 anno_colors = anno_row_colors,
                                 choice = "row", anno_args = anno_args)
  } else {
    row_anno <- NULL # No annotation
  }


  # Colors and color breaks for the heatmap body ------------------
  breaks_and_colors <- get_color_breaks(range(x, na.rm = TRUE),
                                        user_breaks = color_range)

  breaks <- breaks_and_colors[["breaks"]]
  colors <- breaks_and_colors[["colors"]]

  # Color function for heatmap body
  col_fun <- circlize::colorRamp2(breaks = breaks, colors = colors)

  # Create heatmap ------------------------------------------------
  # Arguments that will be passed to ComplexHeatmap::Heatmap
  heatmap_args <- list(matrix = x,
                       col = col_fun,
                       column_title = heatmap_title,
                       border = TRUE, # border around heatmap body
                       # do not trim row names
                       row_names_max_width = max_text_width(rownames(x)),
                       # do not trim column names
                       column_names_max_height = max_text_width(colnames(x)),
                       clustering_distance_rows = clustering_distance,
                       clustering_distance_columns = clustering_distance,
                       clustering_method_rows = clustering_method,
                       clustering_method_columns = clustering_method,
                       cluster_columns = cluster_columns,
                       cluster_rows = cluster_rows,
                       show_row_dend = show_row_dendrogram,
                       show_column_dend = show_column_dendrogram,
                       show_column_names = show_column_names,
                       show_row_names = show_row_names,
                       top_annotation = column_anno,
                       left_annotation = row_anno,
                       heatmap_legend_param = list(
                         title = heatmap_legend_title
                       )
  ) %>%
    modifyList(val = heatmap_args, keep.null = TRUE) # update arguments

  # Create heatmap
  ht <- do.call(what = Heatmap, args = heatmap_args)


  # Save or display heatmap ----------------------------------------
  if (!is.null(filename)) {
    # Begin saving
    save_plot(filename = filename, height = height, width = width)
    # After plot is created, run dev.off
    on.exit(expr = dev.off())
  }
  # Create heatmap
  c(modifyList(list(object = ht,
                    merge_legends = TRUE),
               val = draw_args, keep.null = TRUE)) %>%
    {do.call(what = draw, args = .)}
}

utils::globalVariables(".")


## Helper functions ------------------------------------------------



# Get appropriate color breaks for heatmap body ----
# Set color_range for unique colors within full range of values
get_color_breaks <- function(x, user_breaks = NULL) {

  # If user_breaks is not provided, set to x
  if (is.null(user_breaks)) {
    user_breaks <- x
  }
  # Extend or reduce range for unique colors
  lower <- sort(unique(c(x[1], user_breaks[1])))
  n_lower <- length(lower) # 1 or 2
  upper <- sort(unique(c(x[2], user_breaks[2])))
  n_upper <- length(upper) # 1 or 2

  if (all(sign(x) == +1 | sign(x) == 0)) { # If mix of positive and 0
    # at most 3 breaks and colors
    breaks <- c(0, upper)
    colors <- c("white", rep("red", n_upper))
  } else if (all(sign(x) == -1 | sign(x) == 0)) { # If mix of negative and 0
    # at most 3 breaks and colors
    breaks <- c(lower, 0)
    colors <- c(rep("blue", n_upper), "white")
  } else { # If mix of positive and negative
    # at most 5 breaks and colors
    breaks <- c(lower, 0, upper)
    colors <- c(rep("blue", n_lower), "white",
                rep("red", n_upper))
  }

  return(list(breaks = breaks, colors = colors))
}


# Get colors for each column or row annotation ----
# Use MSnSet.utils::jet.colors for discrete variables and
# scales::viridis_pal for continuous variables.
#'
#' @importFrom scales viridis_pal
#'
get_anno_colors <- function(x) {
  # If x is numeric, use colorRamp2. Otherwise, use a discrete
  # color palette with length equal to the number of unique groups.
  if (is.numeric(x)) {
    seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 5) %>%
      {circlize::colorRamp2(breaks = .,
                            colors = rev(scales::viridis_pal()(length(.))))}
  } else if (is.factor(x)) {
    jet2.colors(nlevels(x))
  } else {
    jet2.colors(length(unique(x)))
  }

}

utils::globalVariables(".")


# Name annotation colors ----
# x is a vector of colors or a color function from get_anno_colors
# or supplied by the user. y is the values from the data
set_anno_names <- function(x, y) {
  if (!is.function(x)) {
    if (is.factor(y)) {
      names(x) <- levels(y)
    } else { # character, logical
      y <- unique(y)
      y <- y[!is.na(y)]
      names(x) <- y
    }
  }

  return(x)
}


# Create row or column annotations ----
annotate_heatmap <- function(y, anno, anno_names = anno,
                             choice = c("column", "row"),
                             anno_colors = list(),
                             anno_args = list()) {
  # Check input
  choice <- match.arg(choice)

  # Data frame of annotation values
  anno_df <- data.frame(y[, anno]) %>%
    setNames(anno)

  # List of annotation colors
  anno_col <- lapply(as.list(anno_df), get_anno_colors) %>%
    # Replace colors if provided by user
    modifyList(val = anno_colors, keep.null = TRUE) %>%
    # Set names on colors
    map2(.y = as.list(anno_df), .f = set_anno_names)


  # Arguments passed to HeatmapAnnotation
  anno_args <- list(df = anno_df,
                    col = anno_col,
                    annotation_label = anno_names,
                    border = FALSE,
                    which = choice,
                    # gp = gpar(col = "black"),
                    annotation_legend_param = list(
                      border = TRUE
                    )
  ) %>%
    modifyList(val = anno_args, keep.null = TRUE) # Update arguments

  # Create HeatmapAnnotation object
  return(do.call(what = HeatmapAnnotation, args = anno_args))
}


# Save non-ggplot objects ----
save_plot <- function(filename, height, width) {
  # File extension
  file_ext <- gsub(".*\\.([a-z])", "\\1", filename)

  # Check file extension
  file_ext <- match.arg(file_ext,
                        choices = c("png", "bmp", "jpg", "tif", "pdf"))

  # Function to save plot
  save_func <- switch(file_ext,
                      "png" = png,
                      "bmp" = bmp,
                      "jpg" = jpeg,
                      "tif" = tiff,
                      "pdf" = pdf)

  message(sprintf("Saving %g x %g inch heatmap as %s",
                  height, width, filename))

  # arguments passed to save_func
  save_args <- list(filename = filename, file = filename,
                    width = width, height = height,
                    res = 300, dpi = 300,
                    units = "in", compression = "lzw") %>%
    # Subset to relevant arguments
    .[names(.) %in% names(formals(save_func))]

  # Begin saving
  do.call(what = save_func, args = save_args)
}

