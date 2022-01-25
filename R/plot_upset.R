#' @title Labeled UpSet Plot
#'
#' @description A wrapper around functions from
#'   \code{\link[ComplexHeatmap]{ComplexHeatmap-package}}. Creates an
#'   \href{https://doi.org/10.1109/TVCG.2014.2346248}{UpSet plot} with labeled
#'   bars to visualize relationships between sets. UpSet plots are preferred to
#'   Venn diagrams, especially when the number of intersections is large.
#'
#' @param ... Input sets (a list is recommended). See
#'   \code{\link[ComplexHeatmap]{make_comb_mat}} for details.
#' @param mode The mode for forming the combination set, see Mode section of
#'   \code{\link[ComplexHeatmap]{make_comb_mat}}.
#' @param top_n_comb Number of largest intersections to display. Defaults to all
#'   intersections.
#' @param scale_set_bars logical; whether to scale set size bars so that their
#'   lengths are proportional to the lengths of the intersection size bars.
#' @param row_labels vector of labels for the rows. Defaults to the names of the
#'   input sets.
#' @param row_names_gp graphical parameters for row labels specified by
#'   \code{\link[grid]{gpar}}.
#' @param annotation_name_gp graphical parameters for barplot titles specified
#'   by \code{\link[grid]{gpar}}.
#' @param heatmap_args list of additional arguments passed to
#'   \code{\link[ComplexHeatmap]{Heatmap}}.
#' @param cell_dim \code{\link[grid]{unit}} object used for both the height and
#'   width of each heatmap cell.
#' @param extend numeric; multiply y-axis upper limit by (1 + `extend`).
#'   Increase to prevent clipping of barplot labels and overall title.
#' @param rot numeric \[-360, +360\]; angle to rotate column barplot
#'   labels.
#' @param vjust numeric \[0, 1\]; vertical justification of column barplot
#'   labels.
#' @param hjust numeric \[0, 1\]; horizontal justification of column barplot
#'   labels.
#' @param bar_label_gp graphical parameters for barplot labels specified by
#'   \code{\link[grid]{gpar}}.
#' @param filename character; name of the file to save the plot. Must end in
#'   ".png", ".bmp", ".jpg", ".tif", or ".pdf". If not provided, the UpSet plot
#'   will be displayed instead of saved.
#' @param height numeric; height of the entire plot in inches.
#' @param width numeric; width of the entire plot in inches.
#'
#' @return An object of class \code{\link[ComplexHeatmap]{HeatmapList-class}} or
#'   nothing (save to file instead).
#'
#' @note \code{\link[base]{split}} is useful for creating input lists.
#'
#' @seealso
#' \href{https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html}{ComplexHeatmap
#' Complete Reference: UpSet Plot}
#'
#' @references Gu, Z., Eils, R., & Schlesner, M. (2016). Complex heatmaps reveal
#'   patterns and correlations in multidimensional genomic data.
#'   *Bioinformatics*, *32*(18), 2847-2849.
#'   \url{https://doi.org/10.1093/bioinformatics/btw313}
#'
#'   Lex, A., Gehlenborg, N., Strobelt, H., Vuillemot, R., & Pfister, H. (2014).
#'   UpSet: Visualization of Intersecting Sets. *IEEE transactions on
#'   visualization and computer graphics*, *20*(12), 1983-1992.
#'   \url{https://doi.org/10.1109/TVCG.2014.2346248}
#'
#' @md
#'
#' @import ComplexHeatmap
#' @import grid
#' @importFrom grDevices bmp jpeg png tiff pdf dev.off
#' @importFrom utils modifyList
#'
#' @export plot_upset
#'
#' @examples
#' # Input list
#' set.seed(99)
#' x <- list(A = sample(letters, 5),
#'           B = sample(letters, 10),
#'           C = sample(letters, 15),
#'           D = sample(letters, 20))
#'
#' # Base plot
#' plot_upset(x)
#'
#' # Change label rotation and center horizontally
#' plot_upset(x, rot = 0, hjust = 0.5)
#'
#' # Top 2 largest intersections
#' plot_upset(x, top_n_comb = 2)
#'
#' # Add overall title
#' plot_upset(x, heatmap_args = list(column_title = "Title"))
#'
#' \dontrun{
#' # Save to 3.5" by 7" PDF
#' plot_upset(x, filename = "upset_plot.pdf")
#' }
#'


plot_upset <- function(...,
                       mode = c("distinct", "intersect", "union"),
                       top_n_comb = Inf,
                       scale_set_bars = FALSE,
                       row_labels = rownames(m),
                       row_names_gp = gpar(fontsize = 10),
                       annotation_name_gp = gpar(fontsize = 10,
                                                 fontface = "bold"),
                       heatmap_args = list(),
                       cell_dim = unit(10, "mm"),
                       extend = 0.2, rot = 45, hjust = 0, vjust = 0,
                       bar_label_gp = gpar(fontsize = 10),
                       filename = character(0), height = 3.5, width = 7) {

  # Make combination matrix and calculate size of intersections and sets
  m <- make_comb_mat(..., mode = mode)
  cs <- comb_size(m) # size of combinations
  ss <- set_size(m) # size of full sets

  # Subset m to top_n_sets
  m <- m[, order(-cs)[1:min(top_n_comb, length(cs))]]
  cs <- comb_size(m) # recalculate cs (do not change ss)

  # Ratio of max set size to max combination size or 1
  bar_ratio <- ifelse(scale_set_bars, max(ss) / max(cs), 1)

  # NOTE: Eventually may need to add comb_order and set_order options
  row_order <- 1:length(ss)
  column_order <- order(-cs)

  # Arguments passed to ComplexHeatmap::Heatmap --------------------------------

  # NOTE:
  # I had to use Heatmap instead of UpSet because do.call failed when I
  # tried to modify the column barplots, despite me being able to modify the
  # row barplots without issue. Both row and column annotations are
  # S4 objects, but it only complained about the column annotation.
  # The layer_fun code is taken from UpSet.

  heatmap_args <- modifyList(
    x = list(
      matrix = m,
      row_labels = row_labels,
      row_names_gp = row_names_gp,
      cluster_rows = FALSE, cluster_columns = FALSE,
      row_names_side = "left",
      row_names_max_width = max_text_width(rownames(x)), # do not trim row names
      show_heatmap_legend = FALSE,
      row_order = row_order,
      column_order = column_order,
      height = length(ss) * cell_dim,
      width = length(cs) * cell_dim,
      rect_gp = gpar(type = "none"),
      layer_fun = function(j, i, x, y, w, h, fill) {
        n_row <- round(1 / as.numeric(h[1])) # num rows
        n_col <- round(1 / as.numeric(w[1])) # num columns
        subm <- matrix(pindex(m, i, j), nrow = n_row, byrow = FALSE)
        for (k in seq_len(n_row)) {
          # Even rows are grey, odd are white (no color)
          if (k %% 2) {
            grid.rect(y = k / n_row, height = 1 / n_row, just = "top",
                      gp = gpar(fill = "#f0f0f0", col = NA))
          }
        }
        # Black or light grey points
        grid.points(x, y, size = unit(min(3, as.numeric(cell_dim)), "mm"),
                    pch = 16,
                    gp = gpar(col = ifelse(pindex(m, i, j), "black", "#cccccc")))
        # For each column k, find indices of entries that are 1, and draw a
        # segment from the min to the max
        for (k in seq_len(n_col)) {
          # if (sum(subm[, k]) >= 2) { # Intersection of 2 or more groups
          i_range <- range(which(subm[, k] == 1)) # segment start and end
          grid.lines(c(k - 0.5, k - 0.5) / n_col, # center on column
                     (n_row - i_range + 0.5) / n_row, # center on row
                     gp = gpar(col = "black", lwd = 2))
          # }
        }
      },
      top_annotation = HeatmapAnnotation(
        "Intersection Size" = anno_barplot(
          x = cs,
          extend = extend,
          border = FALSE,
          gp = gpar(fill = "black"),
          height = unit(1.2, "in"),
          which = "column",
          axis = FALSE),
        # annotation_label changes the axis title
        annotation_name_gp = annotation_name_gp,
        annotation_name_side = "left",
        annotation_name_rot = 90),
      right_annotation = HeatmapAnnotation(
        "Set Size" = anno_barplot(
          x = ss,
          extend = extend,
          border = FALSE,
          gp = gpar(fill = "black"),
          width = unit(1.2, "in") * bar_ratio,
          which = "row",
          axis = FALSE),
        which = "row",
        annotation_name_gp = annotation_name_gp,
        annotation_name_side = "bottom",
        annotation_name_rot = 0
      )),
    val = heatmap_args, keep.null = TRUE) #  modify args

  ## Save to file if filename is provided --------------------------------------
  if (!identical(filename, character(0))) {
    save_plot(filename = filename, height = height, width = width)
    # After plot is created, run dev.off
    on.exit(expr = dev.off())
  }

  # UpSet plot -----------------------------------------------------------------
  ht <- do.call(what = Heatmap, args = heatmap_args)
  draw(ht)

  # Update row and column order if changed by user
  column_order <- heatmap_args$column_order
  row_order <- rev(heatmap_args$row_order) # reverse for labels

  # Add Intersection Size labels
  decorate_annotation(annotation = "Intersection Size", {
    grid.text(cs[column_order],
              x = seq_along(cs),
              y = unit(cs[column_order], "native") + unit(3, "pt"),
              vjust = vjust, hjust = hjust, rot = rot,
              default.units = "native", #just = c("left", "top"),
              gp = bar_label_gp)
  })

  # Add Set Size labels
  decorate_annotation(annotation = "Set Size", {
    grid.text(label = ss[row_order],
              x = unit(ss[row_order], "native") + unit(3, "pt"),
              y = seq(1 / (length(ss) * 2), 1 - 1 / (length(ss) * 2),
                      length.out = length(ss)),
              hjust = 0,
              gp = bar_label_gp)
  })
}



