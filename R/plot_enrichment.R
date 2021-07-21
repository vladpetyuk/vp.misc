#' @title Dotplot of GO or Pathway Enrichment Results
#'
#' @description Create a dotplot of Gene Ontology or Pathway enrichment results
#' from \code{\link[ReactomePA]{enrichPathway}},
#' \code{\link[clusterProfiler]{enrichGO}},
#' or \code{\link[clusterProfiler]{enrichKEGG}}.
#'
#' @param x an \code{\link[DOSE]{enrichResult-class}} object.
#' @param num.categories integer; the number of pathways or gene ontologies to
#'     plot. The default selects the 10 terms with the highest counts, including
#'     ties.
#' @param sentence.case logical; whether to capitalize the first letter of each
#'     term. Certain exceptions such as "mRNA" will not be capitalized.
#'     Default is \code{TRUE}.
#' @param p.adjust.cutoff the significance cutoff for enrichment results. Only
#'     terms with adjusted p-values less than \code{p.adjust.cuttoff} will be
#'     considered for plotting.
#' @param label.wrap.chars the maximum number of characters per line for each
#'     term.
#' @param point.size the size of the plotted points. Default is 4.
#' @param gradient.trans the name of a transformation object for the color
#'     gradient, or the object itself. Default is \code{"log10"}. See the
#'     \code{trans} argument in \code{\link[ggplot2]{continuous_scale}} for
#'     more details.
#' @param ... additional arguments passed to
#'     \code{\link[ggplot2]{scale_color_gradient}}.
#'
#' @importFrom dplyr %>% arrange slice_max mutate filter rename
#' @importFrom scales label_wrap percent_format
#' @importFrom stringr str_to_sentence
#' @importFrom ggplot2 ggplot aes geom_point scale_x_continuous
#'             scale_y_discrete scale_color_gradient theme_minimal
#'
#' @export plot_enrichment
#'
#' @examples
#' # Load packages
#' library(clusterProfiler)
#' library(MSnSet.utils)
#'
#' # Data from clusterProfiler
#' data(gcSample)
#'
#' # GO: BP Enrichment results
#' x <- enrichGO(gene = gcSample[[6]],
#'               OrgDb = "org.Hs.eg.db",
#'               ont = "BP",
#'               pvalueCutoff = 1,
#'               universe = unique(unlist((gcSample))),
#'               qvalueCutoff = 1,
#'               minGSSize = 0)
#'
#' # For the purposes of this example, the terms will not be
#' # filtered based on their adjusted p-values.
#' plot_enrichment(x, p.adjust.cutoff = 1,
#'                 num.categories = 6, label.wrap.chars = 40)


plot_enrichment <- function(x,
                            num.categories = 10,
                            sentence.case = TRUE,
                            p.adjust.cutoff = 0.05,
                            label.wrap.chars = Inf,
                            point.size = 4,
                            gradient.trans = "log10",
                            ...) {
  # Check validity of x
  if (class(x) != "enrichResult") {
    stop(paste0(
      "x must be an object of class `enrichResult`, not ",
      "`",
      class(x),
      "`"
    ))
  }

  # Get the results dataframe.
  x <- as.data.frame(x)

  # Arrange in descending order by count.
  # Select the n descriptions with the highest counts (including ties).
  x <- x %>%
    filter(p.adjust < p.adjust.cutoff) %>%
    rename(TermGenes = Count) %>%
    slice_max(order_by = TermGenes, n = num.categories) %>%
    mutate(
      TotalGenes = as.numeric(gsub(".*\\/(.*)", "\\1", GeneRatio)),
      BgTerm = as.numeric(gsub("(.*)\\/.*", "\\1", BgRatio)),
      BgTotal = as.numeric(gsub(".*\\/(.*)", "\\1", BgRatio)),
      GeneProp = TermGenes / TotalGenes,
      BgProp = BgTerm / BgTotal
    )

  if (nrow(x) == 0) {
    stop("No terms pass the significance threshold set by `p.adjust.cutoff`")
  }

  if (sentence.case) {
    # If a sentence begins with these exceptions, do not capitalize them
    exceptions <- c("mRNA")
    exceptions <- paste0("^", paste(exceptions, collapse = "|^"))
    x <- x %>%
      mutate(Description = ifelse(
        grepl(exceptions, Description),
        Description,
        str_to_sentence(Description)
      ))
  }

  # Convert Description to factor to preserve order when plotting
  x <- x %>%
    mutate(Description = factor(Description, levels = rev(Description)))

  p <- ggplot(x) +
    geom_point(aes(x = GeneProp, y = Description, color = p.adjust),
               size = point.size) +
    scale_x_continuous(labels = percent_format(drop0trailing = TRUE)) +
    scale_y_discrete(name = NULL,
                     labels = label_wrap(label.wrap.chars)) +
    scale_color_gradient(name = "Adjusted\np-value",
                         trans = gradient.trans,
                         ...) +
    theme_minimal()

  return(p)
}
