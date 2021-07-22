#' @title Dotplot of GO or Pathway Enrichment Results
#'
#' @description Create a dotplot of Gene Ontology or Pathway enrichment results
#' from \code{\link[ReactomePA]{enrichPathway}},
#' \code{\link[clusterProfiler]{enrichGO}},
#' or \code{\link[clusterProfiler]{enrichKEGG}}.
#'
#' @param x an \code{\link[DOSE]{enrichResult-class}} object.
#' @param num.categories integer; the number of terms to plot. The default
#' selects the 10 terms with the highest counts, including ties.
#' @param sentence.case logical; whether to capitalize the first letter of each
#'     term. Certain exceptions such as "mRNA" will not be capitalized.
#'     Default is \code{TRUE}.
#' @param p.adjust.cutoff the significance cutoff for enrichment results. Only
#'     terms with adjusted p-values less than \code{p.adjust.cuttoff} will be
#'     considered for plotting.
#' @param label.wrap.chars the maximum number of characters per line for each
#'     term.
#' @param point.size the size of the plotted points. Default is 4.
#' @param ... additional arguments passed to
#'     \code{\link[ggplot2]{scale_color_gradient}}.
#'
#' @importFrom dplyr %>% arrange slice_max mutate filter rename
#' @importFrom scales label_wrap percent_format scientific pretty_breaks
#'             trans_breaks
#' @importFrom stringr str_to_sentence
#' @importFrom ggplot2 ggplot aes geom_point scale_x_continuous
#'             scale_y_discrete scale_color_gradient theme_bw
#'             theme element_line element_blank element_rect
#'
#' @export plot_enrichment
#'
#' @examples
#' library(clusterProfiler)
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
#'                 num.categories = 6, label.wrap.chars = 35,
#'                 low = "darkred", high = "grey85")


plot_enrichment <- function(x,
                            num.categories = 10,
                            sentence.case = TRUE,
                            p.adjust.cutoff = 0.05,
                            label.wrap.chars = Inf,
                            point.size = 4,
                            ...) {
  # Check validity of x
  if (class(x) != "enrichResult") {
    stop(paste0(
      "x must be an object of class `enrichResult`, not ",
      "`", class(x), "`"
    ))
  }

  # Get the results dataframe.
  x <- as.data.frame(x)

  # Arrange in descending order by count.
  # Select the n descriptions with the highest counts (including ties).
  x <- x %>%
    filter(p.adjust < p.adjust.cutoff) %>%
    dplyr::rename(TermGenes = Count) %>%
    slice_max(order_by = TermGenes, n = num.categories) %>%
    mutate(
      TotalGenes = as.numeric(gsub(".*\\/(.*)", "\\1", GeneRatio)),
      BgTerm = as.numeric(gsub("(.*)\\/.*", "\\1", BgRatio)),
      BgTotal = as.numeric(gsub(".*\\/(.*)", "\\1", BgRatio)),
      GeneProp = TermGenes / TotalGenes,
      BgProp = BgTerm / BgTotal
    )

  # Capitalize first letter of each term
  if (sentence.case) {
    # If a term begins with these exceptions, do not capitalize them.
    # Add more exceptions to the vector, as needed.
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

  x_lims <- c(floor(min(x$GeneProp) * 100),
              ceiling(max(x$GeneProp) * 100)) / 100

  # Determine whether the identity or log10 transformation should be used.
  # Also determine the gradient label breaks.
  # If the powers are the same, use the identity transformation
  if (length(unique(gsub(".*e\\-(\\d{+})",
                         "\\1", scientific(x$p.adjust)))) == 1) {
    gradient.trans <- "identity"
    gradient.breaks <- pretty_breaks(n = 4)
  } else {
    gradient.trans <- "log10"
    gradient.breaks <-
      trans_breaks("log10", function(x)
        10 ^ x, n = 4)
  }

  # Plot
  p <- ggplot(x) +
    geom_point(aes(x = GeneProp, y = Description, color = p.adjust),
               size = point.size) +
    scale_x_continuous(labels = percent_format(drop0trailing = TRUE),
                       limits = x_lims) +
    scale_y_discrete(name = NULL,
                     labels = label_wrap(label.wrap.chars)) +
    scale_color_gradient(name = "Adjusted\np-value",
                         trans = gradient.trans,
                         breaks = gradient.breaks,
                         ...) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      panel.border = element_rect(color = NA),
      axis.line = element_line(color = "black")
    )

  return(p)
}


