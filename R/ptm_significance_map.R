#' @title PTM significance map
#'
#' @description Mapping PTMs onto the protein sequence with colors
#'   indicating differential abundance and statistical significance. Data must be
#'   filtered to peptides modified with a single PTM in advance.
#' @param x `data.frame` with required columns "accession", "gene", "site",
#'   as well as "adj.P.Val" and "logFC". Results from additional statistical
#'   models may be used by appending adj.P.Val and logFC with "_ModelNameHere".
#'   For example: adj.P.Val_model1
#'   When using multiple models, desired models must be passed to `p_val_from`.
#'   Recommended to provide the `fData` table of an MSnSet object.
#' @param accession character; name of the protein (or accession, in general
#'   terms) for which to plot peptides. Defaults to "all" which will plot
#'   all proteins with at least one feature with adj.P.Val < `alpha`.
#' @param p_val_from character vector; names appending adj.P.Val & logFC
#'   in x to differentiate statistics obtained from different models. When
#'   p_val_from == NULL (default) it will ba assumed that only one model is
#'   present and the column names are "adj.P.Val" and "logFC".
#' @param alpha numeric; desired p.Value cutoff
#'
#' @return \code{ggplot} object.
#'
#' @importFrom dplyr select filter mutate %>% distinct case_when
#' @importFrom tidyr separate_rows
#' @importFrom tidyselect contains
#' @importFrom stringr str_detect
#' @importFrom Biostrings readAAStringSet
#' @importFrom ggplot2 ggplot geom_rect geom_linerange scale_color_manual ylim theme_minimal theme xlab ggtitle element_text element_blank
#' @importFrom ggrepel geom_label_repel
#' @importFrom glue glue
#' @importFrom BiocGenerics width
#'
#' @export ptm_significance_map
#'
#' @md
#'
#' @examples
#' Need example data

ptm_significance_map <- function(x,
                                 fasta_path,
                                 accession = "all",
                                 p_val_from = NULL,
                                 alpha = 0.05)
{
  if (!any(str_detect(colnames(x), "adj\\.P\\.Val"))) {
    stop("missing adj.P.Val")
  }
  if (!any(str_detect(colnames(x), "logFC"))) {
    stop("missing logFC")
  }
  if (!any(str_detect(colnames(x), "gene"))) {
    stop("missing gene")
  }
  if (!any(str_detect(colnames(x), "accession"))) {
    stop("missing accession")
  }
  if (!any(str_detect(colnames(x), "site"))) {
    stop("missing site")
  }

  # Correct? - BiocGenerics width
  # change p_val_from default to use all models if present?

  # Get protein sequences
  fst <- readAAStringSet(fasta_path, format="fasta", nrec=-1L, skip=0L, use.names=TRUE)
  uniprot_acc_iso <- sub("sp\\|([^|]*)\\|.*","\\1", names(fst))
  names(fst) <- uniprot_acc_iso


  if (is.null(p_val_from)) {
    p_val_from_floop <- "NULL"
  } else {
    p_val_from_floop <- p_val_from
  }

  for (j in p_val_from_floop) {

    if (is.null(p_val_from)) {

      df <- x

    } else {

      df <- x %>%
        select(c(gene, accession, site), contains(p_val_from))

      colnames(df) <- sub(glue("_{j}"), "", colnames(df))

    }

    if (accession == "all"){

      df_accession <- df %>%
        dplyr::filter(adj.P.Val <= alpha)

      accession_i <- unique(df_accession$accession)

    } else {

      if (!is.character(accession)) {
        stop("accession must be character")
      }
      accession_i <- accession
    }

    for (i in accession_i) {

      y <- filter(df, accession == i) %>%
        distinct() %>%
        mutate(positions = sub("^.*_(.*)", "\\1", site),
               positions = gsub("[[:lower:]]","", positions, perl = TRUE)) %>%
        separate_rows(positions, sep = ";") %>%
        mutate(num = gsub("[[:alpha:]]", "", positions),
               aas = gsub("\\d", "", positions),
               label = sub("^.*_(.*)", "\\1", site))

      gene_i <- unique(y[y$accession == i, "gene"])

      # unlisting index
      idx <- unlist(lapply(seq_along(y$num), function(.) rep(., length(y$num[[.]]))))

      z <- data.frame(label = unlist(y$label),
                      pos = as.numeric(unlist(y$num)),
                      aas = unlist(y$aas),
                      adj.P.Val = y$adj.P.Val[idx],
                      logFC = y$logFC[idx],
                      stringsAsFactors = FALSE) %>%
        mutate(change = case_when(
          logFC > 0 & adj.P.Val < alpha ~ "up",
          logFC < 0 & adj.P.Val < alpha ~ "down",
          adj.P.Val > alpha ~ "steady"),
          change = ordered(change, levels=c("up","steady","down")))

      d <- data.frame(x1=1, x2=width(fst[accession_i]), y1=-1, y2=+1)

      p <- ggplot(data = z) +
        geom_rect(data = d,
                  mapping = aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
                  fill="lightgrey", show.legend = F) +
        geom_linerange(mapping = aes(ymin=-1,ymax=+1, x=pos, color=change),
                       linewidth = 1) +
        scale_color_manual(values = c("up"="red", "down"="blue", "steady"="black")) +
        geom_label_repel(aes(x=pos, y=1, label=label, color=change),
                         max.overlaps = 100,
                         nudge_y = 10,
                         ylim = c(2,20),
                         segment.alpha = 0.2,
                         show.legend = FALSE,
                         seed = 0,
                         fill="white",
                         box.padding = 1,
                         point.padding = 0) +
        ylim(-1,+20) +
        theme_minimal() +
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              plot.title = element_text(hjust=0.5, size=20),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
        xlab("amino acid position")
      if (!is.null(p_val_from)) {
        p <- p + ggtitle(glue("{gene_i}, {i} ({j})"))
      } else {
        p <- p + ggtitle(glue("{gene_i}, {i}"))
      }
      print(p)
    }
  }
}

utils::globalVariables(
  c("accession", "gene", "site")
)