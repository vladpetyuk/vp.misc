#' @title Peptide Correlation Analysis (PeCorA)
#'
#' @description Perform Peptide Correlation Analysis (PeCorA) as described in
#'   \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8592057/}.
#'
#' @param m \code{MSnSet} object with pData and fData tables. fData must have
#'   columns "Protein" and "Peptide".
#' @param treatment character; the name of the pData column containing treatment
#'   information. The data can be numeric or categorical.
#' @param proteins character; one or more proteins in
#'   \code{fData(m)[["Protein"]]} to analyze.
#' @param median_mod logical; whether to use the median version of PeCorA. See
#'   \code{\link[MSnSet.utils]{pecora_mod}} for details.
#'
#' @details The MSnSet provided to \code{m} should be normalized as described in
#'   the PeCorA publication. See \code{\link[MSnSet.utils]{pecora_preprocess}}
#'   for details. This function was based on the \code{PeCorA_analysis} function
#'   from the PeCorA R package
#'   (\url{https://github.com/jessegmeyerlab/PeCorA}).
#'
#' @return An \code{MSnSet} where the fData table includes new or updated
#'   columns "pecora_pval_{treatment}" and "pecora_adj_pval_{treatment}".
#'
#' @seealso \code{\link[MSnSet.utils]{pecora_preprocess}},
#'   \code{\link[MSnSet.utils]{pecora_mod}},
#'   \code{\link[MSnSet.utils]{plot_pecora}}
#'
#' @references Dermit, M., Peters-Clarke, T. M., Shishkova, E., & Meyer, J. G.
#'   (2021). Peptide Correlation Analysis (PeCorA) Reveals Differential
#'   Proteoform Regulation. \emph{Journal of proteome research, 20(4)},
#'   1972–1980. \url{https://doi.org/10.1021/acs.jproteome.0c00602}
#'
#' @importFrom dplyr select mutate filter rename %>% distinct
#' @importFrom rlang sym
#' @importFrom tidyr pivot_longer any_of
#'
#' @export pecora_analysis

pecora_analysis <- function(m,
                            treatment,
                            proteins,
                            median_mod = FALSE)
{
  f_data <- fData(m)

  ## Must have "Protein" and "Peptide" columns
  for (col_i in c("Protein", "Peptide")) {
    if (!col_i %in% colnames(f_data)) {
      stop(sprintf("fData missing '%s' column.", col_i))
    }
  }

  pval_label <- paste0("pecora_pval_", treatment)
  padj_label <- paste0("pecora_adj_pval_", treatment)

  ## Remove previous results using the same treatment variable.
  f_data <- select(f_data, -any_of(c(pval_label, padj_label)))

  ## Filter to just those proteins for which we want to run PeCorA
  peptide_mapping <- f_data %>%
    filter(Protein %in% proteins) %>%
    select(Protein, Peptide)

  mat <- exprs(m)[rownames(peptide_mapping), ]

  metadata <- pData(m) %>%
    select(sym(treatment)) %>%
    mutate(Sample = rownames(.))

  if (is.character(metadata[[treatment]])) {
    metadata[[treatment]] <- as.factor(metadata[[treatment]])
  }

  # Prepare input and run PeCorA
  PeCorA_input <- mat %>%
    cbind(peptide_mapping) %>%
    as.data.frame() %>%
    pivot_longer(cols = -c(Peptide, Protein),
                 names_to = "Sample",
                 values_to = "LogRatio") %>%
    merge(metadata, by = "Sample") %>%
    dplyr::rename(Condition = sym(treatment))

  PeCorA_result <- pecora_mod(x = PeCorA_input, median_mod = median_mod)

  # Rename p-value and adjusted p-value columns
  pval_index <- which(colnames(PeCorA_result) == "pvalue")
  adj_pval_index <- which(colnames(PeCorA_result) == "adj_pval")
  colnames(PeCorA_result)[[pval_index]] <- pval_label
  colnames(PeCorA_result)[[adj_pval_index]] <- padj_label

  # Add PeCorA results to fData
  f_data <- f_data %>%
    mutate(feature_name = rownames(f_data)) %>%
    left_join(PeCorA_result, by = c("Peptide", "Protein"))
  rownames(f_data) <- f_data$feature_name

  f_data <- select(f_data, -feature_name)

  fData(m) <- f_data[rownames(exprs(m)), ]

  return(m)
}


utils::globalVariables(
  c("Protein", "Peptide")
)

