#' @title Modified Peptide Correlation Analysis (PeCorA)
#'
#' @description This is a modified version of the PeCorA function. For the
#'   original paper see \url{https://doi.org/10.1021/acs.jproteome.0c00602}.
#'
#'   The modification involves aggregating the intensities of "all other
#'   peptides" within the test described in the PeCorA paper. We aggregate "all
#'   other peptides" by sample using the median. This way, each peptide mapping
#'   to a particular protein is not treated as an 'independent' observation. The
#'   effect is that p-values are much higher (less significant) than in the
#'   original version of PeCorA.
#'
#'   This median modification is turned off by default, and only affects the
#'   p-values and adjusted p-values.
#'
#' @param x PeCorA input table. See the original PeCorA function for details.
#' @param median_mod logical; whether to use the median version of PeCorA
#'   (\code{median_mod = TRUE}). Default is \code{FALSE}, which will use the
#'   original version.
#'
#' @details p-values are adjusted across all peptides for a given protein using
#'   the method of Benjamini and Hochberg.
#'
#'   This function was based on the \code{PeCorA} function from the PeCorA R
#'   package (\url{https://github.com/jessegmeyerlab/PeCorA}).
#'
#' @return A \code{data.frame} with columns "Protein", "Peptide", "pvalue", and
#'   "adj_pval".
#'
#' @importFrom dplyr %>% select distinct group_by n filter ungroup pull mutate
#' @importFrom purrr map
#' @importFrom tibble enframe
#' @importFrom tidyr unnest
#'
#' @noRd

pecora_test_all_proteins <- function(x,
                        median_mod = FALSE)
{
  proteins <- x %>%
    distinct(Protein, Peptide) %>%
    group_by(Protein) %>%
    filter(n() > 1) %>% # at least 2 peptides per protein
    ungroup() %>%
    pull(Protein) %>%
    unique()

  res <- x %>%
    filter(Protein %in% proteins) %>%
    split.data.frame(~ Protein) %>%
    map(pecora_test_single_protein,
        median_mod = median_mod,
        .progress = TRUE) %>%
    enframe(name = "Protein") %>%
    unnest(value) %>%
    group_by(Protein) %>%
    mutate(adj_pval = p.adjust(pvalue, method = "BH")) %>%
    ungroup()

  # Results summary
  signif_peptide <- res$adj_pval <= 0.01
  message(paste("Number of proteins tested:", length(unique(res$Protein))))
  message(paste("Number of peptides tested:", nrow(res)))
  message(paste("Number of uncorrelated peptides (1% FDR):",
                sum(signif_peptide)))
  message(paste("Number of proteins with uncorrelated peptides (1% FDR):",
                length(unique(res$Protein[signif_peptide]))))

  return(res)
}


#' @title The brains behind pecora_mod
#'
#' @description Tests for the interaction between the sample condition of
#'   interest and an indicator variable for a given peptide of interest. Loops
#'   over all unique peptides for a given protein.
#'
#' @param x a single-protein \code{data.frame} with columns "Peptide", "Sample",
#'   and "LogRatio".
#' @param median_mod logical; whether to use the median version of PeCorA. See
#'   \code{?pecora_mod} for details.
#'
#' @return a \code{data.frame} with columns "Peptide" and "pvalue".
#'
#' @importFrom purrr map_dbl
#' @importFrom dplyr %>% mutate group_by ungroup select distinct slice pull
#' @importFrom stats setNames lm
#' @importFrom car Anova
#' @importFrom tibble enframe
#' @importFrom tidyr unnest
#'
#' @noRd

pecora_test_single_protein <- function(x,
                                    median_mod = FALSE)
{
  unique_peptides <- unique(x$Peptide)

  # Single p-value per peptide
  res <- map_dbl(unique_peptides, function(pep_i) {
    curr_df <- mutate(x, current_peptide = Peptide == pep_i)

    if (median_mod) {
      curr_df <- curr_df %>%
        group_by(Sample, current_peptide) %>%
        mutate(LogRatio = median(LogRatio, na.rm = TRUE)) %>%
        ungroup() %>%
        select(Sample, LogRatio, Condition, current_peptide) %>%
        distinct()
    }

    # Slowest part
    pval <- lm(LogRatio ~ Condition * current_peptide,
               data = curr_df) %>%
      car::Anova() %>%
      slice(3) %>% # interaction term
      pull(`Pr(>F)`)

    return(pval)
  }) %>%
    setNames(nm = unique_peptides) %>%
    enframe(name = "Peptide",
            value = "pvalue")

  return(res)
}


utils::globalVariables(
  c("Protein", "Peptide", "Sample", "current_peptide", "LogRatio",
    "Condition", "Pr(>F)", "pvalue")
)

