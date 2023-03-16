#' @title Preprocess Protein Expression Data for PeCorA
#'
#' @description This function optionally standardizes samples and/or
#'   zero-centers peptides using the overall mean of each peptide or the mean of
#'   the \code{control_group} samples.
#'
#' @inheritParams PeCorA
#'
#' @details If both \code{standardize_samples} and \code{center_peptides} are
#'   set to \code{FALSE}, the MSnSet will not be modified, and this step can be
#'   safely skipped.
#'
#'   This function was based on the \code{PeCorA_preprocessing} function from
#'   the PeCorA R package (\url{https://github.com/jessegmeyerlab/PeCorA}).
#'
#' @return An \code{MSnSet} object with modified exprs data.
#'
#' @seealso \code{\link[MSnSet.utils]{normalizeByGlob}},
#'   \code{\link[MSnSet.utils]{PeCorA}}
#'
#' @references Dermit, M., Peters-Clarke, T. M., Shishkova, E., & Meyer, J. G.
#'   (2021). Peptide Correlation Analysis (PeCorA) Reveals Differential
#'   Proteoform Regulation. \emph{Journal of proteome research, 20(4)},
#'   1972–1980. \url{https://doi.org/10.1021/acs.jproteome.0c00602}
#'
#' @importFrom dplyr %>% select filter
#' @importFrom rlang !! sym
#'
#' @noRd

.PeCorA_preprocess <- function(m,
                               treatment,
                               control_group = NULL,
                               standardize_samples = TRUE,
                               center_peptides = TRUE)
{
  mat <- exprs(m)
  metadata <- select(pData(m), !!sym(treatment))
  metadata$Sample <- rownames(metadata)

  if (is.character(metadata[[treatment]])) {
    metadata[[treatment]] <- as.factor(metadata[[treatment]])
  }

  if (standardize_samples) {
    mat <- scale(mat)
  }

  ## If center_peptides is true, then use the supplied control group to
  ## normalize the data
  if (center_peptides) {
    if (is.factor(metadata[[treatment]])) {
      if (is.null(control_group)) {
        stop(paste("control_group must be specified if treatment is",
                   "categorical and center_peptides=TRUE"))
      }
      control_samples <- metadata %>%
        filter(!!sym(treatment) == control_group) %>%
        rownames()

      if (length(control_samples) == 0) {
        stop("treatment column of pData does not contain control_group.")
      }

      peptide_means <- rowMeans(mat[, control_samples, drop = FALSE],
                                na.rm = TRUE)
    } else {
      peptide_means <- rowMeans(mat, na.rm = TRUE)
    }

    mat <- sweep(mat, 1, peptide_means, FUN = '-')
  }

  exprs(m) <- mat

  return(m)
}

# ==============================================================================

#' @title Peptide Correlation Analysis
#'
#' @description This is a modified version of the PeCorA function. For the
#'   original publication see
#'   \url{https://doi.org/10.1021/acs.jproteome.0c00602}.
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
#' @param x PeCorA input table. See the \code{\link[MSnSet.utils]{PeCorA}}
#'   function for details.
#' @inheritParams PeCorA
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
#' @importFrom tidyselect all_of
#' @importFrom rlang !! sym
#'
#' @noRd

.PeCorA_all_proteins <- function(x,
                                 protein_column = "Protein",
                                 peptide_column = "Peptide",
                                 median_mod = TRUE)
{
  proteins <- x %>%
    select(all_of(c(protein_column, peptide_column))) %>%
    distinct() %>%
    group_by(!!sym(protein_column)) %>%
    filter(n() > 1) %>% # at least 2 peptides per protein
    ungroup() %>%
    pull(!!sym(protein_column)) %>%
    unique()

  res <- x %>%
    filter(!!sym(protein_column) %in% proteins) %>%
    split.data.frame(f = .[[protein_column]], drop = TRUE) %>%
    map(.PeCorA_single_protein,
        peptide_column = peptide_column,
        median_mod = median_mod,
        .progress = TRUE) %>%
    enframe(name = protein_column) %>%
    unnest(value) %>%
    group_by(!!sym(protein_column)) %>%
    mutate(adj_pval = p.adjust(pvalue, method = "BH")) %>%
    ungroup()

  # Results summary
  signif_peptide <- res$adj_pval <= 0.01
  message(paste("Number of proteins tested:",
                length(unique(res[[protein_column]]))))
  message(paste("Number of peptides tested:", nrow(res)))
  message(paste("Number of uncorrelated peptides (1% FDR):",
                sum(signif_peptide)))
  message(paste("Number of proteins with uncorrelated peptides (1% FDR):",
                length(unique(res[signif_peptide, protein_column]))))

  return(res)
}

# ==============================================================================

#' @title PeCorA for a single protein
#'
#' @description Tests for the interaction between the sample condition of
#'   interest and an indicator variable for a given peptide of interest. Loops
#'   over all unique peptides for a given protein.
#'
#' @param x a single-protein \code{data.frame} with columns
#'   \code{peptide_column}, "Sample", and "LogRatio".
#' @inheritParams PeCorA
#'
#' @return a \code{data.frame} with columns "Peptide" and "pvalue".
#'
#' @importFrom purrr map_dbl
#' @importFrom dplyr %>% mutate group_by summarise slice pull
#' @importFrom stats setNames lm
#' @importFrom car Anova
#' @importFrom tibble enframe
#' @importFrom tidyr unnest
#' @importFrom rlang !! sym
#'
#' @noRd

.PeCorA_single_protein <- function(x,
                                   peptide_column = "Peptide",
                                   median_mod = TRUE)
{
  unique_peptides <- unique(x[[peptide_column]])

  # Single p-value per peptide
  res <- map_dbl(unique_peptides, function(pep_i) {
    curr_df <- mutate(x, current_peptide = !!sym(peptide_column) == pep_i)

    if (median_mod) {
      curr_df <- curr_df %>%
        group_by(Condition, Sample, current_peptide) %>%
        summarise(LogRatio = median(LogRatio, na.rm = TRUE))
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
    enframe(name = peptide_column,
            value = "pvalue")

  return(res)
}


utils::globalVariables(
  c("Sample", "current_peptide", "LogRatio", "Condition", "Pr(>F)", "pvalue")
)

