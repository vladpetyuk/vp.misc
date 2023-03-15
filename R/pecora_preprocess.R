#' @title Preprocess Protein Expression Data for PeCorA
#'
#' @description This function optionally standardizes samples and/or
#'   zero-centers peptides using the overall mean of each peptide or the mean of
#'   the \code{control_group} samples.
#'
#' @inheritParams pecora_analysis
#' @param control_group character; the control group, if \code{treatment} is
#'   categorical, this is one of the levels.
#' @param standardize_samples logical; whether to standardize samples (columns)
#'   in exprs data.
#' @param center_peptides logical; whether to mean-center each peptide. If the
#'   data in the \code{treatment} column is of type "character" or "factor",
#'   peptides will be zero-centered on the control group mean. Otherwise, the
#'   mean of all samples are used for centering.
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
#'   \code{\link[MSnSet.utils]{pecora_analysis}},
#'   \code{\link[MSnSet.utils]{pecora_mod}},
#'   \code{\link[MSnSet.utils]{plot_pecora}}
#'
#' @references Dermit, M., Peters-Clarke, T. M., Shishkova, E., & Meyer, J. G.
#'   (2021). Peptide Correlation Analysis (PeCorA) Reveals Differential
#'   Proteoform Regulation. \emph{Journal of proteome research, 20(4)},
#'   1972–1980. \url{https://doi.org/10.1021/acs.jproteome.0c00602}
#'
#' @importFrom dplyr select filter
#' @importFrom rlang !! sym
#'
#' @export pecora_preprocess

pecora_preprocess <- function(m,
                              treatment,
                              control_group = NULL,
                              standardize_samples = TRUE,
                              center_peptides = TRUE)
{
  mat <- exprs(m)
  metadata <- select(pData(m), sym(treatment))
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
        stop("control_group must be provided if center_peptides = TRUE")
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

