#' @title Plot PeCorA Results
#'
#' @param m \code{MSnSet} object; output of
#'   \code{\link[MSnSet.utils]{pecora_analysis}}.
#' @param feature character; the feature to plot. One of \code{featureNames(m)}.
#' @param treatment character; pattern used to search for columns of the form
#'   "pecora_adj_pval_{treatment}" in \code{fData(m)}. If \code{NULL} (default),
#'   this will automatically be determined if there is only one column that
#'   starts with "pecora_adj_pval_".
#' @inheritParams pecora_analysis
#'
#' @details This function was inspired by the boxplots presented in the PeCorA
#'   publication (\url{https://doi.org/10.1021/acs.jproteome.0c00602}).
#'
#' @seealso \code{\link[MSnSet.utils]{pecora_preprocess}},
#'   \code{\link[MSnSet.utils]{pecora_analysis}},
#'   \code{\link[MSnSet.utils]{pecora_mod}}
#'
#' @references Dermit, M., Peters-Clarke, T. M., Shishkova, E., & Meyer, J. G.
#'   (2021). Peptide Correlation Analysis (PeCorA) Reveals Differential
#'   Proteoform Regulation. \emph{Journal of proteome research, 20(4)},
#'   1972–1980. \url{https://doi.org/10.1021/acs.jproteome.0c00602}
#'
#' @return A \code{ggplot} object.
#'
#' @import ggplot2
#' @importFrom scales hue_pal
#' @importFrom rlang !! sym
#' @importFrom dplyr %>% filter pull select mutate rename group_by distinct
#' @importFrom tidyr pivot_longer
#' @importFrom stats lm update coef
#'
#' @export plot_pecora

plot_pecora <- function(m,
                        feature,
                        treatment = NULL,
                        median_mod = FALSE)
{
  f_data <- fData(m) %>%
    mutate(feature_name = rownames(.))

  if (!feature %in% rownames(m)) {
    stop(sprintf("'%s' is not in featureNames(m)", feature))
  }

  ## If empty treatment string, look for results in f_data.
  if (is.null(treatment)) {
    index <- grep("pecora_adj_pval", colnames(f_data))

    if (length(index) > 1) {
      stop(paste0("Found multiple treatment variables in the PeCorA results.",
                  "Must specify treatment variable."))
    } else if (length(index) == 0) {
      stop("Did not find any PeCorA results")
    }

    treatment <- sub("^pecora_adj_pval_", "", colnames(f_data)[[index]])
  }

  result_col <- paste0("pecora_adj_pval_", treatment)

  padj <- try(
    f_data %>%
      filter(feature_name == feature) %>%
      pull(result_col) %>%
      signif(3)
  )

  if (inherits(padj, "try-error")) {
    stop(paste("Did not find PeCorA results for", treatment))
  }

  if (is.na(padj)) {
    stop(sprintf("%s missing PeCorA p-value when treatment = %s",
                 feature, treatment))
  }

  chosen_protein <- f_data[feature, ] %>%
    pull(Protein) %>%
    unique()

  features <- filter(f_data, Protein == chosen_protein)

  if (length(padj) == 0) {
    stop(paste(feature, "not found within the results."))
  }

  metadata <- pData(m) %>%
    select(sym(treatment)) %>%
    mutate(Sample = rownames(.))

  if (is.character(metadata[[treatment]])) {
    metadata[[treatment]] <- as.factor(metadata[[treatment]])
  }

  plot_df <- exprs(m)[features$feature_name, ] %>%
    as.data.frame() %>%
    mutate(feature_name = rownames(.)) %>%
    pivot_longer(cols = -feature_name,
                 names_to = "Sample",
                 values_to = "value") %>%
    merge(features, by = "feature_name") %>%
    merge(metadata, by = "Sample") %>%
    dplyr::rename(Condition = sym(treatment)) %>%
    mutate(peptide_group = ifelse(feature_name == feature,
                                  feature, "All other peptides")) %>%
    mutate(peptide_group = factor(peptide_group,
                                  levels = c("All other peptides",
                                             feature))) %>%
    select(Sample, Condition, value, peptide_group)

  if (median_mod) {
    plot_df <- plot_df %>%
      group_by(Sample, peptide_group) %>%
      mutate(value = median(value, na.rm = TRUE)) %>%
      distinct()
  }

  ## Boxplot if condition is a factor. Otherwise, scatterplots when using
  ## numerical data.
  if (is.factor(plot_df$Condition)) {
    p <- ggplot(plot_df, aes(x = Condition, y = value,
                             fill = peptide_group)) +
      geom_boxplot(outlier.shape = NA, na.rm = TRUE)
  } else {
    plot_df <- plot_df %>%
      mutate(alpha = ifelse(peptide_group == "All other peptides",
                            0.005, 0.25))

    lm_df <- filter(plot_df, peptide_group == "All other peptides")
    allothers_lm <- lm(value ~ Condition, data = lm_df)

    lm_df <- filter(plot_df, peptide_group == feature)
    chosen_lm <- update(allothers_lm, data = lm_df)

    p <- ggplot(plot_df, aes(x = Condition, y = value,
                             color = peptide_group, alpha = alpha)) +
      geom_point(na.rm = TRUE) +
      geom_abline(intercept = coef(allothers_lm)[1],
                  slope = coef(allothers_lm)[2],
                  color = "red", size = 1) +
      geom_abline(intercept = coef(chosen_lm)[1],
                  slope = coef(chosen_lm)[2],
                  color = scales::hue_pal()(2)[2], size = 1) +
      guides(alpha = "none")
  }

  p <- p +
    labs(title = feature,
         subtitle = paste("BH-adjusted p-value =", padj),
         x = treatment,
         y = "Log Intensity")

  return(p)
}


utils::globalVariables(
  c("Protein", "Peptide", "pvalue", "peptide_group", "Sample", "Condition")
)

