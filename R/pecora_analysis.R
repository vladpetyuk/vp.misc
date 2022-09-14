#' Performs the PeCorA analysis described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8592057/.
#' Take a look at the linked paper to get an idea of the procedure.
#'
#' Takes an msnset, a list of proteins, and a treatment variable. msnset is assumed to be normalized, for instance it can be the output of
#' pecora_preprocessing, which normalizes as in the original PeCorA paper.
#' The PeCorA function is then used for the analysis. Returns a table containing Peptide, Protein, and significance.
#'
#' @param m msnset object containing a pData AND fData table. fData must have a Protein and Peptide column.
#' @param treatment_string the name of the column in pData containing treatment information. The data can be numerical of categorical.
#' @param proteins character vector indicating which proteins to analyze.
#' @param median_mod logical for whether to use the original PeCorA (median_mod = FALSE), or the median version (median_mod = TRUE). See the help of PeCorA_mod for more details.
#'
#' @return PeCorA results table for the supplied proteins. Saves plots of the significant peptides to the given folder.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom grid grobTree textGrob gpar
#'
#' @export pecora_analysis
#'
#'
pecora_analysis <- function(m, treatment_string, proteins, median_mod = FALSE){
  mat <- exprs(m)
  pval_label <- paste0("pecora_pval_", treatment_string)
  padj_label <- paste0("pecora_adj_pval_", treatment_string)

  ## Remove previous results using the same treatment variable.
  f_data <- fData(m) %>%
    select(setdiff(colnames(fData(m)), c(pval_label, padj_label))) %>%
    mutate(Protein = unlist(as.character(Protein)),
           feature_name = rownames(.))

  ## Must have a Protein column
  if (!("Protein" %in% colnames(f_data))){
    stop("fData should have 'Protein' column.\n")
  }

  ## Filter to just those proteins for which we want to run PeCorA
  peptide_mapping <- f_data %>%
    filter(Protein %in% proteins) %>%
    select(feature_name, Protein)

  metadata <- pData(m) %>%
    select(sym(treatment_string))
  metadata$Sample <- rownames(metadata)
  mat <- mat[peptide_mapping$feature_name, ]

  if (typeof(metadata[[1]]) == "character") {
    metadata[[1]] <- as.factor(metadata[[1]])
  }

  ## Prepare input for PeCorA function.
  PeCorA_input <- mat %>%
    as.data.frame() %>%
    mutate(feature_name = rownames(.)) %>%
    pivot_longer(cols = -feature_name, names_to = "Sample", values_to = "LogRatio") %>%
    mutate(modpep_z = feature_name,
           ms1adj = LogRatio) %>%
    merge(peptide_mapping, by = "feature_name") %>%
    merge(metadata, by = "Sample") %>%
    dplyr::rename(Condition = sym(treatment_string))

  PeCorA_result <- PeCorA_mod(PeCorA_input, median_mod) %>%
    dplyr::rename(feature_name = peptide,
                  Protein = protein)

  pval_index <- which(colnames(PeCorA_result) == "pvalue")[[1]]
  adj_pval_index <- which(colnames(PeCorA_result) == "adj_pval")[[1]]
  colnames(PeCorA_result)[[pval_index]] <- pval_label
  colnames(PeCorA_result)[[adj_pval_index]] <- padj_label

  ## Adding PeCorA results to fData.
  m_out <- m
  f_data <- f_data %>%
    mutate(feature_name = rownames(.)) %>%
    left_join(PeCorA_result %>% select(-Protein), by = "feature_name")
  rownames(f_data) <- f_data$feature_name

  fData(m_out) <- f_data[rownames(exprs(m_out)), ]
  return(m_out)
}


#' Plots a particular peptide using the PeCorA input and output tables.
#'
#' @param m msnset object containing a pData AND fData table. fData must have a Protein and Peptide column.
#' @param pecora_results PeCorA results table from pecora_analysis.
#' @param chosen_protein The protein of interest.
#' @param chosen_peptide The peptide of interest. Should map to the supplied protein.
#'
#' @return Plot of the specified peptide + protein.
#'
#' @import ggplot2
#' @importFrom grid grobTree textGrob gpar
#'
#' @export pecora_plot
#'
#'
pecora_plot <- function(m, chosen_feature, treatment_string = NULL,
                        median_mod = FALSE) {

  ## If empty treatment string, look for results in f_data.
  if (is.null(treatment_string)){
    index <- which(grepl("pecora_adj_pval", colnames(fData(m))))
    if (length(index) > 1){
      stop("Found multiple treatment variables in the PeCorA results.\n
           Must specify treatment variable.")
    } else if (length(index) == 0){
      stop("Did not find any PeCorA results")
    }
    treatment_string <- sub("^pecora_adj_pval_", "", colnames(fData(m))[[index]])
  }

  result_col <- paste0("pecora_adj_pval_", treatment_string)

  p.value <- try(fData(m) %>%
                   filter(feature_name == chosen_feature) %>%
                   pull(!!result_col) %>% round(12))
  if (inherits(p.value, "try-error")){
    stop(paste("Did not find PeCorA results for", treatment_string))
  }

  chosen_protein <- fData(m) %>% filter(feature_name == chosen_feature) %>%
    pull(Protein) %>% unlist() %>% unique()
  features <- fData(m) %>%
    filter(Protein == chosen_protein)

  if (length(p.value) == 0){
    message <- paste0(chosen_feature, " not found within the results.\n")
    stop(message)
  }

  metadata <- pData(m) %>%
    select(sym(treatment_string)) %>%
    mutate(Sample = rownames(.))

  if (typeof(metadata[[1]]) == "character") {
    metadata[[1]] <- as.factor(metadata[[1]])
  }

  plot_df <- exprs(m)[features$feature_name, ] %>%
    as.data.frame() %>%
    mutate(feature_name = rownames(.)) %>%
    pivot_longer(cols = -feature_name, names_to = "Sample", values_to = "value") %>%
    merge(features, by = "feature_name") %>%
    merge(metadata, by = "Sample") %>%
    dplyr::rename(Condition = sym(treatment_string)) %>%
    mutate(peptide_group = case_when(feature_name == chosen_feature ~ chosen_feature,
                                     TRUE ~ "All other peptides")) %>%
    mutate(peptide_group = factor(peptide_group,
                                  levels = c("All other peptides", chosen_feature))) %>%
    select(Sample, Condition, value, peptide_group)

  if (median_mod){
    plot_df <- plot_df %>%
      group_by(Sample, peptide_group) %>%
      mutate(value = median(value, na.rm = T)) %>%
      unique()
  }

  ## Boxplot if condition is a factor. Otherwise, scatterplots when using numerical data.
  if (is.factor(plot_df$Condition)) {

    grob <- grobTree(textGrob(paste0("Adj_pval = ", p.value),
                              x = 0.1,  y = 0.97, hjust = 0,
                              gp = gpar(col = "black", fontsize = 11)))

    p <- ggplot(plot_df, aes(x = Condition, y = value, fill = peptide_group)) +
      geom_boxplot(outlier.shape = NA) +
      ggtitle(chosen_feature) + annotation_custom(grob) +
      xlab(treatment_string) + ylab("Log Intensity") +
      theme(plot.title = element_text(hjust = 0.5))

  } else {

    plot_df <- plot_df %>%
      mutate(alpha = case_when(peptide_group == "All other peptides" ~ 0.005,
                               TRUE ~ 0.25))

    lm_df <- plot_df %>%
      filter(peptide_group == "All other peptides")
    allothers_lm <- lm(lm_df$value ~ lm_df$Condition)

    lm_df <- plot_df %>%
      filter(peptide_group == chosen_feature)
    chosen_lm <- lm(lm_df$value ~ lm_df$Condition)

    grob <- grobTree(textGrob(paste0("Adj_pval = ", p.value),
                              x = 0.1,  y = 0.97, hjust = 0,
                              gp = gpar(col = "black", fontsize = 11)))

    p <- ggplot(plot_df, aes(x = Condition, y = value, color = peptide_group, alpha = alpha)) + geom_point() +
      ggtitle(chosen_feature) + annotation_custom(grob) + xlab(treatment_string) +
      ylab("Log Intensity") + guides(alpha = "none") +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_abline(intercept = allothers_lm$coefficients[[1]],
                  slope = allothers_lm$coefficients[[2]],
                  color = "red", size = 1) +
      geom_abline(intercept = chosen_lm$coefficients[[1]],
                  slope = chosen_lm$coefficients[[2]],
                  color = scales::hue_pal()(2)[[2]], size = 1)
  }

  return(p)
}

#' This is more of a convenience function, used to normalize as was done in the original PeCorA paper.
#' Pass an msnset here to preprocess the data. The result is again an msnset. Can perform two normalization
#' steps, also carried out in the PeCorA paper. Use the parameters peptide_standardize and sample_standardize
#' to select which steps should be performed.
#'
#' If it's desired that both normalization steps be skipped (both set to FALSE), then pecora_preprocess does nothing
#' at all to the data, and this step can be safely skipped.
#'
#' @param m msnset object containing a pData AND fData table. fData must have a Protein and Peptide column.
#' @param treatment_string the name of the column in pData containing treatment information. The data can be numerical of categorical.
#' @param sample_standardize logical for whether exprs data should be standardized so that sample wise we have mean zero and variance one.
#' @param peptide_standardize logical for whether exprs data should be normalized so that logratio data is relative to the control group. When using numerical data, this determines whether to zero center each peptide.
#'
#' @return PeCorA results table for the supplied proteins. Saves plots of the significant peptides to the given folder.
#'
#' @export pecora_preprocess
#'
#'
pecora_preprocess <- function(m, treatment_string, control_group = NULL,
                              sample_standardize = FALSE, peptide_standardize = FALSE){

  mat <- exprs(m)
  metadata <- pData(m) %>%
    select(sym(treatment_string))
  metadata$Sample <- rownames(metadata)

  if (typeof(metadata[[1]]) == "character") {
    metadata[[1]] <- as.factor(metadata[[1]])
  }

  if (sample_standardize) {
    cat("Standardizing by sample.\n")
    sample_means <- apply(mat, 2, mean, na.rm = T)
    sample_sds <- apply(mat, 2, sd, na.rm = T)
    mat <- sweep(mat, 2, sample_means, FUN = '-')
    mat <- sweep(mat, 2, sample_sds, FUN = '/')
  }

  ## If peptide_standardize is true, then use the supplied control group to normalize the data
  if (peptide_standardize & is.factor(metadata[[1]])) {
    cat("Standardizing to control group.\n")
    if (is.null(control_group)) {
      stop("Must define 'control_group' in order to normalize relative to control\n")
    }

    samples_control <- metadata %>%
      filter(!!sym(treatment_string) == control_group) %>%
      rownames()

    if (length(samples_control) == 0){
      stop("No samples in control group. Please check value of 'control_group' and pData\n")
    }

    ## In case samples_control consists of a single sample we need as.matrix() in the apply call.
    peptide_means <- apply(as.matrix(mat[, samples_control]), 1, mean, na.rm = T)
    mat <- sweep(mat, 1, peptide_means, FUN = '-')

  } else if (peptide_standardize) {
    cat("Mean centering each peptide.\n")
    peptide_means <- apply(mat, 1, mean, na.rm = T)
    mat <- sweep(mat, 1, peptide_means, FUN = '-')
  }

  m.processed <- MSnSet(exprs = mat, fData = fData(m), pData = pData(m))
  return(m.processed)
}

