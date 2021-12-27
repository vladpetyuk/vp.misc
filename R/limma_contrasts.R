#' Wrapper for LIMMA
#'
#' A convenience wrapper for testing differential expression/abundance with
#' custom comparisons (contrasts) using the LIMMA package.
#'
#' @param eset eset (or most likely eset subclass) object
#' @param model.str character formulation of the model (e.g. "~ 0 + a + b").
#'     This should be a no-intercept model (it should include 0 or -1 as
#'     terms). See \code{\link[stats]{lm}} for more details.
#' @param coef.str vector of character coefficients of interest from
#'      \code{colnames(pData(eset))}. E.g. "a" or c("a", "b").
#' @param contrasts character vector of contrasts to test. All factor levels in
#'     the contrasts should begin with the name of that factor. It is
#'     recommended to use \code{\link[MSnSet.utils]{pairwise_contrasts}} to
#'     create the contrasts.
#' @param var.group NULL or a string; the column in \code{pData(eset)}
#'     indicating groups that will have different relative quality weights.
#'     The default (\code{NULL}), means that each sample will be weighted
#'     equally.
#' @param plot logical; whether to generate diagnostic plots. The default
#'     (\code{TRUE}) generates a PCA plot using
#'     \code{\link[MSnSet.utils]{plot_pca}}, a barplot of weights applied to
#'     each sample, and a plot of the mean-variance trend.
#' @param label.by \code{NULL} or a string; the column in \code{pData(eset)}
#'     used to label the samples in the PCA plot. If \code{NULL} (default),
#'     \code{sampleNames(eset)} will be used.
#' @param color.by \code{NULL} or a string; the column in \code{pData(eset)}
#'     used to determine how to color samples in the PCA plot. If \code{NULL}
#'     (default), each sample will be colored black.
#' @param adjust.method method for p-value adjustment. Default is \code{"BH"}
#'     (Benjamini-Hochberg).
#' @param remove_name if \code{TRUE} (the default) the name of \code{coef.str}
#'     will be removed from the contrasts. This simply improves readability.
#' @param ... additional arguments passed to \code{\link[limma]{eBayes}}.
#'
#' @return data.frame. Basically output of \code{\link[limma]{topTable}}
#'      function with additional columns \code{feature} and \code{contrast}.
#'
#' @details
#' The PCA plot is used to determine the appropriate value of \code{var.group}.
#' If samples within phenotype groups tend to cluster well and have similar
#' dispersions, the default \code{var.group = NULL} is recommended. If samples
#' within phenotype groups tend to cluster well, but different groups are more
#' or less dispersed, it is recommended to set \code{var.group} to the name of
#' the phenotype column. If one or more samples tend to cluster poorly with
#' samples of the same phenotype, it is recommended to weight by sample. That
#' is, set \code{var.group} to the name of the column in \code{pData(eset)} that
#' uniquely identifies each sample.
#'
#' The plot of the mean-variance trend helps determine whether to fit a trend
#' line to the prior variance (default is \code{trend = FALSE}). It also helps
#' determine if the prior variance should be robustified against outlier sample
#' variances (\code{robust = TRUE}; default is \code{FALSE}). See
#' \code{\link[limma]{eBayes}} for more details.
#'
#' p-values are adjusted across all contrasts (globally). It is recommended to
#' adjust p-values from similar contrasts together. See the limma User's Guide
#' linked below for more information.
#'
#' @seealso
#' \url{https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf}
#'
#'
#' @importFrom MSnbase exprs pData
#' @importFrom dplyr mutate select %>% everything
#' @importFrom stats terms model.matrix p.adjust
#' @importFrom data.table rbindlist
#' @importFrom graphics barplot
#' @import limma
#' @import statmod
#'
#' @export limma_contrasts
#'
#' @examples
#' library(MSnSet.utils)
#' data(cptac_oca)
#'
#' # A no-intercept model is required
#' model.str <- "~ 0 + SUBTYPE + AGE"
#' coef.str <- "SUBTYPE"
#'
#' # Create contrasts
#' contrasts <- pairwise_contrasts(pData(oca.set), fct = coef.str)
#'
#' res <- limma_contrasts(oca.set,
#'                        model.str = model.str,
#'                        coef.str = coef.str,
#'                        contrasts = contrasts)
#'
#' # There is a clear mean-variance trend, so we need to capture that
#' res <- limma_contrasts(oca.set,
#'                        model.str = model.str,
#'                        coef.str = coef.str,
#'                        contrasts = contrasts,
#'                        trend = TRUE)
#'
#' # We need to robustify this trend against the presence of outliers
#' res <- limma_contrasts(oca.set,
#'                        model.str = model.str,
#'                        coef.str = coef.str,
#'                        contrasts = contrasts,
#'                        trend = TRUE,
#'                        robust = TRUE)
#'
#' # The mean-variance trend plot looks good, but we need to check the PCA plot.
#' # We will label by the "Batch" column to better see how samples cluster.
#' res <- limma_contrasts(oca.set,
#'                        model.str = model.str,
#'                        coef.str = coef.str,
#'                        contrasts = contrasts,
#'                        trend = TRUE,
#'                        robust = TRUE,
#'                        color.by = "SUBTYPE",
#'                        label.by = "Batch")
#'
#' # It looks like groups have different dispersions and one or more samples
#' # within a group do not agree with the others, so we will weight by sample.
#' # We need to first create a column for the sample names.
#' oca.set$sample_name <- sampleNames(oca.set)
#'
#' res <- limma_contrasts(oca.set,
#'                        model.str = model.str,
#'                        coef.str = coef.str,
#'                        contrasts = contrasts,
#'                        trend = TRUE,
#'                        robust = TRUE,
#'                        color.by = "SUBTYPE",
#'                        label.by = "Batch",
#'                        var.group = "sample_name")
#'
#' # Count number of significant features in each contrast
#' library(dplyr)
#' # Using a significance cutoff of 0.05
#' filter(res, adj.P.Val < 0.05) %>%
#'   group_by(contrast) %>%
#'   tally()
#'
#' head(res)


limma_contrasts <- function(eset, model.str, coef.str, contrasts,
                            var.group = NULL, plot = TRUE,
                            label.by = NULL, color.by = NULL,
                            adjust.method = "BH", remove_name = TRUE,
                            ...) {

    model.formula <- eval(parse(text = model.str), envir = pData(eset))

    if (attr(terms(model.formula), which = "intercept") != 0) {
        stop(paste("Please specify a no-intercept model for model.str",
                   "See ?lm for more details.", sep = "\n"))
    }

    design <- model.matrix(model.formula)

    # Some samples may be NA, so they will be dropped in the design matrix.
    # Subset eset to only those rows present in the design matrix.
    eset <- eset[, as.numeric(rownames(design))]

    contrast.matrix <- makeContrasts(contrasts = contrasts,
                                     levels = design)

    # If var.group is specified, determine group or sample-level
    # quality weights. Otherwise, set all weights to 1.
    if (!is.null(var.group)) {
        weights <- arrayWeights(eset, design,
                                var.group = pData(eset)[, var.group])
    } else {
        weights <- rep(1, ncol(eset))
        names(weights) <- seq_along(weights)
    }

    # Linear modeling
    fit <- lmFit(eset, design, weights = weights)
    fit.contr <- contrasts.fit(fit, contrast.matrix)
    fit.smooth <- eBayes(fit.contr, ...)

    # Diagnostic plots
    if (plot) {

        # Sample PCA plot
        plot_pca(eset, phenotype = color.by, label = label.by)

        # Barplot of sample weights
        barplot(weights, space = 0,
                xlab = "Sample", ylab = "Weight",
                main = "Sample-specific weights")
        abline(1, 0, lty = "longdash", col = "red")

        # Plot of the mean-variance trend
        plotSA(fit.smooth, cex = 1,
               main = "Mean-variance trend")
    }


    res <- lapply(contrasts,
                  function(contrast_i) {
                      topTable(fit.smooth,
                               coef = contrast_i,
                               number = nrow(eset)) %>%
                          mutate(feature = rownames(.),
                                 contrast = contrast_i)
                  }) %>%
        rbindlist() %>%
        # Adjust p-values across all contrasts (global adjustment)
        mutate(adj.P.Val = p.adjust(P.Value, method = adjust.method)) %>%
        # Reorder columns
        select(feature, contrast, everything())

    # Remove coef.str from the contrasts
    if (remove_name) {
        res$contrast <- gsub(coef.str, "", res$contrast)
    }

    return(res)
}

utils::globalVariables(c(".", "feature", "contrast", "P.Value"))

