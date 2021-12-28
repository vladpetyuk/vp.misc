#' Test Custom Contrasts with LIMMA
#'
#' A convenience wrapper for differential analysis with LIMMA using
#' custom contrasts.
#'
#' @param eset \code{eSet} (or most likely \code{eSet} subclass) object
#' @param model.str character; formulation of the model (e.g.
#'     \code{"~ 0 + a + b"}).
#'     This should be a no-intercept model (it should include 0 or -1 as
#'     terms). See \code{\link[stats]{lm}} for more details.
#' @param coef.str character; coefficient of interest. One of
#'      \code{colnames(pData(eset))}. E.g. \code{"a"}.
#' @param contrasts character; contrasts to test. All factor levels in
#'     the contrasts should begin with the name of that factor. It is
#'     recommended to use \code{\link{paircomp}} with the
#'     \code{name} argument specified to create pairwise contrasts of the form
#'     \code{"a-b"}.
#' @param var.group \code{NULL} or character; the column in \code{pData(eset)}
#'     indicating groups that will have different relative quality weights.
#'     The default (\code{NULL}), means that each sample will be weighted
#'     equally.
#' @param plot logical; whether to generate diagnostic plots. The default
#'     (\code{TRUE}) generates a barplot of weights applied to
#'     each sample and a plot of the mean-variance trend using
#'     \code{\link[limma]{plotSA}}.
#' @param adjust.method method for p-value adjustment. Default is \code{"BH"}
#'     (Benjamini-Hochberg).
#' @param remove_name if \code{TRUE} (the default) \code{coef.str}
#'     will be removed from the contrasts. This simply improves readability, but
#'     it will not work properly if \code{coef.str} is a substring of any group
#'     name.
#' @param ... additional arguments passed to \code{\link[limma]{eBayes}}.
#'
#' @return
#' \code{data.frame}. Basically the output of \code{\link[limma]{topTable}}
#' with additional columns \code{feature} and \code{contrast}.
#'
#' @details
#' A PCA plot is used to determine the appropriate value of \code{var.group}.
#' If samples within phenotype groups tend to cluster well and have similar
#' dispersions, the default \code{var.group = NULL} is recommended. If samples
#' within phenotype groups tend to cluster well, but different groups are more
#' or less dispersed, it may be beneficial to set \code{var.group} to the name
#' of the phenotype column. If one or more samples tend to cluster poorly with
#' samples of the same phenotype, it may be beneficial to weight by sample. That
#' is, set \code{var.group} to the name of the column in \code{pData(eset)} that
#' uniquely identifies each sample. If variation between samples or groups of
#' samples is biological rather than technical, weighting is not recommended.
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
#' \href{https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf}{LIMMA User's Guide},
#' \code{\link{limma_a_b}},
#' \code{\link{limma_gen}},
#' \code{\link{plot_pca}}
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
#' contrasts <- paircomp(pData(oca.set)[, coef.str], name = coef.str)
#'
#' # Basic specification
#' res1 <- limma_contrasts(oca.set,
#'                         model.str = model.str,
#'                         coef.str = coef.str,
#'                         contrasts = contrasts)
#'
#' # Show diagnostic plots and robustify mean-variance trend.
#' res2 <- limma_contrasts(oca.set, model.str = model.str, coef.str = coef.str,
#'                         contrasts = contrasts, plot = TRUE,
#'                         trend = TRUE, robust = TRUE)
#'
#' # Count number of significant (< 0.05) features in each contrast
#' library(dplyr)
#' filter(res2, adj.P.Val < 0.05) %>%
#'   group_by(contrast) %>%
#'   tally()
#'


limma_contrasts <- function(eset, model.str, coef.str, contrasts,
                            var.group = NULL, plot = TRUE,
                            # label.by = NULL, color.by = NULL,
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

    ## Diagnostic plots
    if (plot) {

        # # Sample PCA plot
        # plot_pca(eset, phenotype = color.by, label = label.by)

        # Barplot of sample weights
        barplot(weights, space = 0,
                xlab = "Sample", ylab = "Weight",
                main = "Sample-specific weights")
        abline(1, 0, lty = "longdash", col = "red")

        # Plot of the mean-variance trend
        plotSA(fit.smooth, cex = 1,
               main = "Mean-variance trend")
    }


    # topTable for each contrast
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

