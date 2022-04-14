#' Test Custom Contrasts with LIMMA
#'
#' A convenience wrapper for differential analysis with LIMMA using custom
#' contrasts.
#'
#' @param eset \code{eSet} (or most likely \code{eSet} subclass) object
#' @param model.str character; formulation of the model (e.g. \code{"~ 0 + a +
#'   b"}). This should be a no-intercept model (it should include 0 or -1 as
#'   terms). See \code{\link[stats]{lm}} for more details.
#' @param coef.str character; coefficient of interest. One of
#'   \code{colnames(pData(eset))}. E.g. \code{"a"}.
#' @param contrasts character; contrasts to test. All factor levels in the
#'   contrasts should begin with the name of that factor. It is recommended to
#'   use \code{\link{paircomp}} with the \code{name} argument specified to
#'   create pairwise contrasts of the form \code{"a-b"}.
#' @param stat character; the moderated test statistic to use. Either \code{"t"}
#'   or \code{"F"}.
#' @param var.group \code{NULL} or character; the column in \code{pData(eset)}
#'   indicating groups that will have different relative quality weights. The
#'   default (\code{NULL}), means that each sample will be weighted equally.
#' @param plot logical; whether to generate diagnostic plots. The default
#'   (\code{TRUE}) generates a barplot of weights applied to each sample and a
#'   plot of the mean-variance trend using \code{\link[limma]{plotSA}}.
#' @param col.group character; the column in \code{pData(eset)} used to color
#'   the text of the MDS plot.
#' @param adjust.method method for p-value adjustment. Default is \code{"BH"}
#'   (Benjamini-Hochberg). See \code{\link[stats]{p.adjust}} for details.
#' @param remove_name if \code{TRUE} (the default) \code{coef.str} will be
#'   removed from the contrasts. This simply improves readability, but it will
#'   not work properly if \code{coef.str} is a substring of any group name.
#' @param ... additional arguments passed to \code{\link[limma]{eBayes}}.
#'
#' @return \code{data.frame}. Basically the output of
#'   \code{\link[limma]{topTable}} with additional columns \code{feature} and
#'   \code{contrast} (if \code{stat = "t"}). All columns from \code{fData} are
#'   also included.
#'
#' @details An MDS plot is used to determine the appropriate value of
#'   \code{var.group}. If samples within phenotype groups tend to cluster well
#'   and have similar dispersions, the default \code{var.group = NULL} is
#'   recommended. If samples within phenotype groups tend to cluster well, but
#'   different groups are more or less dispersed, it may be beneficial to set
#'   \code{var.group} to the name of the phenotype column. If one or more
#'   samples tend to cluster poorly with samples of the same phenotype, it may
#'   be beneficial to weight by sample. That is, set \code{var.group} to the
#'   name of the column in \code{pData(eset)} that uniquely identifies each
#'   sample. If variation between samples or groups of samples is biological
#'   rather than technical, weighting is not recommended.
#'
#'   The plot of the mean-variance trend helps determine whether to fit a trend
#'   line to the prior variance (default is \code{trend = FALSE}). It also helps
#'   determine if the prior variance should be robustified against outlier
#'   sample variances (\code{robust = TRUE}; default is \code{FALSE}). See
#'   \code{\link[limma]{eBayes}} for more details.
#'
#'   p-values are adjusted across all contrasts (globally). It is recommended to
#'   adjust p-values from similar contrasts together. See the "Multiple Testing
#'   Across Contrasts" section of the LIMMA User's Guide (linked below) for more
#'   information.
#'
#' @seealso
#' \href{https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf}{LIMMA
#' User's Guide}, \code{\link{limma_a_b}}, \code{\link{limma_gen}}
#'
#' @references Liu, R., Holik, A. Z., Su, S., Jansz, N., Chen, K., Leong, H. S.,
#'   Blewitt, M. E., Asselin-Labat, M. L., Smyth, G. K., & Ritchie, M. E.
#'   (2015). Why weight? Modelling sample and observational level variability
#'   improves power in RNA-seq analyses. *Nucleic acids research*, *43*(15),
#'   e97. \url{https://doi.org/10.1093/nar/gkv412}
#'
#'   Phipson, B., Lee, S., Majewski, I. J., Alexander, W. S., & Smyth, G. K.
#'   (2016). ROBUST HYPERPARAMETER ESTIMATION PROTECTS AGAINST HYPERVARIABLE
#'   GENES AND IMPROVES POWER TO DETECT DIFFERENTIAL EXPRESSION. *The annals of
#'   applied statistics*, *10*(2), 946--963.
#'   \url{https://doi.org/10.1214/16-AOAS920}
#'
#'   Smyth G. K. (2004). Linear models and empirical bayes methods for assessing
#'   differential expression in microarray experiments. *Statistical
#'   applications in genetics and molecular biology*, *3*, Article3.
#'   \url{https://doi.org/10.2202/1544-6115.1027}
#'
#' @md
#'
#' @importFrom MSnbase exprs pData
#' @importFrom stats terms model.matrix p.adjust
#' @importFrom data.table rbindlist
#' @importFrom graphics barplot
#' @import limma
#' @import statmod
#'
#' @export limma_contrasts
#'
#' @examples
#' # library(MSnSet.utils)
#' data(cptac_oca)
#'
#' # A no-intercept model is required. Testing SUBTYPE contrasts
#' # with AGE as covariate
#' model.str <- "~ 0 + SUBTYPE + AGE"
#' coef.str <- "SUBTYPE"
#'
#' # Create contrasts
#' contrasts <- paircomp(pData(oca.set)[, coef.str], name = coef.str)
#'
#' # Moderated t-statistics
#' res1 <- limma_contrasts(oca.set, model.str = model.str,
#'                         coef.str = coef.str,
#'                         contrasts = contrasts)
#'
#' # Moderated F-statistics
#' res2 <- limma_contrasts(oca.set, model.str = model.str,
#'                         coef.str = coef.str,
#'                         contrasts = contrasts, stat = "F")
#'
#' # Count number of significant (< 0.05) features
#' sum(res2$adj.P.Val < 0.05)
#'
#' # Robustify mean-variance trend (based on MDS plot)
#' res3 <- limma_contrasts(oca.set, model.str = model.str,
#'                         coef.str = coef.str,
#'                         contrasts = contrasts,
#'                         trend = TRUE, robust = TRUE)
#'
#' # Count number of significant (< 0.05) features in each contrast
#' table(res3$contrast, res3$adj.P.Val < 0.05)
#'


limma_contrasts <- function(eset, model.str, coef.str, contrasts,
                            stat = c("t", "F"), var.group = character(0),
                            plot = FALSE, col.group = coef.str,
                            adjust.method = "BH", remove_name = TRUE,
                            ...) {

    # test statistic
    stat <- match.arg(arg = stat, choices = c("t", "F"))

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
    if (!identical(var.group, character(0))) {
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
    if (plot) limma_diagnostics(eset, col.group, coef.str, weights, fit.smooth)

    # topTable for each contrast or for all contrasts
    if (stat == "t") { coef <- contrasts } else { coef <- list(contrasts) }
    res <- lapply(
        coef,
        function(contrast_i) {
            x <- topTable(fit.smooth, coef = contrast_i, number = nrow(eset))
            x$feature <- rownames(x) # add feature column
            # Add contrast column
            if (stat == "t") {
                x$contrast <- contrast_i
            }
            x
        })
    res <- rbindlist(res)

    # Adjust p-values across all contrasts (global adjustment)
    res$adj.P.Val <- p.adjust(res$P.Value, method = adjust.method)

    # Remove coef.str from the contrasts
    if (remove_name & stat == "t") {
        res$contrast <- gsub(coef.str, "", res$contrast)
    }

    return(res)
}


## Helper functions ----

# diagnostic plots
limma_diagnostics <- function(eset, col.group, coef.str, weights, fit.smooth) {
    # MDS plot
    col.group <- as.factor(pData(eset)[, col.group])
    levels(col.group) <- jet2.colors(nlevels(col.group))
    col.group <- as.character(col.group)
    plotMDS(exprs(eset), col = col.group,
            labels = pData(eset)[, coef.str, drop = TRUE],
            main = "Multidimensional Scaling Plot")

    # Barplot of sample weights
    barplot(weights, space = 0,
            xlab = "Sample", ylab = "Weight",
            main = "Sample-specific weights")
    abline(1, 0, lty = "longdash", col = "red")

    # Plot of the mean-variance trend
    plotSA(fit.smooth, cex = 1, main = "Mean-variance trend")
}



