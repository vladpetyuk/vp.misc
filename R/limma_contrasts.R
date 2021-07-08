#' Wrapper for Limma
#'
#' A convenience wrapper for limma contrast testing.
#'
#' @param eset eset (or most likely eset subclass) object
#' @param model.str character formulation of the model (e.g. "~ 0 + a + b").
#'     This should be a no-intercept model (it should include 0 or -1 as
#'     terms).
#' @param coef.str vector of character coefficients of interest from
#'      \code{colnames(pData(eset))}. E.g. "a" or c("a", "b").
#' @param contrasts character vector of contrasts to test. All factor levels in
#'     the contrasts should begin with the name of that factor. It is
#'     recommended to use \code{[MSnSet.utils]{pairwise_contrasts}} to create
#'     the contrasts.
#' @param adjust.method method for p-value adjustment. Default is "BH"
#'     (Benjamini-Hochberg).
#' @param remove_name if \code{TRUE} (the default) the name of \code{coef.str}
#'     will be removed from the contrasts. This simply improves readability.
#' @param ... arguments for \code{\link[limma]{lmFit}}
#'
#' @return data.frame. Basically output of \code{\link[limma]{topTable}}
#'      function with additional columns \code{feature} and \code{contrast}.
#'
#' @importFrom MSnbase exprs pData
#' @importFrom limma lmFit topTable eBayes contrasts.fit makeContrasts
#' @importFrom dplyr mutate select %>%
#' @importFrom tidyr everything
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
#'                        covariates = covariates,
#'                        contrasts = contrasts)
#'
#' head(res)


limma_contrasts <- function(eset, model.str, coef.str, contrasts,
                            adjust.method = "BH", remove_name = TRUE,
                            ...) {

    model.formula <- eval(parse(text = model.str), envir = pData(eset))

    if (attr(terms(model.formula), which = "intercept") != 0) {
        stop(paste("Please specify a no-intercept model for model.str",
                   "See lm() documentation for more details.", sep = "\n"))
    }

    design <- model.matrix(model.formula)

    # Some samples may be NA, so they will be dropped in the design matrix.
    # Subset eset to only those rows present in the design matrix.
    eset <- eset[, as.numeric(rownames(design))]

    contrast.matrix <- makeContrasts(contrasts = contrasts,
                                     levels = design)

    fit.smooth <-
        lmFit(exprs(eset), design, ...) %>%
        contrasts.fit(contrast.matrix) %>%
        eBayes()

    res <- lapply(contrasts,
                  function(contrast_i) {
                      topTable(fit.smooth,
                               coef = contrast_i,
                               number = nrow(eset)) %>%
                          mutate(feature = rownames(.),
                                 contrast = contrast_i)
                  })

    # Bind topTable output
    res <- do.call(rbind, res)
    rownames(res) <- NULL

    # Remove the name of coef.str from the contrasts
    if (remove_name) {
        res$contrast <- gsub(coef.str, "", res$contrast)
    }

    # Adjust p-value across all contrasts
    res <- res %>%
        mutate(adj.P.Val = p.adjust(P.Value, method = adjust.method)) %>%
        select(feature, contrast, everything())
}

