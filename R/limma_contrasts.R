#' Wrapper for Limma
#'
#' A convenience wrapper for limma contrast testing.
#'
#' @param eset eset (or most likely eset subclass) object
#' @param model.str character formulation of the model (e.g. "~ 0 + a + b")
#' @param coef.str vector of character coefficients of interest from
#'      \code{colnames(pData(eset))}. E.g. "a" or c("a", "b").
#' @param ref.str (character) Reference string for pairwise comparisons
#' @param covariates character vector of potential covariates that
#'      should be included in the model.
#' @param contrasts character vector of contrasts to test.
#' @param contrast.sep string used to join the level labels of the columns
#'      specified by \code{coef.str}. The default is an underscore, so factor A
#'      with entries \code{c("a1", "a1", "a2")} and factor B with entries
#'      \code{c("b1", "b2", "b3")} would combine to form a factor with entries
#'      \code{c("a1_b1", "a1_b2", "a2_b3")}.
#' @param adjust.method method for p-value adjustment. Default is "BH"
#' (Benjamini-Hochberg).
#' @param ... arguments for \code{\link[limma]{lmFit}}
#' 
#' @return data.frame. Basically output of \code{\link[limma]{topTable}}
#'      function with additional columns \code{feature} and \code{contrast}.
#'      
#' @importFrom Biobase exprs pData
#' @importFrom limma lmFit topTable eBayes contrasts.fit makeContrasts
#' @importFrom dplyr mutate select %>%
#' @importFrom tidyr everything
#' @export limma_contrasts limma_contrasts limma_contrasts_ref
#' @examples
#' 
#' library(MSnSet.utils)
#' data(cptac_oca)
#'
#' coef.str <- c("SUBTYPE", "SURVIVALSTATUS")
#' covariates <- c("AGE", "Batch")
#' 
#' # All combinations of coef.str
#' contrasts <- levels(interaction(pData(oca.set)[, coef.str],
#'                     sep = "_"))
#' 
#' res <- limma_contrasts(oca.set,
#'                        coef.str = coef.str,
#'                        covariates = covariates,
#'                        contrasts = contrasts)
#' 
#' head(res)

limma_contrasts <- function(eset, coef.str, covariates = NULL, contrasts,
                            contrast.sep = "_", adjust.method = "BH", ...) {

    pData(eset)[["coef.new"]] <- interaction(pData(eset)[, coef.str],
                                             sep = contrast.sep)

    # Model without intercept. Must include the "0".
    model.str <- paste("~ coef.new", covariates, "0", sep = " + ")

    # Replace any "+  +" with "+"
    model.str <- gsub("\\+  \\+", "\\+", model.str)

    model.formula <- eval(parse(text = model.str),
                          envir = pData(eset))

    design <- model.matrix(model.formula)

    # R doesn't like when column names begin with anything
    # other than a character, so we need to make sure that
    # each level in contrasts begins with "coef.new". We will
    # remove "coef.new" from the contrasts during the lapply step.
    for (i in 1:nlevels(pData(eset)[["coef.new"]])) {
        contrasts <- gsub(levels(pData(eset)[["coef.new"]])[i],
                          colnames(design)[i], contrasts)
    }

    contrast.matrix <- makeContrasts(contrasts = contrasts,
                                     levels = design)

    # This line is in limma_gen and limma_contrasts,
    # but it doesn't seem to do anything:
    eset <- eset[, as.numeric(rownames(design))]
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
                                 contrast = gsub("coef.new", "", contrast_i))
                  })

    # Bind topTable output
    res <- do.call(rbind, res)
    rownames(res) <- NULL
    
    # Adjust p-value across all contrasts
    res <- res %>%
        mutate(adj.P.Val = p.adjust(P.Value, method = adjust.method)) %>%
        select(feature, contrast, everything())
}

limma_contrasts_ref <- function(eset, model.str, coef.str,
                                ref.str = NULL, adjust.method = "BH", ...) {
    model.formula <- eval(parse(text = model.str), envir = pData(eset))
    design <- model.matrix(model.formula)
    
    eset <- eset[, as.numeric(rownames(design))]
    fit <- lmFit(exprs(eset), design, ...)
    fit.smooth <- eBayes(fit)
    
    if (!is.null(ref.str)) {
        eset[[coef.str]] <- relevel(as.factor(eset[[coef.str]]), ref=ref.str)
    }
    
    contrasts <- paste0(colnames(design)[1], "-", colnames(design))
    contrasts <- contrasts[-1]
    args <- c(as.list(contrasts), list(levels=design))
    contrast.matrix <- do.call(makeContrasts, args)
    
    fit2 <- contrasts.fit(fit.smooth, contrast.matrix, ...)
    fit2.smooth <- eBayes(fit2)
    
    sig <- lapply(colnames(contrast.matrix),
                  function(coef) {
                      sig <- topTable(fit2.smooth,
                                      number = nrow(eset), 
                                      sort.by = "none",
                                      coef = coef)
                      sig$feature <- rownames(sig)
                      rownames(sig) <- NULL
                      sig <- sig %>%
                          mutate(contrast = gsub(coef.str, "", coef)) %>%
                          select(contrast, feature, everything())
                  })
    sig <- do.call(rbind, sig)
    
    res <- res %>%
        mutate(adj.P.Val = p.adjust(P.Value, method = adjust.method))
}


