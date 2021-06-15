#' Wrapper for Limma
#'
#' A convenience wrapper for limma contrast testing.
#'
#' @param eset eset (or most likely eset subclass) object
#' @param coef.str vector of character coefficients of interest from
#'      \code{colnames(pData(eset))}. E.g. "a" or c("a", "b").
#' @param covariates character vector of potential covariates that
#'      should be included in the model.
#' @param contrasts character vector of contrasts to test.
#' @param contrast.sep string used to join the level labels of the columns
#'      specified by \code{coef.str}. The default is an underscore, so factor A
#'      with entries \code{c("a1", "a1", "a2")} and factor B with entries
#'      \code{c("b1", "b2", "b3")} would combine to form a factor with entries
#'      \code{c("a1_b1", "a1_b2", "a2_b3")}.
#' @param ... arguments for \code{\link[limma]{lmFit}}
#' @return data.frame. Basically output of \code{\link[limma]{topTable}}
#'      function with additional columns \code{feature} and \code{contrast}.
#' @importFrom Biobase exprs pData
#' @importFrom limma lmFit topTable eBayes contrasts.fit makeContrasts
#' @importFrom dplyr mutate %>%
#' @importFrom tidyr everything
#' @export limma_contrasts_2
#' @examples
#' library("MSnSet.utils")
#' data("cptac_oca", package = "MSnSet.utils")
#' ee <- oca.set
#' ee <- ee[rowSums(!is.na(exprs(ee))) >= 30, ]
#'
#' coef.str <- c("SUBTYPE", "SURVIVALSTATUS")
#' covariates <- c("AGE", "Batch")
#' # All combinations of coef.str
#' contrasts <- levels(interaction(pData(ee)[, coef.str],
#'                     sep = "_"))
#'
#' res <- limma_contrasts_2(ee,
#'                          coef.str = coef.str,
#'                          covariates = covariates,
#'                          contrasts = contrasts)
#'
#' head(res)

limma_contrasts_2 <- function(eset,
                              coef.str,
                              covariates = NULL,
                              contrasts,
                              contrast.sep = "_",
                              ...) {

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
                                 contrast = gsub("coef.new",
                                                 "",
                                                 contrast_i))
                  }
    )

    # Bind topTable output
    res <- do.call(rbind, res)

    # Adjust p-value across all contrasts
    res$adj.P.Val <- p.adjust(res$P.Value, method = "BH")

    # Filter rows
    res <- res[res$adj.P.Val < 0.05, ]
    res <- res[!is.na(res$feature), ]

    # Reorder columns
    res <- res[, c(7, 8, 1:6)]

    return(res)
}

