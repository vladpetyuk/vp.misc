#' Wrapper for Limma
#'
#' A convenience wrapper for limma
#'
#' @param eset eset (or most likely eset subclass) object
#' @param model.str character formulation of the model (e.g. "~ a + b")
#' @param coef.str character coefficient of interest. E.g. "a".
#'              Either numeric or factor with 2 levels.
#' @param ... arguments for \code{\link[limma]{lmFit}}
#'
#' @return data.frame. Basically output of
#'      \code{\link[limma]{topTable}} function.
#'
#' @note limma_a_b is for either 2-factor ANOVA or linear regression.
#'       limma_gen is a more generic version and will return results
#'       for all the coefficients that match \code{coef.str} pattern.
#'
#' @importFrom Biobase exprs pData
#' @importFrom limma lmFit topTable eBayes
#' @importFrom stats model.matrix
#'
#' @export limma_a_b limma_gen limma_contrasts
#'
#' @examples
#' library("MSnSet.utils")
#' data("cptac_oca")
#' library("limma")
#' ee <- oca.set
#' ee <- ee[rowSums(!is.na(exprs(ee))) >= 30,]
#' model.string <- "~ PLATINUM.STATUS + AGE" # model with covariates
#' coef.string <- "PLATINUM.STATUS" # note, only two levels allowed here
#'
#' res <- limma_a_b(oca.set,
#'                  model.str = "~ PLATINUM.STATUS + AGE",
#'                  coef.str = "PLATINUM.STATUS")
#' head(res)



limma_a_b <- function(eset, model.str, coef.str, ...){

    model.formula <- eval(parse(text=model.str), envir=pData(eset))
    design <- model.matrix( model.formula)

    coef.str <- grep(coef.str, colnames(design), value = TRUE)
    if(length(coef.str) > 1){
        stop("Multiple matches for coefficient of interest!")
        stop("Most likely it is factor with > 2 levels.")
    }

    eset <- eset[,as.numeric(rownames(design))]
    fit <- lmFit(exprs(eset), design, ...)
    fit.smooth <- eBayes(fit)
    sig <- topTable(fit.smooth, number=nrow(eset),
                    sort.by='none', coef=coef.str)
    return(sig)
}


#' @describeIn limma_a_b More generalized approach. The tested variable doesn't
#'                       have to be a factor with 2 levels like in limma_a_b.
limma_gen <- function(eset, model.str, coef.str, ...){

    model.formula <- eval(parse(text=model.str), envir=pData(eset))
    design <- model.matrix( model.formula)

    # this malfunctions if coefficient of interest coincided with
    # suffix of some of the covariates. E.g. Plate and PlateCol.
    # Testing for Plate is problematic
    # coef.str <- grep(coef.str, colnames(design), value = TRUE)

    # a new way
    # If coef.str is a factor or character, do this
    if (!(coef.str %in% colnames(design))) {
        idx <- which(names(attr(design, "contrast")) == coef.str)
        idx <- attr(design, "assign") == idx
        coef.str <- colnames(design)[idx]
    }

    eset <- eset[,as.numeric(rownames(design))]
    fit <- lmFit(exprs(eset), design, ...)
    fit.smooth <- eBayes(fit)
    sig <- topTable(fit.smooth, number=nrow(eset),
                    sort.by='none', coef=coef.str)
    return(sig)
}

