#' Wrapper for Limma
#' 
#' A convenience wrapper for limma
#' 
#' @param eset eset (or most likely eset subclass) object
#' @param model.str character formulation of the model (e.g. "~ a + b")
#' @param coef.str character coefficient of interest. E.g. "a".
#'              Either numeric or factor with 2 levels. 
#' @return data.frame. Basically output of 
#'      \code{\link[limma]{topTable}} function.
#' @importFrom Biobase exprs pData
#' @importFrom limma lmFit topTable eBayes
#' @export limma_a_b
#' @examples
#' library("vp.misc")
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

limma_a_b <- function(eset, model.str, coef.str){
    
    model.formula <- eval(parse(text=model.str), envir=pData(eset))
    design <- model.matrix( model.formula)

    coef.str <- grep(coef.str, colnames(design), value = TRUE)
    if(length(coef.str) > 1){
        stop("Multiple matches for coefficient of interest!")
        stop("Most likely it is factor with > 2 levels.")
    }
    eset <- eset[,as.numeric(rownames(design))]
    fit <- lmFit(exprs(eset), design, method="robust", maxit=10000)
    fit.smooth <- eBayes(fit)
    sig <- topTable(fit.smooth, number=nrow(eset), 
                    sort.by='none', coef=coef.str)
    return(sig)
}


