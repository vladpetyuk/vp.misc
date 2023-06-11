#' Wrapper for LIMMA
#'
#' Convenience wrappers for differential analysis with LIMMA.
#'
#' @param eset \code{eSet} (or most likely eset subclass) object
#' @param model.str character; formulation of the model (e.g. \code{"~ a + b"})
#' @param coef.str character; coefficient of interest. One of
#'      \code{colnames(pData(eset))}. E.g. \code{"a"}.
#'     For \code{limma_a_b}, this must be a character or factor with 2 levels
#'     or a numeric vector. For \code{limma_gen}, this can be a character or
#'     factor with at least 2 levels or a numeric vector.
#' @param ... arguments passed to \code{\link[limma]{lmFit}}
#'
#' @details
#' \code{limma_a_b} is for either 2-factor ANOVA or linear regression.
#' \code{limma_gen} is a more generic version and will return results
#' for all the coefficients that match \code{coef.str} pattern. It also
#' works the same as \code{limma_a_b} for linear regression.
#'
#' If performing linear regression, the \code{logFC} column contains the slopes
#' of the regression lines and the \code{AveExpr} column contains the
#' y-intercepts.
#'
#' \code{limma_a_b} outputs moderated t-statistics and associated p-values.
#' \code{limma_gen} outputs moderated F-statistics and associated p-values
#' for a set of contrasts. Also, the columns that begin with the
#' \code{coef.str} report the log2 fold-changes. Contrasts are always relative
#' to a reference level (if \code{coef.str} is a factor) or the first unique
#' entry of \code{coef.str} in terms of alphanumeric order.
#'
#' \code{\link[MSnSet.utils]{limma_contrasts}} is another LIMMA wrapper for
#' testing custom contrasts. It has a separate documentation page.
#'
#' @return
#' \code{data.frame}.
#' Basically output of the \code{\link[limma]{topTable}} function.
#'
#' @importFrom Biobase exprs pData
#' @importFrom limma lmFit topTable eBayes
#' @importFrom stats model.matrix
#'
#' @export limma_a_b
#'
#' @seealso
#' \code{\link{limma_contrasts}},
#' \code{\link[limma]{topTable}},
#' \code{\link[limma]{lmFit}}
#'
#' @examples
#' library("MSnSet.utils")
#' data("cptac_oca")
#'
#' # Filter to features detected in at least 30 samples
#' ee <- oca.set
#' ee <- ee[rowSums(!is.na(exprs(ee))) >= 30, ]
#'
#' # Note: only two levels are allowed in this case
#' res1 <- limma_a_b(ee,
#'                   model.str = "~ PLATINUM.STATUS + AGE",
#'                   coef.str = "PLATINUM.STATUS")
#' head(res1)
#'
#' # Linear regression
#' res2 <- limma_a_b(ee, model.str = "~ AGE", coef.str = "AGE")
#' head(res2)
#'
#' # Multiple levels - "Immunoreactive" is used as the reference level
#' res3 <- limma_gen(ee,
#'                   model.str = "~ SUBTYPE",
#'                   coef.str = "SUBTYPE")
#' head(res3)
#'



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

#' @export
#' @describeIn limma_a_b More generalized approach. The tested variable can have
#'                       more than 2 groups, unlike \code{limma_a_b}.
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
    idx <- which(names(attr(design, "contrast")) == gsub("`","",coef.str))
    idx <- attr(design, "assign") == idx
    coef.str <- colnames(design)[idx]
  }

  eset <- eset[,as.numeric(rownames(design))]
  fit <- lmFit(exprs(eset), design, ...)
  fit.smooth <- eBayes(fit)
  sig <- topTable(fit.smooth, number=nrow(eset),
                  sort.by='none', coef=coef.str)
  se <- (sqrt(fit.smooth$s2.post) * fit.smooth$stdev.unscaled)[,coef.str,drop=FALSE]
  colnames(se) <- paste("SE.", colnames(se), sep = "")
  result <- cbind(sig, se)
  return(result)
}

