#' Normalization of LC-MS/MS Data
#'
#' Fits a nonparametric trend into relative
#' abundance data. Any significant (at this point just any)
#' trend is corrected.
#' Converts one ExpressionSet/MSnSet to another ExpressionSet/MSnSet.
#'
#' @param eset ExpressionSet/MSnSet object
#' @param property character the column in the fData that the
#'              relative intensities regress against.
#' @param method character either "lowess" or "loess" at this point.
#' @param ... passed to \code{lowess} or \code{loess}.
#'
#' @note So far the only property I have in mind is elution time in
#'      label-free LC-MS/MS data.
#'
#' @importFrom Biobase exprs fData sampleNames exprs<-
#' @importFrom stats lowess loess predict
#'
#' @export normalize_by_feature_property


normalize_by_feature_property <- function(eset, property,
                                          method=c("lowess","loess"), ...){
    # elution time is a typical problem in a label-free LC-MS/MS
    method <- match.arg(method)
    o <- order(fData(eset)[[property]])
    if(method == "lowess"){
        for(sample in sampleNames(eset)){
            m <- lowess(fData(eset)[[property]], exprs(eset)[,sample], ...)
            exprs(eset)[,sample] <- exprs(eset)[,sample] - m$y[order(o)]
        }
    }else if(method == "loess"){
        x <- cbind(fData(eset)[,property,drop=F], exprs(eset))
        for(sample in sampleNames(eset)){
            xx <- x[,c(property,sample)]
            colnames(xx) <- c('x','y')
            m <- loess(y ~ x, data=xx, ...)
            exprs(eset)[,sample] <- exprs(eset)[,sample] - predict(m, xx)
        }
    }else{
        warning("method is not defined")
    }
    return(eset)
}
