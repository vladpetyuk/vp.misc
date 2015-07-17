
#' Basic Normalization Routines
#' 
#' Normalizes data based on assumption that on average peptide deviation 
#' is zero.  The approach is typically used for global shotgun proteomics
#' data. Converts one MSnSet to another MSnSet.
#' 
#' @param m ExpresionSet or MSnSet object
#' @param method character. "medpolish" - uses median polish algorithm.
#'          "once" means just correct for median differences in the samples.
#'          Note, doing it just "once" may de-center the features (rows).
#' 
#' @export normalizeByGlob
#' @examples
#' data(srm_msnset)
#' image_msnset(msnset)
#' msnset2 <- normalizeByGlob(msnset)
#' image_msnset(msnset2)

normalizeByGlob <- function(m, method=c("medpolish","once")){
    #
    method <- match.arg(method)
    #
    if(method == "once"){
        mc <- m[complete.cases(exprs(m)),]
        sample.bias <- apply(exprs(mc), 2, median, na.rm=T)
        exprs(m) <- sweep(exprs(m), 2, sample.bias, '-')
    }else if(method == "medpolish"){
        out <- medpolish(exprs(m), eps = .Machine$double.eps, 
                         maxiter = 100, na.rm = T, trace.iter=F)
        exprs(m) <- out$residuals
    }
    return(m)
}
