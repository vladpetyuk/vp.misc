


#' Normalization of SRM Data.
#' 
#' Normalizes data based on assumption that on average peptide deviation 
#' is zero.  The approach is typically used for global shotgun proteomics
#' data. Converts one MSnSet to another MSnSet. Uses median polish algorithm.
#' 
#' @param m MSnSet object
#' 
#' @export normalizeByGlob

normalizeByGlob <- function(m){
    #
    # sample.bias <- apply(exprs(m), 2, median, na.rm=T)
    # exprs(m) <- sweep(exprs(m), 2, sample.bias, '-')
    out <- medpolish(exprs(m), eps = .Machine$double.eps, 
                     maxiter = 100, na.rm = T, trace.iter=F)
    exprs(m) <- out$residuals
    return(m)
    
}

# normalizeByRef(m, MBP.1)
# normalizeByRef(m, TH.2)
# normalizeByRef(m, (MBP.1+TH.2)/2)

# x <- medpolish(exprs(m), eps = .Machine$double.eps, maxiter = 100, na.rm = F)
# plot(apply(x$residuals, 2, median, na.rm=T))
# plot(apply(x$residuals, 1, median, na.rm=T))
