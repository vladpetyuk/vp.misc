


#' Normalization of SRM Data
#'
#' Subtracts expression data of one condition from another.
#' Converts one MSnSet to another MSnSet.
#'
#' @param m MSnSet object
#' @param exprsn expression with proteins as terms.
#'          It is evaluated and used as reference.
#'
#' @importFrom Biobase exprs exprs<-
#'
#' @export normalizeByRef

normalizeByRef <- function(m, exprsn){
    #
    ref.call <- substitute(exprsn)
    ref <- eval(ref.call, envir = data.frame(t(exprs(m))))
    exprs(m) <- sweep(exprs(m), 2, ref, '-') # must be log scale
    return(m)

}

# normalizeByRef(m, MBP.1)
# normalizeByRef(m, TH.2)
# normalizeByRef(m, (MBP.1+TH.2)/2)

