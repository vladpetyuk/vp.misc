

#' Concatenate MSnSet Objects by Features
#'
#' Combines MSnset objects by concatenating features. Samples must
#' be the same and maintained the same. Pheno data also must be the
#' same.
#'
#' @param ... MSnSet objects
#' @param tags vector of characters attached to feature names of
#'          corresponding MSnSet objects
#'
#' @return MSnSet object with concatenated features
#'
#' @importFrom Biobase pData featureNames exprs
#' @importFrom MSnbase MSnSet
#'
#' @importMethodsFrom MSnbase featureNames<-
#'
#' @export concatenate_by_features
#'
#' @examples
#' data(srm_msnset)
#' dim(msnset)
#' featureNames(msnset)
#' msnset3 <- msnset2 <- msnset
#' # Z-score transform
#' exprs(msnset2) <- sweep(exprs(msnset2), 1, apply(exprs(msnset2), 1, sd), '/')
#' # retaining just the sign
#' exprs(msnset3) <- ifelse(exprs(msnset3) > 0, 1, -1)
#' m <- concatenate_by_features(msnset, msnset2, msnset3,
#'                              tags=c("_ori","_zsc","_sgn"))
#' dim(m)
#' featureNames(m)


concatenate_by_features <- function(..., tags=NULL){
    mL <- list(...)
    if(is.null(tags)) tags <- seq_along(mL)
    mL <- lapply(seq_along(mL), function(x){
        featureNames(mL[[x]]) <-
            paste(featureNames(mL[[x]]), tags[x], sep='')
        return(mL[[x]])
    })
    m3 <- Reduce( function(x1, x2) {
        e2 <- rbind(exprs(x1), exprs(x2))
        stopifnot(identical(pData(x1), pData(x2)))
        p2 <- pData(x1)
        f2 <- data.frame(feature=rownames(e2), row.names=rownames(e2))
        x3 <- MSnSet(exprs=e2, fData=f2, pData=p2)
    }, mL)
    return(m3)
}

