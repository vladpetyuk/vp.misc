


#' Customized Distance Function
#'
#' Based on \code{\link[stats]{dist}}. Additionally
#' it can handle correlation-based distances and
#' has argument for imputation of NA distances. To
#' avoid confusion with \code{stats::dist} all the
#' default settings are such that it performs equivalently
#' to the aforementioned one.
#'
#' @param x a numeric matrix
#' @param method the distance measure to be used.
#'          This must be one of "euclidean", "maximum", "manhattan",
#'                              "canberra", "binary", "minkowski",
#'                              "pearson", "kendall", "spearman".
#'          Any unambiguous substring can be given.
#' @param na.avg logical if yes, then NA distance will
#'                      be imputed with average distance
#'                      value. Default is FALSE.
#' @param ... arguments of \code{\link[stats]{dist}}
#'
#' @note In the current set-up
#'      corrlation values -1, 0 and 1 correspond to
#'      distances 1, 0.5 and 1, respectively.\cr
#'      R=1  -> dist=0 \cr
#'      R=0  -> dist=0.5 \cr
#'      R=-1 -> dist=1 \cr
#'
#' @importFrom stats as.dist cor
#'
#' @examples
#' x <- matrix(c(1:5, 5:9, 9:13), ncol=3)
#' stats::dist(x)
#' dist(x)
#' dist(x, method="pearson")
#'
#' @export

dist <- function(x, method="euclidean", na.avg=FALSE, ...){
    cor.dists <-c("pearson", "kendall", "spearman")
    tru.dists <- c("euclidean", "maximum", "manhattan",
                   "canberra", "binary", "minkowski")
    method <- match.arg(method, c(cor.dists, tru.dists))
    if(method %in% cor.dists){
        res <- (1-cor( t(x),
                       method=method,
                       use="pairwise.complete.obs" ))/2
    }else{
        res <- stats::dist(x, ...)
    }
    if(na.avg){
        res[is.na(res)] <- mean(res, na.rm=T)
    }
    return(as.dist(res))
}









# make like dist.nona that same as dist, except if dist NA it assumes avg dist.
# also it intentially masks stats::dist and adds method pearson.


# dist_nona <- function(x, ...){
#     distMat <- dist(x, ...)
#     distMat[is.na(distMat)] <- mean(distMat, na.rm=T)
#     return(as.dist(distMat))}
# }




# pearson.dist <- function(x){
#     corMat <- (1-cor( t(x),
#                       method="pearson",
#                       use="pairwise.complete.obs" ))/2
#     corMat[is.na(corMat)] <- mean(corMat, na.rm=T)
#     return(as.dist(corMat))}
#
# eucl.dist <- function(x){
#     distMat <- dist(x)
#     distMat[is.na(distMat)] <- mean(distMat, na.rm=T)
#     return(as.dist(distMat))}



