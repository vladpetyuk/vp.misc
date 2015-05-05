
#' Customized K-NN Imputation
#' 
#' K-NN imputation adopted to LC-MS proteomics data.
#' The main reason for missing data in LC-MS datasets is
#' low abundance of the protein/peptides. Therefore this 
#' K-NN imputation algorithm explicitely relies on this
#' assumption.
#' 
#' The algorithm. For each row in the exprs matrix, 
#' that contain missing values, perform the following steps:
#' \enumerate{
#'  \item impute missing values with the lowest values in the row (feature)
#'  \item find K features (with no missing values) with 
#'          highest Spearman correlation
#'  \item scale the K-neighbors, so that median intensity ratio is 1
#'  \item impute missing values with mean value of scaled K-neighbors
#' }
#' 
#' @param x MSnSet or ExpressionSet object
#' @param K number of nearest neighbors
#' @param show.diagnostics logical indicating if to 
#'                  plot the results of imputation for each feature
#' 
#' @note The algorithm assumes that the data is not log-transformed. Thus, if
#'          the data is log-transform - exponentiate.
#' @examples
#' suppressPackageStartupMessages(library("MSnbase"))
#' data(naset)
#' image(naset[1:50,])
#' x <- knn_imputation(naset)
#' image(x[1:50,])
#' 
#' @export


knn_imputation <- function( x, K=10, show.diagnostics=F){   
    M <- exprs(x)
    M.original <- M
    # first extract entries with no missing values
    neighborhood <- M[ !apply(is.na(M),1,any), ]
    for( row_i in 1:nrow(M)){
        if( any(is.na(M[row_i,]))){
            imputed <- impute_one( M[row_i,], neighborhood, K)
            M[row_i,] <- imputed$vec
            if( show.diagnostics)
                plot_diagnostics( M.original[row_i,], 
                                  imputed$vec, 
                                  imputed$knn.transformed)
        }
    }
    exprs(x) <- M
    return(x)
}



get_knn <- function( vec, neighborhood, K=10){
   # function that return K nearest neighbors
   # vec - vector
   # neighborhood - set to select nearest neighbors. This is a selection with no NAs
   # K - number of nearest neighbors to pull out
   #
   # impute vector, so there are no missing values
   # note! this assumes that missingness is due to low abundance
   vec.imputed = vec
   vec.imputed[is.na(vec.imputed)] = min(vec.imputed, na.rm=T)
   #
   get_distance <- function(x){
      distMatr = cor( as.numeric(x), as.numeric(vec.imputed), 
                      method="spearman", use="pairwise.complete.obs")
   }
   #
   distances = apply( neighborhood, 1, get_distance)
   knn.names = names(rev(sort(distances))[1:K])
   knn = neighborhood[knn.names,]
   return(knn)
}



impute_one <- function(vec, neighborhood, K=10){
   #
   knn = get_knn(vec, neighborhood, K)
   #
   # linTrans  <- t(apply( knn, 1,
                         # function(x) lm( as.numeric(vec) ~ as.numeric(x) - 1)$coefficients))     
   linTrans = 1/apply( t(knn)/as.numeric(vec),2,median, na.rm=T) # scaling, intercept is zero
   #
   dSetTrans <- knn*linTrans #[,2] + linTrans[,1] - this is in case of a+b*x model
   na_i = is.na(vec)
   vec[na_i] = apply( dSetTrans, 2, mean, na.rm=T)[na_i]
   #
   return(list(vec=vec, knn=knn, knn.transformed=dSetTrans))   
}



plot_diagnostics <- function( vec, vec.imputed, knn.transformed ){
   #
   matplot( t(knn.transformed), type="b")
   matlines( t(vec.imputed), type="b", lwd=5, col="black")
   matpoints( t(vec), type="b", pch=19, cex=2, col="green")          
   vec.imputed.only = vec.imputed
   vec.imputed.only[!is.na(vec)] = NA
   matpoints( t(vec.imputed.only), type="b", pch=19, cex=2, col="red")          
}


