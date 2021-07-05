#' Principal Component Reduction
#' 
#' Reduce an ExpressionSet object down to its principal components.  
#' 
#' @param eset Expression-Set object or similar
#' 
#' @examples
#' data("cptac_oca")
#' oca.pc <- pc_reduction(oca.set)
#' 
#' res <- limma_a_b(oca.pc, 
#'                  model.str = "~ PLATINUM.STATUS + AGE", 
#'                  coef.str = "PLATINUM.STATUS")
#' 
#' @export

pc_reduction <- function(eset) {
  
  eset_cc <- eset[complete.cases(exprs(eset)),]
  
  mat <- t(exprs(eset_cc))
  
  # z-score transform
  mat <- sweep(mat, 1, rowMeans(mat), "-")
  mat <- sweep(mat, 1, apply(mat, 1, sd), "/")
  
  pca <- prcomp(x = mat, scale. = F)
  loadings <- pca$rotation # loading scores
  scores <- t(pca$x) # scores
  
  loadings <- split(loadings, col(loadings))
  names(loadings) <- colnames(pca$rotation)
  
  # Set variable loading names
  for (i in seq_along(loadings)) {
    names(loadings[[i]]) <- rownames(pca$rotation)
  }
  
  feat_data <- data.frame(loading_scores = I(loadings))
  
  pheno_data <- pData(eset)[colnames(scores),]
  
  eset_pca <- MSnSet(scores, feat_data, pheno_data)
  
  return(eset_pca)
  
}

