#' PCA Plot
#' 
#' A convenience function for plotting PCA scatter plot
#' for samples in ExpressionSet/MSnSet object.
#' 
#' @param eset eset (or most likely eset subclass) object
#' @param phenotype character one of the \code{colnames(pData(eset))}
#' @param show.ellipse logical determining to plot 95\% CI based on 
#'          Hotelling's T-test or not.
#' @return plot
#' @importFrom Biobase exprs pData
#' @importFrom ade4 dudi.pca
#' @importFrom ggplot2 ggplot geom_point coord_fixed theme_bw
#'                      guides guide_legend stat_ellipse aes
#' @export plot_pca_v1
#' 
#' @examples
#' data(srm_msnset)
#' plot_pca_v1(msnset, phenotype = "subject.type", show.ellispe = F)
#' plot_pca_v1(msnset, phenotype = "subject.type", show.ellispe = T)
#' plot_pca_v1(msnset)

plot_pca_v1 <- function(eset, phenotype=NULL, show.ellispe=TRUE){
    
    # get rid of NA values
    stopifnot(sum(complete.cases(exprs(eset))) > 1)
    eset <- eset[complete.cases(exprs(eset)),]

    # PCA itself
    px <- dudi.pca(exprs(eset), scannf = F, scale=T, center=T, nf = 2)
    px <- px$co
    colnames(px) <- c("PC1", "PC2")
    
    if(!is.null(phenotype)){
        colorBy <- pData(eset)[[phenotype]]
    }else{
        colorBy <- 1
        phenotype <- ""
        show.ellispe <- FALSE
    }
    ggdata <- data.frame(px, colorBy)
    
    # visualize
    p <- 
    ggplot(ggdata) +
        geom_point(aes(x=PC1, y=PC2, color=factor(colorBy)), 
                   size=5, shape=20, show_guide = TRUE) +
        coord_fixed() +
        guides(color=guide_legend(phenotype),
               fill=guide_legend(phenotype)) +
        theme_bw()
    if(show.ellispe){
        p <- p +
        stat_ellipse(aes(x=PC1, y=PC2, fill=factor(colorBy)),
                     geom="polygon", type="norm", 
                     level=0.5, alpha=0.1, show_guide = TRUE)
    }
    plot(p)
}


plot_pca_v1 <- function(eset, phenotype=NULL, show.ellispe=TRUE, names=FALSE){
    
    # get rid of NA values
    stopifnot(sum(complete.cases(exprs(eset))) > 1)
    eset <- eset[complete.cases(exprs(eset)),]
    
    # PCA itself
    # scale=T, center=T to arguments
    px <- dudi.pca(exprs(eset), scannf = F, scale=T, center=T, nf = 2)
    px <- px$co
    colnames(px) <- c("PC1", "PC2")
    
    if(!is.null(phenotype)){
        colorBy <- pData(eset)[[phenotype]]
    }else{
        colorBy <- 1
        phenotype <- ""
        show.ellispe <- FALSE
    }
    ggdata <- data.frame(px, colorBy)
    
    # visualize
    p <- 
        ggplot(ggdata) +
        geom_point(aes(x=PC1, y=PC2, color=factor(colorBy)), 
                   size=5, shape=20, show_guide = TRUE) +
        coord_fixed() +
        guides(color=guide_legend(phenotype),
               fill=guide_legend(phenotype)) +
        theme_bw()
    if(show.ellispe){
        p <- p +
            stat_ellipse(aes(x=PC1, y=PC2, fill=factor(colorBy)),
                         geom="polygon", type="norm", 
                         level=0.5, alpha=0.1, show_guide = TRUE)
    }
    plot(p)
}



# data(srm_msnset)
# plot_pca_v1(msnset, phenotype = "subject.type", show.ellispe = F)
# plot_pca_v1(msnset, phenotype = "subject.type", show.ellispe = T)
# plot_pca_v1(msnset)
 
 