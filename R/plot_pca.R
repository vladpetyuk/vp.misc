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
#' @importFrom ade4 "dudi.pca"
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
                   size=5, shape=20, show.legend = TRUE) +
        coord_fixed() +
        guides(color=guide_legend(phenotype),
               fill=guide_legend(phenotype)) +
        theme_bw()
    if(show.ellispe){
        p <- p +
        stat_ellipse(aes(x=PC1, y=PC2, fill=factor(colorBy)),
                     geom="polygon", type="norm", 
                     level=0.5, alpha=0.1, show.legend = TRUE)
    }
    return(p)
}


#' @describeIn plot_pca_v1 Alternative PCA
#' @importFrom pcaMethods pca
#' @importFrom RColorBrewer "brewer.pal"
#' @export plot_pca_v2
#' 
#' @examples
#' data(srm_msnset)
#' plot_pca_v2(msnset, phenotype = "subject.type")
#' plot_pca_v2(msnset)

plot_pca_v2 <- function(eset, phenotype=NULL, names=FALSE){
    
    # library("pcaMethods")
    # library("RColorBrewer")
    
    pcaResults <- pca(exprs(eset))
    plot(pcaResults@loadings*1.4, pch='')
    colz <- 'black'
    if(!is.null(phenotype)){
        coloring <- as.factor(pData(eset)[[phenotype]])
        colScheme <- brewer.pal(max(3, nlevels(coloring)), "Set1")
        text.col <- colScheme[seq_len(nlevels(coloring))]
        legend("topleft", 
               legend=levels(coloring), 
               text.col=text.col,
               bg="#EEEEEE80")
        colz <- colScheme[coloring]
        # title(sprintf("PCA colored by %s", "WHIM sample type"))
    }
    
    text(pcaResults@loadings,
         sampleNames(eset), 
         cex=0.75, font=2, col=colz, adj=c(0.25,0))
    grid()
}

# # I don't like this way of combining
# #
# #' @describeIn plot_pca_v1
# #' @importFrom made4 ord plotarrays
# # ' @importFrom ade4 dudi.coa dudi.pca
# # ' @export plot_pca_v3
# #' @examples 
# #'
# #' plot_pca_v3(msnset, type='pca', phenotype="subject.type")
# #' plot_pca_v3(msnset, type='coa', phenotype="subject.type")
# 
# plot_pca_v3 <- function(eset, phenotype=NULL, ...){
#     phenotype <- as.factor(pData(eset)[[phenotype]])
#     ord.res <- ord(exprs(eset), classvec=phenotype, ...)
#     plotarrays(ord.res, ...)
# }
# 


