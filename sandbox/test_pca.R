
# plot_pca_v1


library("vp.misc")
# data(package="vp.misc")
data(srm_msnset)



plot_pca_v1 <- function(msnset, phenotype=NULL, show.ellispe=TRUE){
    
    library("ade4")
    px <- dudi.pca(exprs(msnset), scannf = F, scale=T, center=T, nf = 2)
    px <- px$co
    colnames(px) <- c("PC1", "PC2")
    
    if(!is.null(phenotype)){
        colorBy <- pData(msnset)[[phenotype]]
    }else{
        colorBy <- 1
        phenotype <- ""
        show.ellispe <- FALSE
    }
    ggdata <- data.frame(px, colorBy)
    
    library("ggplot2")
    p <- 
    ggplot(ggdata) +
        geom_point(aes(x=PC1, y=PC2, color=factor(colorBy)), 
                   size=5, shape=20, show_guide = TRUE) +
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


plot_pca_v1(msnset, phenotype = "subject.type", show.ellispe = F)
plot_pca_v1(msnset)


