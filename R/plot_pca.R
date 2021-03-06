#' PCA Plot
#' 
#' A convenience function for plotting PCA scatter plot
#' for samples in ExpressionSet/MSnSet object.
#' 
#' @param eset eset (or most likely eset subclass) object
#' @param phenotype character one of the \code{colnames(pData(eset))}
#' @param show.NA logical Should the datapoints for which phenotype is 
#'          unknown be shown? Default is TRUE
#' @param show.ellipse logical determining to plot 95\% CI based on 
#'          Hotelling's T-test or not.
#' @param legend.title.width integer Wrapping up too long legend titles.
#'          Passed to stringr::str_wrap as width argument. Default is 20.
#' @return plot
#' @importFrom Biobase exprs pData
#' @importFrom ade4 "dudi.pca"
#' @importFrom ggplot2 ggplot geom_point coord_fixed theme_bw
#'                      guides guide_legend stat_ellipse aes
#' @importFrom stringr str_wrap
#' @export plot_pca_v1
#' 
#' @examples
#' data(srm_msnset)
#' plot_pca_v1(msnset, phenotype = "subject.type", show.ellispe = F)
#' plot_pca_v1(msnset, phenotype = "subject.type", show.ellispe = T)
#' plot_pca_v1(msnset)

plot_pca_v1 <- function(eset, phenotype=NULL, show.ellispe=TRUE, 
                        show.NA=TRUE, legend.title.width=20){
    
    # handling coloring by phenotype
    if (!is.null(phenotype)) {
        colorBy <- pData(eset)[[phenotype]]
        if(!show.NA){
            idx <- !is.na(colorBy)
            eset <- eset[,idx]
            colorBy <- colorBy[idx]
        }
    }
    else {
        colorBy <- ""
        phenotype <- ""
        show.ellispe <- FALSE
    }
    
    # get rid of NA values
    stopifnot(sum(complete.cases(exprs(eset))) > 1)
    eset <- eset[complete.cases(exprs(eset)),]

    # PCA itself
    px <- dudi.pca(exprs(eset), scannf = F, scale=T, center=T, nf = 2)
    px <- px$co # this is plotting loadings rather then scores!!!
    colnames(px) <- c("PC1", "PC2")
    
    # visualize
    ggdata <- data.frame(px, colorBy)
    p <- 
    ggplot(ggdata) +
        geom_point(aes(x=PC1, y=PC2, color=colorBy), 
                   size=5, shape=20, show.legend = TRUE) +
        coord_fixed() +
        theme_bw()
    
    # Ugly engtanglement of if/else statements. Needs to be improved.
    phenotype_str <- str_wrap(phenotype, legend.title.width)
    if(show.ellispe){
        p <- p +
            stat_ellipse(aes(x=PC1, y=PC2, fill=colorBy),
                         geom="polygon", type="norm", 
                         level=0.5, alpha=0.1, show.legend = TRUE) +
            guides(color=guide_legend(phenotype_str),
                   fill=guide_legend(phenotype_str))
    }else{
        if(is.numeric(colorBy)){
            p <- p + guides(color=guide_colorbar(phenotype_str))
        } else if (!identical(colorBy,"")){
            p <- p + guides(color=guide_legend(phenotype_str))
        }
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





#' @describeIn plot_pca_v1 Alternative PCA
#' @importFrom ggrepel geom_label_repel
#' @export plot_pca_v3
#' 
#' @examples
#' data(srm_msnset)
#' plot_pca_v3(msnset, phenotype = "subject.type")
#' plot_pca_v3(msnset, phenotype = "subject.type", label = "sample.id")
#' plot_pca_v3(msnset)

plot_pca_v3 <- function(eset, phenotype=NULL, label=NULL, point_size = 5, show.ellispe=TRUE, 
                        show.NA=TRUE, legend.title.width=20, ...){
    
    # handling coloring by phenotype
    if (!is.null(phenotype)) {
        colorBy <- pData(eset)[[phenotype]]
        if(!show.NA){
            idx <- !is.na(colorBy)
            eset <- eset[,idx]
            colorBy <- colorBy[idx]
        }
    }
    else {
        colorBy <- ""
        phenotype <- ""
        show.ellispe <- FALSE
    }
    
    # get rid of NA values
    stopifnot(sum(complete.cases(exprs(eset))) > 1)
    eset <- eset[complete.cases(exprs(eset)),]
    
    # PCA itself
    z <- t(exprs(eset))
    #  zero-center and scale because we more care about correlation 
    #  rather Eucledian dist
    z <- sweep(z, 1, rowMeans(z), FUN = "-")
    z <- sweep(z, 1, apply(z,1,sd), FUN = "/")
    pca1 <- prcomp(z, scale. = F)
    # create data frame with scores
    scores <- as.data.frame(pca1$x)
    
    exp_var <- 100 * summary(pca1)$importance[2,][c(1,2)] # from ggord
    axes <- paste0("PC", c(1,2))
    axes <- paste0(axes, " (", round(exp_var, 2), "%)")
    
    
    # visualize
    ggdata <- data.frame(scores[,c(1,2)], colorBy) # first two PCs
    p <- 
        ggplot(ggdata) +
        geom_point(aes(x=PC1, y=PC2, color=colorBy), 
                   size=point_size, shape=20, show.legend = TRUE) +
        coord_fixed() +
        xlab(axes[1]) + ylab(axes[2]) +
        theme_bw()
    
    
    # Add labels by 'label'
if (!is.null(label)) {
      custom_args <- list(mapping = aes(x = PC1, y = PC2, fill = colorBy, label = pData(eset)[[label]]),
                          fontface = 'bold', color = 'white',
                          box.padding = 0.25, point.padding = 0.25,
                          segment.color='grey50',
                          label.size = 0.01,
                          segment.alpha = 0.50,
                          size = 2.5)
      user_args <- list(...)
      custom_args[names(user_args)] <- user_args
      
      p <- p +
        do.call(geom_label_repel, custom_args)    
    }
    
    # Ugly engtanglement of if/else statements. Needs to be improved.
    phenotype_str <- str_wrap(phenotype, legend.title.width)
    if(show.ellispe){
        p <- p +
            stat_ellipse(aes(x=PC1, y=PC2, fill=colorBy),
                         geom="polygon", type="norm", 
                         level=0.5, alpha=0.1, show.legend = TRUE) +
            guides(color=guide_legend(phenotype_str),
                   fill=guide_legend(phenotype_str))
    }else{
        if(is.numeric(colorBy)){
            p <- p + guides(color=guide_colorbar(phenotype_str))
        } else if (!identical(colorBy,"")){
            p <- p + guides(color=guide_legend(phenotype_str))
        }
    }
    return(p)
}




#' @describeIn plot_pca_v1 Alternative PCA with option to select pair of principal components
#' @importFrom ggrepel geom_label_repel
#' @export plot_pca_v4
#' 
#' @examples
#' data(srm_msnset)
#' plot_pca_v4(msnset, phenotype = "subject.type", pc.index=c(1,3))
#' plot_pca_v4(msnset, phenotype = "subject.type", label = "sample.id", pc.index=c(1,3))
#' plot_pca_v4(msnset, pc.index=c(1,3))

plot_pca_v4 <- function (eset, phenotype = NULL, label = NULL, point_size = 5, 
                         show.ellispe = TRUE, show.NA = TRUE, legend.title.width = 20,
                         pc.index = c(1, 2),
                         ...) 
{
  if (!is.null(phenotype)) {
    colorBy <- pData(eset)[[phenotype]]
    if (!show.NA) {
      idx <- !is.na(colorBy)
      eset <- eset[, idx]
      colorBy <- colorBy[idx]
    }
  }
  else {
    colorBy <- ""
    phenotype <- ""
    show.ellispe <- FALSE
  }
  stopifnot(sum(complete.cases(exprs(eset))) > 1)
  eset <- eset[complete.cases(exprs(eset)), ]
  z <- t(exprs(eset))
  z <- sweep(z, 1, rowMeans(z), FUN = "-")
  z <- sweep(z, 1, apply(z, 1, sd), FUN = "/")
  pca1 <- prcomp(z, scale. = F)
  scores <- as.data.frame(pca1$x)
  exp_var <- 100 * summary(pca1)$importance[2, ][pc.index]
  axes <- paste0("PC", pc.index)
  axes <- paste0(axes, " (", round(exp_var, 2), "%)")
  ggdata <- data.frame(scores[, pc.index], colorBy)
  p <- ggplot(ggdata) + geom_point(aes(x = ggdata[,1], y = ggdata[,2], color = colorBy), 
                                   size = point_size, shape = 20, show.legend = TRUE) + 
    coord_fixed() + xlab(axes[1]) + ylab(axes[2]) + theme_bw()
  if (!is.null(label)) {
    custom_args <- list(mapping = aes(x = ggdata[,1], y = ggdata[,2], 
                                      fill = colorBy, label = pData(eset)[[label]]), fontface = "bold", 
                        color = "white", box.padding = 0.25, point.padding = 0.25, 
                        segment.color = "grey50", label.size = 0.01, segment.alpha = 0.5, 
                        size = 2.5)
    user_args <- list(...)
    custom_args[names(user_args)] <- user_args
    p <- p + do.call(geom_label_repel, custom_args)
  }
  phenotype_str <- str_wrap(phenotype, legend.title.width)
  if (show.ellispe) {
    p <- p + stat_ellipse(aes(x = ggdata[,1], y = ggdata[,2], fill = colorBy), 
                          geom = "polygon", type = "norm", level = 0.5, alpha = 0.1, 
                          show.legend = TRUE) + guides(color = guide_legend(phenotype_str), 
                                                       fill = guide_legend(phenotype_str))
  }
  else {
    if (is.numeric(colorBy)) {
      p <- p + guides(color = guide_colorbar(phenotype_str))
    }
    else if (!identical(colorBy, "")) {
      p <- p + guides(color = guide_legend(phenotype_str))
    }
  }
  return(p)
}






# # I don't like this way of combining
# #
# #' @describeIn plot_pca_v1
# #' @importFrom made4 ord plotarrays
# # ' @importFrom ade4 dudi.coa dudi.pca
# # ' @export plot_pca_v4
# #' @examples 
# #'
# #' plot_pca_v4(msnset, type='pca', phenotype="subject.type")
# #' plot_pca_v4(msnset, type='coa', phenotype="subject.type")
# 
# plot_pca_v4 <- function(eset, phenotype=NULL, ...){
#     phenotype <- as.factor(pData(eset)[[phenotype]])
#     ord.res <- ord(exprs(eset), classvec=phenotype, ...)
#     plotarrays(ord.res, ...)
# }
# 


