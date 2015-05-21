


#' Heatmap
#' 
#' Customized Heatmap. TO FILL WHAT EXACTLY IS CUSTOM ABOUT IT.
#' 
#' @param x MSnSet object
#' @param dist.fun distance function "eucledian" or"pearson"
#' @param col.ramp color mapping function. default gplots::bluered
#' @param z.transf logical perform Z-transform or not.
#' @param breaks color key breaks
#' @param linkage see ?hclust
#' @param nLabLim max limit of the row/column labels to show. default 30
#' @param ... further arguments to \code{gplots::heatmap.2}
#'      
#' @importFrom gplots heatmap.2 bluered
#' @export heatmap.3
#' 
#' @examples
#' data(srm_msnset)
#' set.seed(0)
#' clrz <- sample(colors(), 17)
#' heatmap.3(cor(exprs(msnset)), 
#'           dist.fun="pearson",
#'           linkage="average", 
#'           nLabLim=50,
#'           ColSideColors=c('red3','yellow3','green3')[as.factor(pData(msnset)$subject.type)],
#'           RowSideColors=clrz[as.factor(pData(msnset)$match.group)])

heatmap.3 <- function( x, 
                       dist.fun=c("eucledian","pearson"),
                       col.ramp=bluered,
                       # column.factors=NULL,
                       z.transf=c(FALSE, TRUE),
                       breaks=seq(-1,+1,length=100),
                       linkage=c("average", "ward.D", "ward.D2", "single", 
                                 "complete", "mcquitty", "median",
                                 "centroid"),
                       nLabLim=30,
                       ...)
{
    # library( "gplots" )
    # library( "RColorBrewer" )
    
    if(isTRUE(z.transf))
        x = x/apply(x,1,sd,na.rm=TRUE)
    
    # column coloring
    # if NULL, then ColSideColors has to be missing in the call
    #     if(!is.null(column.factors)){
    #         # match by column name first
    #         design = sapply( column.factors, grepl, colnames(x))
    #         # assign the condition name
    #         conditionsToColumns = apply( design, 1, function(xx){
    #             if(!any(xx)){
    #                 return(NA)
    #             }else{
    #                 return(column.factors[xx])
    #             }})
    #         cols = as.factor( conditionsToColumns )
    #         colScheme = brewer.pal( max(3,nlevels(cols)), "Set1")   
    #         ColSideColors = colScheme[cols]
    #     }else{
    #         ColSideColors = NULL
    #     }
    
    # selecting distance type
    dist.fun <- match.arg(dist.fun)
    if(dist.fun == "eucledian"){
        distfun=function(x, ...) dist(x, 
                                      method = "euclidean", 
                                      ...)
    }else if(dist.fun == "pearson"){
        distfun=function(x) as.dist((1-cor( t(x), 
                                            method="pearson",
                                            use="pairwise.complete.obs" ))/2)
    }
    
    # heatmap itself
    linkage <- match.arg(linkage)
    heatmap.expression = expression(
        heatmap.2(  as.matrix(x), 
                    trace="none", 
                    #                   dendrogram="row",
                    distfun=distfun,
                    hclustfun=function(xx,...)
                    {hclust(xx,method=linkage,...)},
                    col=col.ramp(length(breaks)-1),
                    labRow = if(nrow(x)>nLabLim) "" else NULL,
                    labCol = if(ncol(x)>nLabLim) "" else NULL,
                    symkey=T, 
                    breaks=breaks,
                    na.color="white",#gray(0.5),
                    ...
        ))
    # updating ColSideColors
    # browser()
    #     heatmap.expression[[1]]$ColSideColors <- list(...)$ColSideColors
    #     heatmap.expression[[1]]$RowSideColors <- list(...)$RowSideColors
    eval(heatmap.expression)
    
}

