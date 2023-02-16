
#' Visualize MSnSet
#'
#' Visualize MSnSet
#'
#' @param m MSnSet object
#' @param valueName name of the value to be displayed. Default - "value".
#' @param facetBy character of pheno data column containing factor value.
#'          Default is NULL.
#' @param sOrderBy character of pheno data telling how to order samples
#' @param fOrderBy character of feature data telling how to order features
#' @param valRange number for the pseudocolor range from -valRange to +valRange.
#'          Default is NULL. In that case the scale goes from 0.025 to 0.975
#'          quantile.
#' @param maxNRows maximum number of rows to display. Default is 50.
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_raster scale_fill_gradientn aes facet_grid
#' @importFrom ggplot2 theme element_blank element_rect xlab ylab labs
#' @importFrom ggplot2 element_text scale_fill_gradient2
#' @importFrom scales rescale
#' @importFrom gplots bluered
#' @importFrom grid unit
#' @importFrom stats quantile reorder
#'
#' @export image_msnset


image_msnset <- function(m, valueName="value",
                         facetBy=NULL,
                         sOrderBy=NULL,
                         fOrderBy=NULL,
                         valRange=NULL,
                         maxNRows=50){

    # convertion to long format
    mlong <- melt(exprs(m),
                  varnames=c("feature id", "sample name"),
                  value.name='value')
    mlong[['feature id']] <- as.character(mlong[['feature id']])
    mlong[['sample name']] <- as.character(mlong[['sample name']])
    mlong <- merge(mlong, fData(m), by.x="feature id", by.y=0) # slowish
    mlong <- merge(mlong, pData(m), by.x="sample name", by.y=0) # slow
    x <- mlong # not nice
    # order by feature name just for image purpose
    x[['feature id']] <-
        ordered(x[['feature id']],
                levels=rev(sort(unique(x[['feature id']]))))

    #----------------------------
    if(is.null(valRange))
        valRange <- mean(abs(quantile(x$value, c(0.025, 0.975), na.rm = TRUE)))
    x$value[x$value > +valRange] <- +valRange
    x$value[x$value < -valRange] <- -valRange
    qn01 <- rescale(c(c(-valRange,+valRange),range(x$value,na.rm=TRUE)))
    #----------------------------

    if(!is.null(sOrderBy))
        x[['sample name']] <- reorder(x[['sample name']], x[[sOrderBy]])

    if(!is.null(fOrderBy))
        x[['feature id']] <- reorder(x[['feature id']], x[[fOrderBy]])

    if(!is.null(facetBy)) x$facetBy <- x[[facetBy]]

    p <- ggplot(x, aes(x=`sample name`, y=`feature id`, fill=value)) +
        geom_raster() +
        scale_fill_gradientn(
            colours=bluered(100),
            values = c(0, seq(qn01[1], qn01[2], length.out = 98), 1),
            limits = c(-valRange,+valRange)) +
        # scale_fill_gradientn(colours=bluered(100)) +
        # scale_fill_gradient2(low="blue", high="red", na.value="black", name="")
        theme(
            axis.text.x=element_text(angle=+90),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(1.5, "lines"),
            panel.border = element_rect(linetype = "dashed",
                                        size=1, colour = "black", fill=NA),
            legend.key.height = unit(2, "lines")
        )
    if(nlevels(x$`feature id`) > maxNRows)
    {
        p <- p + theme(axis.text.y = element_blank())
        p <- p + theme(axis.ticks.y = element_blank())
    }
    if(!is.null(facetBy))
        p <- p + facet_grid( . ~ facetBy, scales='free', space='free')
    return(p)
}

utils::globalVariables(c("sample name", "feature id", "value"))



image_msnset_old <- function(m, valueName="value", facetBy=NULL,
                         sOrderBy=NULL, valRange=NULL){

    # convertion to long format
    mlong <- melt(exprs(m),
                  varnames=c("feature id", "sample name"),
                  value.name='value')
    mlong[['feature id']] <- as.character(mlong[['feature id']])
    mlong[['sample name']] <- as.character(mlong[['sample name']])
    mlong <- merge(mlong, fData(m), by.x="feature id", by.y=0) # slowish
    mlong <- merge(mlong, pData(m), by.x="sample name", by.y=0) # slow
    x <- mlong # not nice
    # order by feature name just for image purpose
    x[['feature id']] <-
        ordered(x[['feature id']],
                levels=rev(sort(unique(x[['feature id']]))))

    #----------------------------
    qn <- mean(abs(quantile(x$value, c(0.025, 0.975), na.rm = TRUE)))
    qn <- c(-qn, +qn)
    qn01 <- rescale(c(qn,range(x$value,na.rm=TRUE)))
    #----------------------------

    # what does that do? reordering by run order?
    #     x[['sample name']] <-
    #         with(x, reorder(`sample name`, SampleNum))
    if(!is.null(sOrderBy))
        x[['sample name']] <- reorder(x[['sample name']], x[[sOrderBy]])

    if(!is.null(facetBy)) x$facetBy <- x[[facetBy]]

    p <- ggplot(x, aes(x=`sample name`, y=`feature id`, fill=value)) +
        geom_raster() +
        scale_fill_gradientn(
            colours=bluered(100),
            values = c(0, seq(qn01[1], qn01[2], length.out = 98), 1)) +
        theme(
            axis.text.x=element_text(angle=+90),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(1.5, "lines"),
            panel.border = element_rect(linetype = "dashed",
                                        size=1, colour = "black", fill=NA),
            legend.key.height = unit(2, "lines")
        )
    if(nlevels(x$`feature id`) > 50) # todo: to be an argument
    {
        p <- p + theme(axis.text.y = element_blank())
        p <- p + theme(axis.ticks.y = element_blank())
    }
    if(!is.null(facetBy))
        p <- p + facet_grid( . ~ facetBy, scales='free', space='free')
    invisible(p)
}

utils::globalVariables(c("sample name", "feature id", "value"))


# further arguments to provide: x_orderBy, x_splitBy

# data(hndc)
# exprs(m) <- log2(exprs(m))
# exprs(m) <- sweep(exprs(m), 1, rowMeans(exprs(m)))
# library("reshape2")
