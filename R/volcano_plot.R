#' Volcano Plot
#' 
#' A convenience function for creating volcano plot.
#' 
#' @param logFC a numeric vector of log2 of fold change for each feature
#' @param significance a numeric vector of significance values 
#'          (p-value, q-value or adjusted p-value)
#' @param feature_names a character vector with the names of the features to
#'          diplay or NULL (default)
#' @param threshold a numeric value where to draw a 
#'          significance threshold or NULL (default)
#' @param top_n_names number of top significant features 
#'          to show the names (default 5)
#' @param rep.fact see \code{\link[FField]{FFieldPtRep}} for details. 
#'          Default here is 500.
#' @param adj.lmt see \code{\link[FField]{FFieldPtRep}} for details. 
#'          Default here is 3.
#' @param adj.max see \code{\link[FField]{FFieldPtRep}} for details. 
#'          Default here is 40.
#' @return plot
#' @importFrom FField FFieldPtRep
#' @importFrom scales trans_new log_breaks pretty_breaks
#' @importFrom ggplot2 ggplot geom_point aes scale_y_continuous theme
#'              theme_bw xlim ylim geom_hline geom_segment geom_text
#' @import qvalue
#' @export volcano_plot
#' 
#' @examples
#' library("vp.misc")
#' data("cptac_oca")
#' library("limma")
#' ee <- oca.set
#' ee <- ee[rowSums(!is.na(exprs(ee))) >= 30,]
#' model.string <- "~ PLATINUM.STATUS + AGE" # model with covariates
#' coef.string <- "PLATINUM.STATUS" # note, only two levels allowed here
#' 
#' res <- limma_a_b(oca.set, 
#'                  model.str = "~ PLATINUM.STATUS + AGE", 
#'                  coef.str = "PLATINUM.STATUS")
#' head(res)
#' 
#' #library(qvalue)
#' library(ggplot2)
#' # res <- subset(res, !is.na(res$P.Value))
#' #res$q.value <- qvalue(res$P.Value)$qvalue
#' volcano_plot(res$logFC, 
#'              res$P.Value, 
#'              feature_names = rownames(res), 
#'              threshold = 0.05, 
#'              top_n_names = sum(res$P.Value < 0.001, na.rm=TRUE), 
#'              rep.fact=10000, adj.lmt=1000, adj.max=1000) + 
#'    theme(text=element_text(size = 20)) +
#'    ylab("p-value") +
#'    xlab("log2 fold change")

volcano_plot <- function(logFC, 
                         significance, 
                         feature_names=NULL, 
                         threshold = NULL, 
                         top_n_names=5, 
                         rep.fact=500, adj.lmt=3, adj.max=40){
    
    # for y scale transform
    log10_rev_trans <- trans_new(
        "log10_rev",
        function(x) -log10(x),
        function(x) 10^(-x),
        breaks = function(x) {
            y <- log_breaks(10)(x)
            rev(y)},
        domain = c(1e-100, Inf)
    )
    
    # range on confidence
    # will step like 1, 0.3, 0.1, 0.03, 0.01, ...
    min_ = min(significance, na.rm = T)
    breaks <-signif(10^(pretty_breaks()(0:floor(log10(min_)*2))/2),1)
    
    res <- data.frame(logFC, significance)
    p <- ggplot(res, aes(x=logFC, y=significance)) +
        geom_point() + 
        scale_y_continuous(trans=log10_rev_trans, 
                           breaks=breaks, 
                           limits=c(1,min(breaks))) +
        xlim(-max(range(logFC)),+max(range(logFC))) +
        theme_bw()
    if(!is.null(threshold)) 
        p <- p + geom_hline(yintercept=threshold, col='red', linetype='dashed')
    
    # showing names of top
    scale_to <- function(x){
        x <- x - min(x, na.rm=TRUE)
        x <- x/max(x, na.rm=TRUE)
        return(x)
    }
    scale_from <- function(x.t, x.o){
        x.t <- x.t * (max(x.o, na.rm=T) - min(x.o, na.rm=T))
        x.t <- x.t + min(x.o, na.rm=T)
        return(x.t)
    }
    if(!is.null(feature_names) & 
       length(top_n_names) == 1 &
       top_n_names > 0){
        i <- order(significance)[seq_len(top_n_names)]
        res_names <- data.frame(x=logFC[i], 
                                y=significance[i], 
                                lbl=as.character(feature_names[i]))
        # fixing crowding with FField::FFieldPtRep
        xt <- scale_to(res_names$x)*100
        yt <- scale_to(-log10(res_names$y))*100
        coords <- FFieldPtRep(cbind(xt,yt), 
                              rep.fact=rep.fact, 
                              adj.lmt=adj.lmt, 
                              adj.max=adj.max)/100
        res_names$xff <- scale_from(coords$x, res_names$x)
        res_names$yff <- 10^(-scale_from(coords$y, -log10(res_names$y)))
        # end of fixing crowding
        p + geom_point(mapping=aes(x=x, y=y), data=res_names, color='red') +
            geom_text(mapping=aes(x=xff, y=yff, label=lbl), data=res_names) +
            geom_segment(mapping=aes(x=x, y=y, xend=xff, yend=yff), 
                         data=res_names, color='grey')
        
    }
    
}



