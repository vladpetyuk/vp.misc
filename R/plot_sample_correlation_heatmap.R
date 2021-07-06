#' Heatmap for sample correlations
#'
#' Wrapper functions to generate correlation matrix plots
#' for MSnSet objects. There is a parameter for the phenotype to
#' order the columns (samples) by.
#' 
#' @param msnset (MSnSet)
#' @param phenotype (character) name of the column to order by (usually plex)
#' 
#' @param ... additional params to the heatmap call
#'
#' @return (list) plot object
#'
#' @importFrom MSnbase pData exprs 
#' @importFrom dplyr %>% select 
#' @importFrom pheatmap pheatmap
#'
#' @examples
#' 
#' library(MSnSet.utils)
#' data(cptac_oca)
#'
#' m <- oca.set
#' 
#' # Using base heatmap
#' plot_sample_correlation_heatmap(m, "Batch")
#' 
#' # Using `pheatmap` (pretty heatmap)
#' plot_sample_correlation_heatmap(m, "Batch")
#' 

#' @export
plot_sample_correlation_heatmap <- function(m, phenotype, ...) {
    m <- m[,order(pData(m)[,phenotype])]
    x <- cor(exprs(m), use = "complete.obs")
    x <- x[nrow(x):1, ]
    heatmap(x, Rowv = NA, Colv = NA, ...)
}

#' @export
plot_sample_correlation_pheatmap <- function(m, phenotype, ...) { 
    m <- m[,order(pData(m)[,phenotype])]
    x <- cor(exprs(m), use = "complete.obs")
    pheno = pData(m) %>%
        select(!!phenotype)
    pheatmap(x, 
             cluster_rows=F, cluster_cols=F,
             annotation_col = pheno, annotation_row = pheno)
}
