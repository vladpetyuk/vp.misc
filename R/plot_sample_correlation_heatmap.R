#' Heatmap for sample correlations
#'
#' Wrapper functions to generate correlation matrix plots
#' for MSnSet objects. There is a parameter for the phenotype to
#' order the columns (samples) by.
#'
#' @param m an MSnSet object
#' @param phenotype (character) name of the column to order by (usually plex)
#'
#' @param ... additional arguments passed to \code{\link[stats]{heatmap}}
#'
#' @return (list) plot object
#'
#' @export plot_sample_correlation_heatmap
#'
#' @importFrom MSnbase pData exprs
#' @importFrom dplyr %>% select
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
#' plot_sample_correlation_pheatmap(m, "Batch")


plot_sample_correlation_heatmap <- function(m, phenotype, ...) {
    m <- m[,order(pData(m)[,phenotype])]
    x <- cor(exprs(m), use = "complete.obs")
    x <- x[nrow(x):1, ]
    heatmap(x, Rowv = NA, Colv = NA,
            scale = "none", symm = TRUE,
            ...)
}



#' @describeIn plot_sample_correlation_heatmap Using the \code{pheatmap} package
#' @importFrom pheatmap pheatmap
#' @export plot_sample_correlation_pheatmap
plot_sample_correlation_pheatmap <- function(m, phenotype, ...) {
    m <- m[,order(pData(m)[,phenotype])]
    x <- cor(exprs(m), use = "complete.obs")
    pheno = pData(m) %>%
        select(!!phenotype)
    pheatmap(x,
             cluster_rows = FALSE, cluster_cols = FALSE,
             annotation_col = pheno, annotation_row = pheno,
             ...)
}
