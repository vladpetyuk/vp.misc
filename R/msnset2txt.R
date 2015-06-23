#' Exporting MSnSet as text file
#' 
#' Writes MSnSet as three tables: expression, features and phenotypes.
#' The format is tab-delimited text.
#' 
#' @param m MSnSet object
#' @param prefix character. All three text file names will 
#'          start with this prefix. By default, the R object name.
#' @param sig.dig number of significant digits in expression data.
#'          Default is 3.
#' 
#' @export
msnset2text <- function(m, prefix=NULL, sig.dig=3){
    if(is.null(prefix)) 
        prefix <- deparse(substitute(m))
    # take care of expression
    out.exprs <- signif(exprs(m), digits = sig.dig)
    write.table(out.exprs, file = paste(prefix, "expression.txt", sep="_"),
                quote = FALSE, sep = '\t', row.names = TRUE)
    # take care of features
    out.features <- fData(m)
    write.table(out.features, file = paste(prefix, "features.txt", sep="_"),
                quote = FALSE, sep = '\t', row.names = TRUE)
    # take care of phenotypes
    out.pheno <- pData(m)
    write.table(out.pheno, file = paste(prefix, "samples.txt", sep="_"),
                quote = FALSE, sep = '\t', row.names = TRUE)
    invisible(NULL)
}

