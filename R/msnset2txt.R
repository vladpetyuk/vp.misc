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
#' @param drop currently either "features" or "samples"
#'
#' @importFrom utils write.table
#'
#' @export msnset2txt
#'
#' @examples
#' \dontrun{
#' data(srm_msnset)
#' msnset2txt(msnset)
#' }


msnset2txt <- function(m, prefix=NULL, sig.dig=3, drop=NULL){
    if(is.null(prefix))
        prefix <- deparse(substitute(m))
    #
    out.exprs <- signif(exprs(m), digits = sig.dig)
    out.features <- fData(m)
    out.pheno <- pData(m)
    #
    if(is.null(drop)){
        write.table(out.exprs,
                    file = paste(prefix, "expression.txt", sep="_"),
                    quote = FALSE, sep = '\t', row.names = TRUE, na='',
                    col.names = NA)
        write.table(out.features,
                    file = paste(prefix, "features.txt", sep="_"),
                    quote = FALSE, sep = '\t', row.names = TRUE, na='',
                    col.names = NA)
        write.table(out.pheno,
                    file = paste(prefix, "samples.txt", sep="_"),
                    quote = FALSE, sep = '\t', row.names = TRUE, na='',
                    col.names = NA)
    }else if(drop == "samples"){
        write.table(cbind(out.features, out.exprs),
                    file = paste(prefix, ".txt", sep=""),
                    quote = FALSE, sep = '\t', row.names = FALSE, na='')
    }else if(drop == "features"){
        write.table(cbind(out.pheno, t(out.exprs)),
                    file = paste(prefix, ".txt", sep=""),
                    quote = FALSE, sep = '\t', row.names = FALSE, na='')
    }else{
        message("invalid drop argument")
    }
    #
    invisible(NULL)
}

