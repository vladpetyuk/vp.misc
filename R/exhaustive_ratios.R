#' Exhaustive Ratios
#' 
#' For a group of features size N it exhaustively
#' generates N*(N-1)/2 ratios and returns another MSnSet.
#' 
#' There are few scenarios when this function can be useful.
#' \itemize{
#'      \item SRM data. Generate full exhaustive ratios to normalized per
#'          tissue leakage.
#'      \item Top-Down data. Generate ratios of proteoforms within a given gene.
#'      \item any data. Ratios within a given pathway to assess dysregulation.
#' }
#' 
#' @param m MSnSet object
#' @param INDEX character Name of one of the columns in fData(m). 
#'              If NULL (defaults), then ratios generated assuming all
#'              features (e.g. peptides or proteforms) belong to one
#'              grouping factor.
#' 
#' @export
exhaustive_ratios <- function(m, INDEX=NULL){
    x <- exprs(m)
    if(is.null(INDEX)){
        fctrs <- factor(rep(1, nrow(m)))
    }else{
        fctrs <- as.factor(fData(m)[[INDEX]])
        fctrs.to.retain <- names(which(table(fctrs) > 1))
        idx <- fctrs %in% fctrs.to.retain
        x <- x[idx,]
        fctrs <- droplevels(fctrs[idx])
    }
    #
    out <- lapply(levels(fctrs), function(fi){
        xi <- x[fctrs == fi,]
        out2 <- lapply(seq_len(nrow(xi)-1), function(i){
            name.i <- rownames(xi)[i]
            # if it has been log-trans, then '-'
            # so not really ratios in log scale
            xr <- sweep(xi[-(1:i),,drop=F], 2, xi[i,], '-')
            rownames(xr) <- paste(rownames(xi[-(1:i),,drop=F]), name.i, sep='-')
            return(xr)
        })
        out2 <- Reduce(rbind, out2)
    })
    out <- Reduce(rbind, out)
    # pData is the same as in m
    # fData is a bit tricky.  Howerver, it is parsed out from "out".
    Peptide.1 <- sapply(strsplit(rownames(out), '-'), '[[', 1)
    Peptide.2 <- sapply(strsplit(rownames(out), '-'), '[[', 2)
    fData <- data.frame(Peptide.1, Peptide.2, 
                        Peptide.ID=rownames(out), 
                        row.names=rownames(out), stringsAsFactors=FALSE)
    m <- MSnSet(out, fData, pData(m))
    return(m)
}


# library("LewyBodies.SN.TopDown")
# data(iq_human_db)
# (m)     
# head(fData(m))
# table(table(fData(m)$ProteinName))
# image.lewy.iq(m[fData(m)$ProteinName == 'SYUA_HUMAN'])
# 
# debugonce(exhaustive_ratios)
# m2 <- exhaustive_ratios(m, INDEX="ProteinName")
# m3 <- exhaustive_ratios(m)
# 
# x <- eset_lm(m2, 
#              "y ~ subject.type + match.group", 
#              "y ~ match.group")



