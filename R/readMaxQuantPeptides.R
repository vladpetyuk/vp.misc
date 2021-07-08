

#' Reading MaxQuant Output
#'
#' Reads in "peptides.txt" output (or its compressed version) as MSnSet.
#'
#'
#' @note The "peptides.txt" file can be compressed (gzip) to save space.
#' @param path character path to the folder containing "peptides.txt file".
#'          It should be in "{raw files folder}/combined/txt" of MaxQuant output.
#' @param msms_id_only logical If TRUE retains values associated with MS/MS peptide
#'          identifications only. Otherwise reports IDs identified both by MS/MS
#'          and match between the runs (MBR). Default it FALSE - that is reporting
#'          data from both type of identifications.
#' @param verbose numeric controls the text output
#'
#' @note The quantitative data reported by this function is in "Intensity" columns of
#'       the "peptides.txt" file.
#'
#' @return \code{MSnSet} object
#' @importFrom MSnbase MSnSet
#' @importFrom plyr ddply
#' @export readMaxQuantPeptides
#' @examples
#'
#' # label-free data
#' m <- readMaxQuantPeptides(system.file("extdata/MaxQuantPep",
#'                                        package="MSnSet.utils"))
#' exprs(m) <- log2(exprs(m))
#' exprs(m) <- sweep(exprs(m), 1, rowMeans(exprs(m), na.rm=TRUE), '-')
#' image_msnset(m)
#'
readMaxQuantPeptides <- function (path=".",
                                  msms_id_only = FALSE,
                                  verbose = 1){

    smmr <- readMaxQuantSummary(path)

    fpath <- list.files(path = path, pattern = "peptides.txt", full.names = TRUE)
    stopifnot(length(fpath) == 1)
    x <- read.delim(fpath, check.names = FALSE, stringsAsFactors = FALSE)

    if (verbose > 0) {
        print("MaxQuant columns")
        print(colnames(x))
    }
    id.cols <- c("Sequence", "Amino acid before", "Amino acid after",
                 "Proteins", "Leading razor protein",
                 "Unique (Groups)", "Unique (Proteins)",
                 "id", "Protein group IDs", "MS/MS Count")

    # in peptides it is only intensity
    quant.cols <- paste("Intensity", rownames(smmr), sep = " ")
    stopifnot(all(quant.cols %in% colnames(x)))
    cols_missing <- setdiff(id.cols, colnames(x))

    if (length(cols_missing) > 0) {
        warn_msg <- paste(cols_missing, collapse = ", ")
        warn_msg <- sprintf("Missing columns: %s.", warn_msg)
        warning(warn_msg)
        id.cols <- intersect(colnames(x), id.cols)
    }

    pref <- paste("Intensity", "\\s+", sep = "")
    colnames(x) <- sub(pref, "", colnames(x))
    quant.cols <- sub(pref, "", quant.cols)

    # removing proteins with reversed sequences
    not.rev <- x$Reverse != "+"
    x <- x[not.rev, ]

    # creating MSnSet object
    x.exprs <- as.matrix(x[,quant.cols])
    x.exprs[x.exprs == 0] <- NA
    #
    if(msms_id_only){
        it <- x[,grep("Identification type", colnames(x), value = T)]
        x.exprs[it != "By MS/MS"] <- NA
    }

    x$feature.name <- with(x, paste(`Leading razor protein`,`Sequence`,sep="@"))
    rownames(x.exprs) <- x$feature.name
    x.pdata <- data.frame(sample.name = colnames(x.exprs),
                          dataset.name = smmr[colnames(x.exprs),],
                          stringsAsFactors = FALSE)
    rownames(x.pdata) <- colnames(x.exprs)

    x.fdata <- x[,id.cols]

    rownames(x.fdata) <- x$feature.name
    ans <- MSnbase::MSnSet(exprs = x.exprs, fData = x.fdata,
                           pData = x.pdata)
    return(ans)
}
