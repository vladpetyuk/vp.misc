

#' Reading MaxQuant Output
#' 
#' Reads "proteinGroups.txt" output as MSnSet. 
#' 
#' So far uses only LFQ data as quantitative data. As feature data, only first
#' 12 columns and iBAQ columns are used. Features are based on the genes.
#' In case of multiple proteins matching the same gene, only the one with 
#' the highest iBAQ value is returned.
#' @note The "proteinGroups.txt" file can be compressed to save space.
#' @param path character path to the folder containing "proteinGroups.txt file".
#'          It should be in "{raw files folder}/combined/txt"
#' @return \code{MSnSet} object
#' @importFrom Biobase lcPrefix
#' @importFrom MSnbase MSnSet
#' @importFrom plyr ddply
#' @export readMaxQuantProtGroups
#' @examples
#' m <- readMaxQuantProtGroups(system.file("extdata/MaxQuant",
#'                                          package="vp.misc"))
#' library("vp.misc")
#' exprs(m) <- log2(exprs(m))
#' exprs(m) <- sweep(exprs(m), 1, rowMeans(exprs(m), na.rm=TRUE), '-')
#' image_msnset(m)
#' 
readMaxQuantProtGroups <- function(path){
    # no options in the current version
    # use genes for IDs
    # and LFQ for quant data
    
    # fpath <- file.path(path, "proteinGroups.txt")
    # I assume the file can be compressed for the sake of space
    # thus ther is a bit relaxed search for file
    fpath <- 
        list.files(path = path, 
                   pattern = "proteinGroups.txt",
                   full.names = TRUE)
    stopifnot(length(fpath) == 1)
    
    x <- read.delim(fpath, check.names = FALSE, stringsAsFactors = FALSE)
    
    # feature IDs in first 12 columns
    id.cols <- c("Protein IDs", 
                 "Majority protein IDs", 
                 "Peptide counts (all)",
                 "Peptide counts (razor+unique)", 
                 "Peptide counts (unique)", 
                 "Protein names", 
                 "Gene names", 
                 "Fasta headers", 
                 "Number of proteins", 
                 "Peptides", 
                 "Razor + unique peptides",
                 "Unique peptides")
    ibac.col <- "iBAQ" # to resolve gene ambiguity
    quant.cols <- grep("^LFQ", colnames(x), value = TRUE)
    x <- x[,c(id.cols, ibac.col, quant.cols)]
    
    # trim the quant.cols name
    pref <- Biobase::lcPrefix(quant.cols)
    colnames(x) <- sub(pref, '', colnames(x))
    quant.cols <- sub(pref, '', quant.cols)
    
    #
    #.. GENE LEVEL
    # Let's get rid of "CON" "REV" and empty gene names, then check for redundancy.
    not.con <- !grepl('CON__', x$`Majority protein IDs`)
    not.rev <- !grepl('REV__', x$`Majority protein IDs`)
    not.empty <- x$`Gene names` != ''
    x <- x[not.con & not.rev & not.empty,]
    gns <- sapply(strsplit(x$`Gene names`, split = ';'), '[', 1)
    x$feature.name <- gns
    # retain genes with higher iBAQ
    x <- plyr::ddply(.data = x, .variables = ~ feature.name,
                     .fun = function(d){d[which.max(d$iBAQ),]})
    
    # to MSnSet
    x.exprs <- as.matrix(x[,quant.cols])
    x.exprs[x.exprs == 0] <- NA
    rownames(x.exprs) <- x$feature.name
    x.pdata <- data.frame(sample.name = colnames(x.exprs), 
                          stringsAsFactors = FALSE)
    rownames(x.pdata) <- colnames(x.exprs)
    x.fdata <- x[,c(id.cols, ibac.col)]
    rownames(x.fdata) <- x$feature.name
    ans <- MSnbase::MSnSet(exprs = x.exprs, 
                           fData = x.fdata,
                           pData = x.pdata)
    if (validObject(ans)) 
        return(ans)
}


