

#' Reading MaxQuant Output
#'
#' Reads in "proteinGroups.txt" output (or its compressed version) as MSnSet.
#'
#' So far uses only LFQ data as quantitative data. As feature data, only first
#' 12 columns and iBAQ columns are used. Features are based on the genes.
#' In case of multiple proteins matching the same gene, only the one with
#' the highest iBAQ value is returned.
#' @note The "proteinGroups.txt" file can be compressed to save space.
#' @param path character path to the folder containing "proteinGroups.txt file".
#'          It should be in "{raw files folder}/combined/txt"
#' @param quantType character Defines pattern what type of column to
#'                  use for quantificaiton. 
#'                  E.g. "LFQ intensity" or "Ratio H/L normalized".
#' @param verbose numeric controls the text output
#' 
#' @note It looks like the convention for column naming in MaxQuant is
#'       type of quantification followed by sample name separated by space.
#'       E.g. "LFQ intensity Sample1" or "Ratio H/L normalized Control33".
#'       In \code{quantType} you need to specify the full first component
#'       in the name that defines the type of quantification.
#'       
#' @note iBAQ option must be enabled in MaxQuant analysis because it is used
#'       to resolve ambiguity between two proteins matching one gene.
#'
#' @return \code{MSnSet} object
#' @importFrom MSnbase MSnSet
#' @importFrom plyr ddply
#' @export readMaxQuantProtGroups
#' @examples
#' 
#' # label-free data
#' m <- readMaxQuantProtGroups(system.file("extdata/MaxQuant",
#'                                          package="vp.misc"),
#'                             quantType="LFQ intensity")
#' exprs(m) <- log2(exprs(m))
#' exprs(m) <- sweep(exprs(m), 1, rowMeans(exprs(m), na.rm=TRUE), '-')
#' image_msnset(m)
#' 
#' # O18/O16 data
#' m <- readMaxQuantProtGroups(system.file("extdata/MaxQuant_O18",
#'                                          package="vp.misc"),
#'                             quantType="Ratio H/L normalized")
#' exprs(m) <- log2(exprs(m))
#' exprs(m) <- sweep(exprs(m), 1, rowMeans(exprs(m), na.rm=TRUE), '-')
#' image_msnset(m)
#' 
#' 
readMaxQuantProtGroups <- function(path, quantType, verbose=1){
    # no options in the current version
    # use genes for IDs
    
    #.. get dataset names
    # dataset names are in the summary.txt
    fpath <-
        list.files(path = path,
                   pattern = "summary.txt",
                   full.names = TRUE)
    stopifnot(length(fpath) == 1)
    smmr <- read.delim(fpath, check.names = FALSE, stringsAsFactors = FALSE)
    warning("Asumming \"summary.txt\" file contains 2n+1 rows (excluding header).
            Where n is the number of datasets.
            Please double check.")
    smmr <- smmr[seq_len((nrow(smmr)-1)/2),]
    smmr <- data.frame(dataset.name = smmr[,"Raw file"],
                       row.names = smmr[,"Experiment"],
                       stringsAsFactors = FALSE)
    # safety check to make sure there is nothing odd
    # about the summary file. Experiment names should mutually match.
    # stopifnot(all.equal(rownames(smmr), quant.cols))
    
    

    # fpath <- file.path(path, "proteinGroups.txt")
    # I assume the file can be compressed for the sake of space
    # thus there is a bit relaxed search for file
    fpath <-
        list.files(path = path,
                   pattern = "proteinGroups.txt",
                   full.names = TRUE)
    stopifnot(length(fpath) == 1)

    x <- read.delim(fpath, check.names = FALSE, stringsAsFactors = FALSE)
    if(verbose > 0){
        print("MaxQuant columns")
        print(colnames(x))
    }

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
    ibac.col <- "iBAQ" # to resolve gene ambiguity in situations when
                       # we need to select one gene per protein id.
                       # that ".+" is pretty much a hack to ensure that 
                       # I do not read in columns that have nothing to do 
                       # with samples
    quant.cols <- paste(quantType, rownames(smmr), sep=' ')
    # safety check
    stopifnot(all(quant.cols %in% colnames(x)))
    # quant.cols <- grep(paste("^",quantType,".+",sep=''), colnames(x), value = TRUE)
    x <- x[,c(id.cols, ibac.col, quant.cols)]

    # trim the quant.cols name and retain only sample names
    pref <- paste(quantType,"\\s+",sep="")
    colnames(x) <- sub(pref, '', colnames(x))
    quant.cols <- sub(pref, '', quant.cols)

    #
    #.. GENE LEVEL
    # Let's get rid of ("CON" not anymore) "REV" and empty gene names, 
    # then check for redundancy.
    # not.con <- !grepl('CON__', x$`Majority protein IDs`)
    not.rev <- !grepl('REV__', x$`Majority protein IDs`)
    not.empty <- x$`Gene names` != ''
    # x <- x[not.con & not.rev & not.empty,]
    x <- x[not.rev & not.empty,]
    gns <- sapply(strsplit(x$`Gene names`, split = ';'), '[', 1)
    x$feature.name <- gns
    # retain genes with higher iBAQ
    x <- plyr::ddply(.data = x, .variables = ~ feature.name,
                     .fun = function(d){d[which.max(d$iBAQ),]})

    #.. Denote potential contaminants
    contaminants <- grepl('CON__', x$`Majority protein IDs`)
    x$`Majority protein IDs` <- sub('CON__','',x$`Majority protein IDs`)
    x$isContaminant <- contaminants
    id.cols <- c(id.cols, 'isContaminant')

    # to MSnSet
    x.exprs <- as.matrix(x[,quant.cols])
    x.exprs[x.exprs == 0] <- NA
    rownames(x.exprs) <- x$feature.name
    x.pdata <- data.frame(sample.name = colnames(x.exprs),
                          dataset.name = smmr[colnames(x.exprs),],
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

