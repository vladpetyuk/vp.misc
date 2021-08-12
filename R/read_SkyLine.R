#' Read SkyLine Results
#'
#' Reads SkyLine ouput file "Peptide Ratio Results.csv" as MSnSet
#'
#' @param file character, path to the file
#'
#' @return \code{\link[MSnbase]{MSnSet}} object
#'
#' @importFrom plyr ddply
#' @importFrom MSnbase MSnSet
#' @importFrom utils read.csv
#' @importFrom reshape2 acast
#'
#' @export readSkyLinePRR
#'
#' @examples
#' \dontrun{#to be filled}


readSkyLinePRR <- function(file){
    x <- read.csv(file, na.strings = '#N/A', stringsAsFactors = FALSE)
    xf <- plyr::ddply(x, .variables=c("Peptide.Sequence", "Protein.Name"),
                plyr::summarise, mean.rt = mean(Peptide.Retention.Time, na.rm=T))
    rownames(xf) <- xf$Peptide.Sequence
    xe <- acast(x, Peptide.Sequence ~ Replicate.Name, value.var = "Ratio.To.Standard")
    xp <- data.frame(Replicate.Name=colnames(xe),
                     row.names=colnames(xe),
                     stringsAsFactors = FALSE)
    msnset <- MSnSet(exprs = xe, fData = xf, pData = xp)
    if(validObject(msnset))
        return(msnset)
    else
        return(NULL)
}

utils::globalVariables("Peptide.Retention.Time")

