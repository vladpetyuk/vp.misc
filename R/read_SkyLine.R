#' Read SkyLine Results
#'
#' Reads SkyLine output file "Peptide Ratio Results.csv" as MSnSet
#'
#' @param file character, path to the file
#'
#' @return \code{\link[MSnbase]{MSnSet}} object
#'
#' @importFrom dplyr %>% group_by summarise
##' @importFrom tidyr pivot_wider
#' @importFrom MSnbase MSnSet
#' @importFrom utils read.csv
#' @importFrom reshape2 acast
#'
#' @export readSkyLinePRR
#'
#' @examples
#' \dontrun{#to be filled}


readSkyLinePRR <- function(file) {
    x <- read.csv(file, na.strings = '#N/A', stringsAsFactors = FALSE)
    # xf <- plyr::ddply(x, .variables=c("Peptide.Sequence", "Protein.Name"),
    #             plyr::summarise, mean.rt = mean(Peptide.Retention.Time, na.rm=T))
    xf <- x %>%
        group_by(Peptide.Sequence, Protein.Name) %>%
        summarise(mean.rt = mean(Peptide.Retention.Time, na.rm = TRUE)) %>%
        as.data.frame()

    rownames(xf) <- xf$Peptide.Sequence
    # This is one of two times the reshape2 package is used.
    # Instead, dplyr version below should work.
    xe <- acast(x, Peptide.Sequence ~ Replicate.Name, value.var = "Ratio.To.Standard")
    # xe <- x %>%
    #     pivot_wider(id_cols = c("Peptide.Sequence", "Replicate.Name"),
    #                 values_from = "Ratio.To.Standard")

    xp <- data.frame(Replicate.Name=colnames(xe),
                     row.names=colnames(xe),
                     stringsAsFactors = FALSE)
    msnset <- MSnSet(exprs = xe, fData = xf, pData = xp)
    if(validObject(msnset))
        return(msnset)
    else
        return(NULL)
}

utils::globalVariables(c("Peptide.Sequence", "Protein.Name",
                         "Peptide.Retention.Time"))

