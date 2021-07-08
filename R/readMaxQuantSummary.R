

#' Reading MaxQuant Output
#'
#' Reads in "summary.txt" output.
#'
#'@param path character path to the folder containing "summary.txt file".
#'          It should be in "{raw files folder}/combined/txt"
#'
#' @return data.frame
#' @export readMaxQuantSummary
#' @examples
#'
#' path <- system.file("extdata/MaxQuant", package="MSnSet.utils")
#' smmr <- readMaxQuantSummary(path)


readMaxQuantSummary <- function (path) {
    fpath <- list.files(path = path, pattern = "summary.txt",
                        full.names = TRUE)
    stopifnot(length(fpath) == 1)
    smmr <- read.delim(fpath, check.names = FALSE, stringsAsFactors = FALSE)

    smmr <- smmr[seq_len(nrow(smmr) - 1), ] # drop last row
    if (any(smmr$Experiment == "")) {
        warning("Dropping rows from \"summary.txt\" with missing experiment label. Please double check.")
        smmr <- smmr[which(smmr$Experiment != ""), ]
    }
    warning(paste("Found", nrow(smmr), "datasets"))
    smmr <- data.frame(dataset.name = smmr[, "Raw file"],
                       row.names = smmr[, "Experiment"],
                       stringsAsFactors = FALSE)
    return(smmr)
}