#' Reading MaxQuant Output at Peptide MS/MS Counts Level
#'
#' Reads in "evidence.txt" output (or its compressed version) as MSnSet.
#'
#'
#' @note The "evidence.txt" file can be compressed (gzip) to save space.
#' @param path character path to the folder containing "peptides.txt file".
#'          It should be in "{raw files folder}/combined/txt" of MaxQuant output.
#'
#' @return \code{MSnSet} object
#'
#' @importFrom MSnbase MSnSet
#' @importFrom utils read.delim
#'
#' @export readMaxQuantPeptideCounts
#'
#' @examples
#'
#' # to be added


readMaxQuantPeptideCounts <- function (path="."){

  fpath <- list.files(path = path, pattern = "evidence.txt", full.names = TRUE)
  x <- read.delim(fpath, check.names = FALSE, stringsAsFactors = FALSE)

  x_counts <- x %>%
    group_by(Sequence, `Leading razor protein`, Experiment) %>%
    summarise(spectral_count = sum(`MS/MS count`, na.rm=TRUE)) %>%
    rename(protein = `Leading razor protein`) %>%
    filter(!grepl("^CON", protein)) %>%
    filter(!grepl("^Contaminant", protein)) %>%
    filter(!grepl("^REV", protein))

  # let's make MSnSet out of it
  x_data <- pivot_wider(x_counts,
                        values_from = "spectral_count",
                        values_fill = 0,
                        names_from = "Experiment")
  x_data_m <- as.matrix(x_data[,-c(1,2)])
  rownames(x_data_m) <- x_data$Sequence

  f_data <- x_data[,c(1,2)] %>%
    as.data.frame() %>%
    `rownames<-`(.$Sequence)

  p_data <- data.frame(Experiment = colnames(x_data_m),
                       row.names = colnames(x_data_m))

  m <- MSnSet(exprs = x_data_m, fData = f_data, pData = p_data)

  return(m)

}


utils::globalVariables(c("Sequence", "Leading razor protein", "Experiment",
                         "MS/MS count", "protein"))

