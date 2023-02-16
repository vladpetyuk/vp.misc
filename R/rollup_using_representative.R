
#' Selecting one feature to represent a group of features
#'
#' The typical use of the function would be selecting representative peptide
#' per protein. Another use is selecting representative FAIMS CV value per
#' proteoform. This function is alternative to reference rollup (rrollup).
#' Instead of scaling peptides' intensities, here we select one peptide to
#' represent a protein.
#'
#' @param msnset msnset (or most likely eset subclass) object
#' @param roll_by character. A column name in pData(msnset).
#'                           Controls grouping of the features.
#' @param prioritize_by character. Feature variables to prioritize the selection.
#'                          Note, ordering is always in the descending order.
#'
#' @return MSnSet object
#' @importFrom Biobase fData featureNames
#' @importFrom dplyr %>% group_by desc across all_of ungroup filter
#' @export rollup_using_representative
#'
#' @examples
#' # Not run
#' \dontrun{
#' fData(m)$proteoform <- with(fData(m), paste(Gene, pcGroup, sep="_"))
#' fData(m)$count <- rowSums(!is.na(exprs(m)))
#' m2 <- rollup_using_representative(m)
#' }

rollup_using_representative <- function(msnset,
                                        roll_by,
                                        prioritize_by){
  #
  sel_cols <- c(roll_by, prioritize_by)
  x <- fData(msnset)[,sel_cols]
  x$idx <- 1:nrow(x)
  y <- x %>%
    group_by(!!sym(roll_by)) %>%
    arrange(desc(across(all_of(prioritize_by))), .by_group = TRUE) %>%
    ungroup() %>%
    filter(!duplicated(!!sym(roll_by)))
  idx <- x$idx %in% y$idx
  msnset <- msnset[idx,]
  featureNames(msnset) <- fData(msnset)[,roll_by]
  return(msnset)
}

