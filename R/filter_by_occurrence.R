
#' Retain features passing occurrence threshold
#'
#' @param m ExpresionSet or MSnSet object
#' @param occurrence Minimum proportion of observations
#'
#' @importFrom Biobase exprs
#'
#' @return (MSnSet) MSnSet object
#'
#' @export filter_by_occurrence
#'
#' @examples #
#' data(srm_msnset)
#' MSnSet.utils::image_msnset(msnset)
#' msnset2 <- filter_by_occurrence(msnset, occurrence = 0.6)
#' MSnSet.utils::image_msnset(msnset2)


filter_by_occurrence <- function(m, occurrence = 0.5){

  m[rowSums(!is.na(exprs(m))) >= ncol(m) * occurrence,]

}
