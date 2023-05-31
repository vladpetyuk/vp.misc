
#' Filter rows for < 50% missingness
#'
#' Calculates for a given row the number of columns which contain non-NA values
#' and keep rows with > 50% non-NA. Converts one MSnSet to another MSnSet.
#'
#' @param m ExpresionSet or MSnSet object
#'
#' @importFrom Biobase exprs<-
#'
#' @export occurence_filter
#'
#' @examples #
#' data(srm_msnset)
#' MSnSet.utils::image_msnset(msnset)
#' msnset2 <- occurence_filter(msnset)
#' MSnSet.utils::image_msnset(msnset2)


occurence_filter <- function(m){

  m <- m[rowSums(!is.na(exprs(m))) > ncol(m) * 0.5,]

    return(m)
}
