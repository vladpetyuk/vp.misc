
#' Log2 transform and center around zero
#'
#' Log2 transform applied to expression crosstab of ExpressionSet or MSnSet
#' object followed by centering the median around zero. Converts one
#' MSnSet to another MSnSet.
#'
#' @param m ExpresionSet or MSnSet object
#'
#' @importFrom Biobase exprs<- exprs
#'
#' @return (MSnSet) MSnSet object
#'
#' @export log2_zero_center
#'
#' @examples
#' data(srm_msnset)
#' msnset2 <- log2_zero_center(msnset)


log2_zero_center <- function(m){

  exprs(m) <- log2(exprs(m))

  if (any(is.infinite(exprs(m))) == TRUE) {
    stop("After transformation, infinite values are present.")
  }

  exprs(m) <- sweep(exprs(m),
                    MARGIN = 1,
                    STATS = apply(exprs(m), 1, median, na.rm = TRUE),
                    FUN = "-")

  return(m)
}
