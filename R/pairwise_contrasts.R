#' Construct Vector of Contrasts
#'
#' Construct a vector of contrasts of the form \code{a-b}, where \code{a} and
#' \code{b} are factor levels.
#'
#' @param x matrix or data frame that contains the column specified
#'     by \code{fct}.
#' @param fct character or factor variable from which contrasts will
#'     be generated. If \code{fct} is not a factor, it will be converted
#'     to one.
#' @param ref.str level(s) of \code{fct} that all other levels should
#'     be compared against. This is a character vector.
#' @param make_names if \code{TRUE} (the default), the name of \code{fct}
#'     will be appended to the beginning of each level of \code{fct} in
#'     the contrasts. This is required for
#'     \code{\link[MSnSet.utils]{limma_contrasts}}.
#' @param sep string used to separate terms in the contrasts. The default,
#'     \code{"-"}, generates contrasts of the form \code{"a-b"}.
#'
#' @return character vector of contrasts.
#'
#' @importFrom gtools combinations
#' @importFrom dplyr filter %>%
#'
#' @export pairwise_contrasts
#'
#' @examples
#' # Load package
#' library(MSnSet.utils)
#'
#' # Data frame for testing
#' x <- data.frame(v1 = c("a", "b", NA, "c", "d"),
#'                 v2 = c("e", "e", "f", "f", "g")) %>%
#' x$v3 <- ifelse(is.na(x$v1) | is.na(x$v2), NA,
#'                paste(x$v1, x$v2, sep = "_"))
#'
#' pairwise_contrasts(x, "v1", make_names = FALSE)
#' # This results in an error
#' pairwise_contrasts(x, "v1", ref.str = c("hello", "a"))
#' pairwise_contrasts(x, "v2")
#' pairwise_contrasts(x, "v2", ref.str = "g",
#'                    make_names = FALSE)
#' pairwise_contrasts(x, "v3", ref.str = c("b_e", "a_e"),
#'                    make_names = FALSE)


pairwise_contrasts <- function(x, fct, ref.str = NULL,
                               make_names = TRUE, sep = "-") {

  x[, fct] <- as.factor(x[, fct])

  if(!is.null(ref.str)) {
    # Error if ref.str is not in the levels of fct
    if(any(!(ref.str %in% levels(x[, fct])))) {
      wrong_ref <- ref.str[which(!(ref.str %in% levels(x[, fct])))]
      stop(paste0("ref.str = c(\"",
                  paste(wrong_ref,
                        collapse = "\", \""), "\")",
                  ifelse(length(wrong_ref) > 1, " are ", " is "),
                  "not in the levels of fct."))
    } else {
      grid_mat <- levels(x[, fct]) %>%
        expand.grid(., .) %>%
        # Filter to remove cases where
        # 1) the difference between the means is 0
        # 2) the first coefficient is a reference
        # 3) The second coefficient is not a reference
        filter(Var1 != Var2,
               !(Var1 %in% ref.str),
               Var2 %in% ref.str)
    }
  } else {
    grid_mat <- combinations(n = nlevels(x[, fct]),
                             r = 2, # Pairwise contrasts
                             v = levels(x[, fct]))
  }

  # make_names = TRUE appends the factor name to the
  # beginning of each level of fct.
  if (make_names) {
    contrasts <- paste0(fct, grid_mat[, 1], sep, fct, grid_mat[, 2])
  } else {
    contrasts <- paste0(grid_mat[, 1], sep, grid_mat[, 2])
  }

  return(contrasts)
}

