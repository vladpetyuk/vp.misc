#' Construct Pairwise Comparisons
#'
#' Construct a vector of comparisons of the form \code{namea-nameb},
#' where \code{a} and \code{b} are unique elements of the factor or character
#' vector \code{x}, and \code{name} is a character string.
#'
#' @param x character vector or factor.
#' @param name character appended before each element of x when constructing
#'     the contrasts. Default is \code{""}.
#' @param ref character; value of x that will be used to compare all other
#'     values against. Creates reference comparisons instead of all pairwise
#'     comparisons.
#' @param sep character used to separate terms in the comparisons. The default,
#'     \code{"-"}, generates comparisons of the form \code{"a-b"}.
#'
#' @details
#' \code{paircomp} is primarily intended to generate contrasts for
#' \code{\link[MSnSet.utils]{limma_contrasts}}. If using \code{paircomp} for
#' this purpose, \code{x} will be a column in \code{pData(eset)} and \code{name}
#' will be the name of that column.
#'
#' @return character vector of comparisons
#'
#' @export paircomp
#'
#' @examples
#' library(MSnSet.utils)
#'
#' # Character vector for comparisons
#' x <- LETTERS[1:5]
#'
#' # All pairwise comparisons
#' paircomp(x)
#'
#' # Compare against specific level
#' paircomp(x, ref = "B")
#'
#' # Compare against specific level and append name
#' paircomp(x, name = "test", ref = "B")


paircomp <- function(x, name = "", ref = "", sep = "-") {

  # Check x
  if (!(is.character(x) | is.factor(x))) {
    stop("x must be a factor or character vector")
  }

  ## Prepare x
  if (is.character(x)) {
    # If x is a character vector, use the unique non-missing values
    x <- unique(x[!is.na(x)])
  } else {
    # If x is a factor, use the levels
    x <- levels(x)

    if (ref == "") {
      # If ref not provided, set to first level
      ref <- x[1]
    } else {
      # Check that provided ref is valid
      if (!(ref %in% x)) {
        stop("ref is not in levels(x)")
      }
    }
  }

  ## Create pairwise comparisons
  ctr <- c() # Vector to store comparisons

  if (ref == "") {
    # If no reference - all pairwise comparisons
    for (i in 1:(length(x) - 1)) {
      a <- x[i]
      b <- x[(i + 1):length(x)]

      ctr <- c(ctr, paste0(name, a, sep, name, b))
    }
  } else {
    # Compare to reference
    for (i in 1:length(x)) {
      a <- x[i]
      b <- ref

      ctr_i <- ifelse(a != b, paste0(name, a, sep, name, b), NA)
      ctr <- c(ctr, ctr_i)
    }
    ctr <- ctr[!is.na(ctr)] # remove NA comparison
  }

  return(ctr)
}


## Old code:
# pairwise_contrasts <- function(x, fct, ref.str = NULL,
#                                make_names = TRUE, sep = "-") {
#
#   x[, fct] <- as.factor(x[, fct])
#
#   if(!is.null(ref.str)) {
#     # Error if ref.str is not in the levels of fct
#     if(any(!(ref.str %in% levels(x[, fct])))) {
#       wrong_ref <- ref.str[which(!(ref.str %in% levels(x[, fct])))]
#       stop(paste0("ref.str = c(\"",
#                   paste(wrong_ref,
#                         collapse = "\", \""), "\")",
#                   ifelse(length(wrong_ref) > 1, " are ", " is "),
#                   "not in the levels of fct."))
#     } else {
#       grid_mat <- levels(x[, fct]) %>%
#         expand.grid(., .) %>%
#         # Filter to remove cases where
#         # 1) the difference between the means is 0
#         # 2) the first coefficient is a reference
#         # 3) The second coefficient is not a reference
#         filter(Var1 != Var2,
#                !(Var1 %in% ref.str),
#                Var2 %in% ref.str)
#     }
#   } else {
#     # grid_mat <- combinations(n = nlevels(x[, fct]),
#     #                          r = 2, # Pairwise contrasts
#     #                          v = levels(x[, fct]))
#
#     for (i in 1 : (length(y) - 1)) {
#       a <- y[i]
#       b <- y[(i + 1) : length(y)]
#       res <- c(res, sprintf("%s%s-%s%s", fct, a, fct, b))
#     }
#   }
#
#   # make_names = TRUE appends the factor name to the
#   # beginning of each level of fct.
#   if (make_names) {
#     contrasts <- paste0(fct, grid_mat[, 1], sep, fct, grid_mat[, 2])
#   } else {
#     contrasts <- paste0(grid_mat[, 1], sep, grid_mat[, 2])
#   }
#
#   return(contrasts)
# }
#
# utils::globalVariables(c(".", "Var1", "Var2"))

