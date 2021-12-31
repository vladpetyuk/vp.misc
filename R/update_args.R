#' @title Update a List of Arguments
#'
#' @description Combine two lists. Used to update arguments in
#' \code{\link{plot_volcano}}, \code{\link{plot_pca}}, and
#' \code{\link{complex_heatmap}}. All list elements must be named.
#'
#' @param old_list old list
#' @param new_list new list
#'
#' @returns a list
#'
#' @importFrom purrr map2
#'
#'
#' @export
#'
#' @examples
#'
#' old_list <- list(color = "red", mapping = list(x = letters[1:3], y = 1:2))
#' new_list <- list(color = "blue", mapping = list(x = letters[4:6]))
#'
#' update_args(old_list, new_list)
#' # returns list(color = "blue", mapping = list(x = letters[4:6], y = 1:2))
#'


update_args <- function(old_list, new_list) {

  # Names only in one list
  old_only <- old_list[setdiff(names(old_list), names(new_list))]
  new_only <- new_list[setdiff(names(new_list), names(old_list))]

  # Names in common
  common_names <- intersect(names(old_list), names(new_list))

  # Subset args to common names
  if (!is.null(common_names)) {
    old_list_sub <- old_list[common_names]
    new_list_sub <- new_list[common_names]

    # Combine lists with common elements
    common_args <- map2(.x = old_list_sub, .y = new_list_sub,
                        .f = function(x, y) {
                          # If x or y is NULL, just return y
                          if (is.null(y) | is.null(x)) {
                            y
                          } else {
                            # Combine lists. y must be first!
                            c(y, x)
                          }
                        })

    common_args <- lapply(common_args, function(m) {
      if (is.list(m)) {
        # If m is a list, get the positions of the
        # first occurrences of each name and subset m
        y <- names(m)
        m[match(unique(y), y)]
      } else {
        # If m is a vector, select first element (from new_list_sub)
        m[1]
      }
    })

  } else {
    # Nothing in common
    common_args <- NULL
  }

  # Updated args
  return(c(common_args, old_only, new_only))
}

