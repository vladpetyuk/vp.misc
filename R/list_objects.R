#' Improved list of objects
#'
#' Create a \code{data.frame} containing information about objects in the global
#' environment.
#'
#' @return
#' \code{data.frame}
#' \describe{
#'   \item{\code{Type}}{the object class}
#'   \item{\code{Size}}{the object size (in bytes)}
#'   \item{\code{Rows}}{the number of rows, if applicable}
#'   \item{\code{Columns}}{the number of columns, if applicable}
#' }
#'
#' @param n maximum number of rows to output. Default is 10.
#' @param ... arguments passed to \code{.ls.objects}.
#'
#' @export lsos
#'


lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}


.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
        fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}

utils::globalVariables(c("object.size"))

