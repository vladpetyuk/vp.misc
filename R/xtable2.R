#' 
#' Customized xtable
#'
#' The main customization pertains to alternating
#' row shading. This one will work only in LaTeX-type documents.
#' See \code{\link[xtable]{xtable}} for details.
#'
#' @param x typically a \code{data.frame}. 
#' @return an object of class \code{"xtable"}
#' @importFrom xtable xtable
#' @seealso \code{\link[xtable]{xtable}}
#' @name xtable2
#' 
xtable2 <- function(x, col.lab.rot=FALSE, ...){
    rws <- seq(1, (nrow(x)-1), by = 2)
    col <- rep("\\rowcolor[gray]{0.95}", length(rws))
    add.to.row = list(pos = as.list(rws), command = col)
    include.colnames = TRUE
    if(col.lab.rot){
        add.to.row$pos <- c(0, add.to.row$pos)
        cols <- paste0("\\rothead{",colnames(x),"}", collapse = " & ")
        cols <- paste0("\\rothead{} & ", cols, "\\\\")
        add.to.row$command <- c(cols, add.to.row$command)
        include.colnames = FALSE
    }
    print(xtable(x),
          booktabs = TRUE,
          add.to.row = add.to.row,
          include.colnames = include.colnames,
          ...)
}