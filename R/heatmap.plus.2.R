#' Tweak of "heatmap.plus"
#'
#' This function is basically
#' \href{https://github.com/cran/heatmap.plus/blob/master/R/heatmap.plus.R}{heatmap.plus::heatmap.plus}
#' from \href{https://github.com/cran/heatmap.plus/blob/master}{heatmap.plus}
#' package, but it includes a pair of additional arguments controlling the
#' proportion of row and column color sidebars.
#'
#' @param x numeric matrix of the values to be plotted
#' @param Rowv determines if and how the *row* dendrogram should be computed and
#'   reordered. Either a dendrogram or a vector of values used to reorder the
#'   row dendrogram or \code{NA} to suppress any row dendrogram (and reordering)
#'   or by default, \code{NULL}.
#' @param Colv determines if and how the *column* dendrogram should be
#'   reordered. Has the same options as the Rowv argument above and additionally
#'   when \code{x} is a square matrix, \code{Colv = "Rowv"} means that columns
#'   should be treated identically to the rows.
#' @param distfun function used to compute the distance (dissimilarity) between
#'   both rows and columns. Defaults to \code{dist}.
#' @param hclustfun function used to compute the hierarchical clustering when
#'   \code{Rowv} or \code{Colv} are not dendrograms. Defaults to \code{hclust}.
#' @param reorderfun function(d,w) of dendrogram and weights for reordering the
#'   row and column dendrograms. The default uses \code{reorder.dendrogram}.
#' @param add.expr expression that will be evaluated after the call to
#'   \code{image}. Can be used to add components to the plot.
#' @param symm logical indicating if \code{x} should be treated symmetrically;
#'   can only be true when \code{x} is a square matrix.
#' @param revC logical indicating if the column order should be reversed for
#'   plotting, such that e.g., for the symmetric case, the symmetry axis is as
#'   usual.
#' @param scale character indicating if the values should be centered and scaled
#'   in either the row direction or the column direction, or none. The default
#'   is \code{"row"} if \code{symm} false, and \code{"none"} otherwise.
#' @param na.rm logical indicating whether \code{NA}'s should be removed.
#' @param margins numeric vector of length 2 containing the margins (see
#'   \code{par(mar= *)}) for column and row names, respectively.
#' @param ColSideColors (optional) character matrix with number of columns
#'   matching number of columns in \code{x}. Each column is plotted as a column
#'   similar to heatmap()'s ColSideColors. colnames() will be used for labels if
#'   present.
#' @param RowSideColors (optional) character matrix with number of rows matching
#'   number of rows in \code{x}. Each column is plotted as a row similar to
#'   heatmap()'s ColSideColors. colnames() will be used for labels if present.
#' @param propColSide numeric controls how much room to allocate for the side
#'   color bars. E.g. 0.5 is even split between the heatmap and the bars.
#'   Default is 0.1. That is, 10\% bars, 90\% heatmap.
#' @param propRowSide numeric see above
#' @param cexRow positive number, used as \code{cex.axis} in for the row axis
#'   labeling. The default currently only uses number of rows.
#' @param cexCol see above.
#' @param labRow character vectors with row labels to use; defaults to
#'   rownames().
#' @param labCol see above.
#' @param main main title; defaults to none.
#' @param xlab x-axis title; defaults to none.
#' @param ylab y-axis title; defaults to none.
#' @param keep.dendro logical indicating if the dendrogram(s) should be kept as
#'   part of the result (when \code{Rowv} and/or \code{Colv} are not \code{NA}).
#' @param verbose logical indicating if information should be printed.
#' @param ... additional arguments passed to \code{\link[graphics]{image}}.
#'
#'
#' @importFrom graphics image par layout axis mtext frame title
#' @importFrom stats reorder hclust as.dendrogram order.dendrogram sd
#'
#'
#' @export heatmap.plus.2
#'
#' @examples
#' z = matrix(rnorm(30),nrow=5,ncol=6);
#' rlab = matrix(as.character(c(1:5,2:6,3:7,4:8)),nrow=5,ncol=4);
#' clab = matrix(as.character(c(1:6,6:1)),nrow=6,ncol=2);
#' colnames(rlab) = LETTERS[1:dim(rlab)[2]];
#' colnames(clab) = 1:dim(clab)[2];
#' heatmap.plus.2(z, ColSideColors=clab,
#'                   RowSideColors=rlab,
#'                   propColSide=0.1,
#'                   propRowSide=0.2);


heatmap.plus.2 <- function(x,
                           Rowv = NULL,
                           Colv = if (symm) "Rowv" else NULL,
                           distfun = dist,
                           hclustfun = hclust,
                           reorderfun = function(d, w) reorder(d, w),
                           add.expr,
                           symm = FALSE,
                           revC = identical(Colv, "Rowv"),
                           scale = c("row", "column", "none"),
                           na.rm = TRUE,
                           margins = c(5, 5),
                           ColSideColors,
                           RowSideColors,
                           propColSide=0.1,
                           propRowSide=0.1,
                           cexRow = 0.2 + 1/log10(nr),
                           cexCol = 0.2 + 1/log10(nc),
                           labRow = NULL,
                           labCol = NULL,
                           main = NULL,
                           xlab = NULL,
                           ylab = NULL,
                           keep.dendro = FALSE,
                           verbose = getOption("verbose"),
                           ...)
{
    scale <- if (symm && missing(scale)) "none"
    else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("'x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("'margins' must be a numeric vector of length 2")
    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (is.null(Rowv))
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv))
        Colv <- colMeans(x, na.rm = na.rm)
    if (doRdend) {
        if (inherits(Rowv, "dendrogram"))
            ddr <- Rowv
        else {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            if (!is.logical(Rowv) || Rowv)
                ddr <- reorderfun(ddr, Rowv)
        }
        if (nr != length(rowInd <- order.dendrogram(ddr)))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1:nr
    if (doCdend) {
        if (inherits(Colv, "dendrogram"))
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr != nc)
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            ddc <- ddr
        }
        else {
            hcc <- hclustfun(distfun(if (symm) x else t(x)))
            ddc <- as.dendrogram(hcc)
            if (!is.logical(Colv) || Colv)
                ddc <- reorderfun(ddc, Colv)
        }
        if (nc != length(colInd <- order.dendrogram(ddc)))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1:nc
    x <- x[rowInd, colInd]
    labRow <- if (is.null(labRow))
        if (is.null(rownames(x)))
            (1:nr)[rowInd]
    else rownames(x)
    else labRow[rowInd]
    labCol <- if (is.null(labCol))
        if (is.null(colnames(x)))
            (1:nc)[colInd]
    else colnames(x)
    else labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    lmat <- rbind(c(NA, 3), 2:1)
    #     lwid <- c(if (doRdend) 1 else 0.05, 4)
    #     lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 4)
    lwid <- c(if (doRdend) 1 else 0.5, 4)
    lhei <- c((if (doCdend) 1 else 0.5) + if (!is.null(main)) 0.2 else 0, 4)
    if (!missing(ColSideColors)) {
        if (!is.matrix(ColSideColors))
            stop("'ColSideColors' must be a matrix")
        if (!is.character(ColSideColors) || dim(ColSideColors)[1] != nc)
            stop("'ColSideColors' dim()[2] must be of length ncol(x)")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        #lhei <- c(lhei[1], 0.2, lhei[2])
        lhei <- c(lhei[1], lhei[2]*propColSide, lhei[2]*(1-propColSide))
    }
    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors))
            stop("'RowSideColors' must be a matrix")
        if (!is.character(RowSideColors) || dim(RowSideColors)[1] != nr)
            stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)
        #lwid <- c(lwid[1], 0.2, lwid[2])
        lwid <- c(lwid[1], lwid[2]*propRowSide, lwid[2]*(1-propRowSide))
    }
    lmat[is.na(lmat)] <- 0
    if (verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei, "; lmat=\n")
        print(lmat)
    }
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        rsc = RowSideColors[rowInd, ]
        rsc.colors = matrix()
        rsc.names = names(table(rsc))
        rsc.i = 1
        for (rsc.name in rsc.names) {
            rsc.colors[rsc.i] = rsc.name
            rsc[rsc == rsc.name] = rsc.i
            rsc.i = rsc.i + 1
        }
        rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
        image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
        if (length(colnames(RowSideColors)) > 0) {
            axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1),
                 colnames(RowSideColors), las = 2, tick = FALSE)
        }
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        csc = ColSideColors[colInd, ]
        csc.colors = matrix()
        csc.names = names(table(csc))
        csc.i = 1
        for (csc.name in csc.names) {
            csc.colors[csc.i] = csc.name
            csc[csc == csc.name] = csc.i
            csc.i = csc.i + 1
        }
        csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
        image(csc, col = as.vector(csc.colors), axes = FALSE)
        if (length(colnames(ColSideColors)) > 0) {
            axis(2, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1),
                 colnames(ColSideColors), las = 2, tick = FALSE)
        }
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    if (!symm || scale != "none") {
        x <- t(x)
    }
    if (revC) {
        iy <- nr:1
        ddr <- rev(ddr)
        x <- x[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
              c(0, nr), axes = FALSE, xlab = "", ylab = "", ...)
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
         cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
         cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    par(mar = c(margins[1], 0, 0, 0))
    if (doRdend)
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else
        frame()
    par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2]))
    if (doCdend)
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    else if (!is.null(main))
        frame()
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    invisible(list(rowInd = rowInd,
                   colInd = colInd,
                   Rowv = if (keep.dendro && doRdend) ddr,
                   Colv = if (keep.dendro && doCdend) ddc))
}

