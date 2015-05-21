


#' Tweaks of plot.annHeatmap
#' 
#' Tweaks of \code{\link[Heatplus]{plot.annHeatmap}} method for 
#' \code{annHeatmap} class.  Additional argument for defining color 
#' for missing values.  There is also a hack of 
#' \code{\link[Heatplus]{picketPlot}} that enforce axis labels 
#' to be always perpendicular.
#' 
#' @param na.color character color of NA values
#' @export plot_annHeatmap2
#' 

plot_annHeatmap2 <- function (x, na.color='lightgrey', widths, heights, ...) 
{
    doRClusLab = !is.null(x$cluster$Row$labels)
    doCClusLab = !is.null(x$cluster$Col$labels)
    omar = rep(0, 4)
    if (doRClusLab) 
        omar[1] = 2
    if (doCClusLab) 
        omar[4] = 2
    par(oma = omar)
    if (!missing(widths)) 
        x$layout$width = widths
    if (!missing(heights)) 
        x$layout$height = heights
    with(x$layout, layout(plot, width, height, respect = TRUE))
    nc = ncol(x$data$x2)
    nr = nrow(x$data$x2)
    doRlab = !is.null(x$labels$Row$labels)
    doClab = !is.null(x$labels$Col$labels)
    mmar = c(1, 0, 0, 2)
    if (doRlab) 
        mmar[x$labels$Row$side] = x$labels$Row$nrow
    if (doClab) 
        mmar[x$labels$Col$side] = x$labels$Col$nrow
    with(x$data, {
        par(mar = mmar)
        image(1:nc, 1:nr, t(x2), axes = FALSE, 
                xlim = c(0.5, nc + 0.5), ylim = c(0.5, nr + 0.5), 
                xlab = "", ylab = "", 
                col = col, breaks = breaks, ...)
        mmat <- ifelse(is.na(t(x2)), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, 
                xlab = "", ylab = "", 
                col = na.color, add=TRUE)
    })
    with(x$labels, {
        if (doRlab) 
            axis(Row$side, 1:nr, las = 2, line = -0.5, tick = 0, 
                  labels = Row$labels, cex.axis = Row$cex)
        if (doClab) 
            axis(Col$side, 1:nc, las = 2, line = -0.5, tick = 0, 
                  labels = Col$labels, cex.axis = Col$cex)
    })
    with(x$dendrogram$Col, if (status == "yes") {
        par(mar = c(0, mmar[2], 3, mmar[4]))
        Heatplus:::cutplot.dendrogram(dendro, h = x$cluster$Col$cuth, cluscol = x$cluster$Col$col, 
                                 horiz = FALSE, axes = FALSE, xaxs = "i", leaflab = "none", 
                                 lwd = x$dendrogram$Col$lwd)
    })
    with(x$dendrogram$Row, if (status == "yes") {
        par(mar = c(mmar[1], 3, mmar[3], 0))
        Heatplus:::cutplot.dendrogram(dendro, h = x$cluster$Row$cuth, cluscol = x$cluster$Row$col, 
                                 horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none", 
                                 lwd = x$dendrogram$Row$lwd)
    })
    if (!is.null(x$annotation$Col$data)) {
        par(mar = c(1, mmar[2], 0, mmar[4]), xaxs = "i", yaxs = "i")
        picketPlot2(x$annotation$Col$data[x$data$colInd, , drop = FALSE], 
                      grp = x$cluster$Col$grp, 
                      grpcol = x$cluster$Col$col, control = x$annotation$Col$control, 
                      asIs = TRUE)
    }
    if (!is.null(x$annotation$Row$data)) {
        par(mar = c(mmar[1], 0, mmar[3], 1), xaxs = "i", yaxs = "i")
        picketPlot2(x$annotation$Row$data[x$data$rowInd, , drop = FALSE], 
                      grp = x$cluster$Row$grp, 
                      grpcol = x$cluster$Row$col, control = x$annotation$Row$control, 
                      asIs = TRUE, horizontal = FALSE)
    }
    if (x$legend) {
        if (x$layout$legend.side %in% c(1, 3)) {
            par(mar = c(2, mmar[2] + 2, 2, mmar[4] + 2))
        }
        else {
            par(mar = c(mmar[1] + 2, 2, mmar[3] + 2, 2))
        }
        Heatplus::doLegend(x$data$breaks, col = x$data$col, x$layout$legend.side)
    }
    invisible(x)
}




picketPlot2 <- function (x, grp = NULL, grpcol, grplabel = NULL, 
                         horizontal = TRUE, asIs = FALSE, control = list()) 
{
    cc = Heatplus:::picketPlotControl()
    cc[names(control)] = control
    x = Heatplus::convAnnData(x, asIs = asIs)
    nsamp = nrow(x)
    npanel = ncol(x)
    bpanel = apply(x, 2, function(y) all(y[is.finite(y)] %in% 
                                             c(0, 1)))
    panelw = nsamp * (cc$boxw + 2 * cc$hbuff)
    panelh = cc$boxh + 2 * cc$vbuff
    totalh = sum(panelh * ifelse(bpanel, 1, cc$numfac))
    LL = cbind(0, 0)
    UR = cbind(panelw, totalh)
    xbase = seq(cc$hbuff, by = cc$boxw + 2 * cc$hbuff, length = nsamp)
    xcent = xbase + cc$boxw/2
    if (!is.null(grp)) {
        grp = as.integer(factor(grp, levels = unique(grp)))
        tt = table(grp)
        gg = length(tt)
        grpcoord = c(0, cumsum(tt/sum(tt)) * panelw)
        grp0 = cbind(grpcoord[1:gg], rep(0, gg))
        grp1 = cbind(grpcoord[2:(gg + 1)], rep(totalh, gg))
        if (missing(grpcol)) {
            grpcol = Heatplus::RainbowPastel
        }
        if (is.function(grpcol)) 
            grpcol = grpcol(gg)
    }
    panels = list()
    voff = 0
    for (i in 1:npanel) {
        if (bpanel[i]) {
            x0 = xbase
            x1 = x0 + cc$boxw
            y0 = voff + cc$vbuff
            y1 = y0 + cc$boxh
            fill = ifelse(x[, i, drop = FALSE] == 1, "black", 
                          "transparent")
            fill[is.na(fill)] = cc$nacol
            label = colnames(x)[i]
            labcc = if (!is.null(label)) 
                (y0 + y1)/2
            else NULL
            panels[[i]] = list(ll = cbind(x0, y0), ur = cbind(x1, y1), 
                               fill = fill, label = label, labcc = labcc)
            voff = voff + panelh
        }
        else {
            xv = x[, i]
            rr = range(xv, na.rm = TRUE)
            yval = voff + cc$vbuff * cc$numfac + 
                ((xv - rr[1])/(rr[2] - rr[1])) * cc$boxh * cc$numfac
            if ((cc$degree > 0) & (cc$span > 0)) {
                yy = predict(loess(yval ~ xcent, span = cc$span, 
                                   degree = cc$degree))
            }
            else {
                yy = rep(NA, length(xcent))
            }
            label = colnames(x)[i]
            labcc = if (!is.null(label)) 
                mean(range(yval, na.rm = TRUE))
            else NULL
            axlab = pretty(range(xv, na.rm = TRUE))
            axcc = voff + cc$vbuff * cc$numfac + 
                ((axlab - rr[1])/(rr[2] - rr[1])) * cc$boxh * cc$numfac
            panels[[i]] = list(raw = cbind(xcent, yval), smo = cbind(xcent, yy), 
                               label = label, labcc = labcc, axlab = axlab, 
                               axcc = axcc)
            voff = voff + panelh * cc$numfac
        }
    }
    if (!is.null(grp) & !is.null(grplabel)) {
        mids = (grpcoord[1:gg] + grpcoord[2:(gg + 1)])/2
        labelnum = length(grplabel)
        if (labelnum < gg) {
            warning("more groups than labels (filling up with blanks)")
            grplabel = c(grplabel, rep(" ", gg - labelnum))
        }
        else if (gg < labelnum) {
            warning("more labels than groups (ignoring the extras)")
            grplabel = grplabel[1:gg]
        }
    }
    h2v = function(cc) cbind(cc[, 2] - totalh, cc[, 1])
    if (horizontal) {
        grpaxis = 1
        labaxis = 2
        covaxis = 4
        las = 1
    }
    else {
        grpaxis = 4
        labaxis = 3
        covaxis = 1
        las = 3
        LL = h2v(LL)
        UR = h2v(UR)
        if (!is.null(grp)) {
            grp0 = h2v(grp0)
            grp1 = h2v(grp1)
        }
        for (i in 1:npanel) {
            panels[[i]][[1]] = h2v(panels[[i]][[1]])
            panels[[i]][[2]] = h2v(panels[[i]][[2]])
            panels[[i]]$labcc = panels[[i]]$labcc - totalh
            panels[[i]]$axcc = panels[[i]]$axcc - totalh
        }
    }
    plot(rbind(LL, UR), type = "n", xaxt = "n", yaxt = "n", xlab = "", 
         ylab = "")
    if (!is.null(grp)) {
        rect(grp0[, 1], grp0[, 2], grp1[, 1], grp1[, 2], col = grpcol, 
             border = "transparent")
    }
    for (i in 1:npanel) {
        if (bpanel[i]) {
            with(panels[[i]], rect(ll[, 1], ll[, 2], ur[, 1], 
                                   ur[, 2], col = fill, border = "transparent"))
        }
        else {
            with(panels[[i]], points(raw[, 1], raw[, 2], pch = cc$pch, 
                                     cex = cc$cex.pch, col = cc$col.pch))
            if ((cc$degree > 0) & (cc$span > 0)) {
                with(panels[[i]], lines(smo[, 1], smo[, 2]))
            }
            with(panels[[i]], axis(covaxis, at = axcc, labels = axlab, las=2))
        }
        if (!is.null(panels[[i]]$label)) {
            axis(labaxis, at = panels[[i]]$labcc, labels = panels[[i]]$label, 
                 las = las, tick = FALSE, font = 2, col = par("bg"), 
                 col.axis = par("fg"))
        }
    }
    if (!is.null(grp) & !is.null(grplabel)) {
        axis(grpaxis, grpcoord, labels = FALSE, tcl = -1.5)
        axis(grpaxis, mids, labels = grplabel, font = 2, cex.axis = cc$cex.label, 
             tick = FALSE)
    }
    invisible(panels)
}

