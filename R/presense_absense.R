#' Presence/Absence Testing
#'
#' A function for testing the significance in differences
#' of presence/absence patterns.
#'
#' The recommended way of using the test is two-group comparison.
#' All except "binom" options will work for larger number of groups.
#'
#'
#' @param eset eset (or most likely eset subclass) object
#' @param grouping character defining the column in phenoData
#' @param test character: Fisher exact test, \eqn{\chi^{2}}{\chi^2}, G-test or binomial.
#' @param ... other aruments to pass to the test functions
#'              (e.g. simulate.p.value = TRUE)
#'
#' @return data.frame
#'      \describe{
#'          \item{\code{p.value}}{p-value}
#'          \item{\code{effect}}{difference between second and first group
#'                              (according to the factor levels).
#'                              In case of larger number of levels,
#'                              it is the absolute value of the max diffence}
#'          \item{\code{...}}{proportions of present for each group.
#'                              Column names are level names.}
#'      }
#' @note see also \url{http://www.ncbi.nlm.nih.gov/pubmed/20831241}
#'
#' @importFrom Biobase exprs pData featureNames rowMax rowMin
#' @importFrom stats fisher.test chisq.test binom.test
#'
#' @export eset_presence_absence
#'
#' @examples
#' library("MSnbase")
#' Nsam = 20
#' Npep = 5
#' set.seed(99)
#' M <- matrix(rbinom(Nsam * Npep, 1, 0.5), ncol = Nsam)
#' pd <- data.frame(group = gl(2, Nsam/2))
#' fd <- data.frame(otherfdata = letters[1:Npep])
#' x0 <- MSnSet(M, fd, pd)
#'
#' eset_presence_absence(x0, 'group', test = 'fisher')
#' eset_presence_absence(x0, 'group', test = 'chisq', simulate.p.value = TRUE)
#' eset_presence_absence(x0, 'group', test = 'g')
#' eset_presence_absence(x0, 'group', test = 'binom')


eset_presence_absence <- function(eset,
                                  grouping,
                                  test = c('fisher', 'chisq', 'g', 'binom'),
                                  ...) {

    test <- match.arg(test, choices = c('fisher', 'chisq', 'g', 'binom'))

    if(is.character(grouping) & length(grouping) == 1)
        grouping <- pData(eset)[[grouping]]

    stopifnot(is.factor(grouping))

    ## Presence count
    # Sum by row within each group
    pc <- t(apply(exprs(eset), 1, tapply, INDEX = grouping, sum))
    total.pc <- colSums(pc) # Overall group totals
    N <- table(grouping) # Size of each group

    # Get the p-values
    if(test != 'binom'){
        TEST.FUN <- switch(test,
                           fisher = fisher.test,
                           chisq = chisq.test,
                           g = g.test)
        out <- apply(pc, 1, function(x) {
            cm <- rbind(x, total.pc - x)
            res <- TEST.FUN(cm, ...)
            res$p.value
        })
    } else { # binomial test
        if(nlevels(grouping) > 2)
            warning(paste("binom option will test only the difference between",
                          "the first two levels of grouping"))
        p <- pc[, 1] / N[1]
        out <- sapply(seq_len(nrow(pc)), function(i){
            res <- binom.test(pc[i,2], N[2], p[i])
            res$p.value})
    }

    # Estimate the effect
    est <- sweep(pc, 2, N, '/')
    colnames(est) <- levels(grouping)

    if(nlevels(grouping) > 2)
        effect <- abs(rowMax(est) - rowMin(est))

    if(nlevels(grouping) == 2)
        effect <- est[, 2] - est[, 1]

    # Data frame with column for p-values, effects, and proportion
    # of presence for each group
    out <- data.frame(p.value = out, effect = effect, est)
    rownames(out) <- featureNames(eset)

    return(out)
}


# The g.test function is not exported by Biostrings.
# It would be best to copy the function to this package, but
# we need to check the Biostrings LICENSE first. For now,
# do this to remove the note about :::
g.test <- utils::getFromNamespace("g.test", "Biostrings")



