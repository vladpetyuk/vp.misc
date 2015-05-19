#' Nested Linear Model Analysis
#' 
#' A function for analysis of sigfinicance of factors (akin ANOVA)
#' or continuous variable (regression) using nested linear models.
#' Shuffling-based FDR estimation is available as an option.
#' 
#' @param eset eset (or most likely eset subclass) object
#' @param form.alt character formulation of alternative model
#' @param form.nul character formulation of NULL model
#' @param facs data frame with the factors in its columns.
#'          If NULL, then pData(eset) will be used.
#' @param norm.coef vector with sample-to-sample normalization 
#'          coefficients in log2 scale.
#' @param N number of shuffles for FDR estimation. Default is NULL, that is
#'          no shuffling-based FDR estimation. Warning! Will be N-times slower.
#'          This is critical for Windows, since paralellization is not impletemted.
#'          In case of enabling shuffling-based FDR estimation, 
#'          the recommended value is at least N=1000.
#' 
#' @return data.frame
#'      \describe{
#'          \item{\code{effect}}{either max contrast in case of factor 
#'                               or slope in case of continuous variable}
#'          \item{\code{F.stat}}{F statistic}
#'          \item{\code{p.value}}{p-value}
#'      }
#' @importFrom Biobase exprs pData featureNames
#' @export eset_lm
#' 
#' @examples
#' data(srm_msnset)
#' head(varLabels(msnset))
#' out <- eset_lm(msnset, "y ~ subject.type", "y ~ 1")
#' head(out)
#' # now with shuffling
#' out <- eset_lm(msnset, "y ~ subject.type", "y ~ 1", N = 100)
#' head(out)
#' head(out[order(out$p.value),])


eset_lm <- function (eset, form.alt, form.nul, facs = NULL, norm.coef = NULL, N=NULL) 
{
    if(is.null(N))
        return(eset_lm_one_pass(eset, form.alt, form.nul, facs, norm.coef))
    else
        return(eset_lm_shuffled(eset, form.alt, form.nul, facs, norm.coef, N))
}




eset_lm_one_pass <- function (eset, form.alt, form.nul, facs = NULL, norm.coef = NULL) 
{
    # Check models. They should be nested and differ only in one term.
    stopifnot(length(setdiff(all.vars(as.formula(form.alt)),
                             all.vars(as.formula(form.nul)))) == 1)
    #
    data <- exprs(eset)
    if (is.null(facs)) 
        facs <- pData(eset)
    if (is.null(norm.coef)) 
        norm.coef <- rep(0, ncol(data)) # log2 scale
    Np <- nrow(data)
    res <- data.frame(effect = numeric(Np), 
                         F.stat = numeric(Np), 
                         p.value = numeric(Np))
    # make parallel
    if(.Platform$OS.type == "unix")
        mc.cores <- parallel::detectCores()
    else
        mc.cores <- 1
    res <- parallel::mclapply(seq_len(nrow(eset)), function(i){
        lm_for_one(data[i,], form.alt, form.nul, 
                   facs, norm.coef)}, 
        mc.cores=mc.cores)
    res <- Reduce(rbind, res)
    
    #     for (i in seq_len(Np)) 
    #         result[i, ] <- lm.for.one(data[i,], form.alt, form.nul, 
    #                                   facs, norm.coef)
    
    rownames(res) <- rownames(data)
    res <- as.data.frame(res) # coerce to specific type in case of mat/df ambig
    res$fdr <- p.adjust(res$p.value, method="fdr")
    return(res)
}



eset_lm_shuffled <- function(eset, form.alt, form.nul, 
                             facs = NULL, norm.coef = NULL, N=1000)
{
    # Check models. They should be nested and differ only in one term.
    stopifnot(length(setdiff(all.vars(as.formula(form.alt)),
                             all.vars(as.formula(form.nul)))) == 1)
    #
    # normal test here
    res <- eset_lm(eset, form.alt, form.nul, facs, norm.coef)
    #
    sh.res.mat <- matrix(NA, ncol=N, nrow=nrow(eset))
    for(i in seq_len(N)){
        # shuffle here
        sh.idx <- sample(seq_len(ncol(eset)))
        # nested LM test
        res.sh <- eset_lm(eset[,sh.idx], form.alt, form.nul, facs, norm.coef)
        # count
        sh.res.mat[,i] <- sapply(res$p.value, function(x) sum(res.sh$p.value <= x))
    }
    #
    trufls <- sapply(res$p.value, function(x) sum(res$p.value <= x))  # T + F
    fls <- rowSums(sh.res.mat)/N  # F
    res$fdr.sh <- fls/trufls
    res$fdr.sh[out$fdr.sh > 1] <- 1
    rownames(res) <- featureNames(eset)
    return(res)
}



lm_for_one <- function(ints, form.alt, form.nul, facs, off) {
    #
    data <- data.frame(y = ints, facs, off)
    selected.rows <- rownames(model.frame(as.formula(form.alt), data=data))
    data <- data[selected.rows,]
    #
    mod.alt <- lm(as.formula(form.alt), offset = off, data = data)
    mod.nul <- lm(as.formula(form.nul), offset = off, data = data)
    anstat  <- anova(mod.alt, mod.nul, test = "F")
    p.value <- anstat[2, 'Pr(>F)']
    F.stat <- anstat[2, 'F']
    # now extract the effect
    # detemine is this is multifactorial ANOVA
    dif.term <- setdiff(all.vars(as.formula(form.alt)), 
                        all.vars(as.formula(form.nul)))
    if(is.factor(data[[dif.term]]) & nlevels(data[[dif.term]]) > 2){
        # basically I want to pull out maximum contrast
        col.pos <- grep(dif.term, names(coef(mod.alt)))
        effect <- coef(mod.alt)[col.pos]
        if("(Intercept)" %in% names(coef(mod.alt)))
            effect <- c(0, effect)
        effect <- diff(range(c(effect)))
    }else{
        col.pos <- grep(dif.term, names(coef(mod.alt)))
        effect <- coef(mod.alt)[col.pos]
    }
    return(c(effect=effect, F.stat=F.stat, p.value=p.value))
}



