#' Nested Linear Model Analysis
#' 
#' A function for analysis of sigfinicance of factors (akin ANOVA)
#' or continuous variable (regression) using nested linear models.
#' The function modelled after `msmsTests::msms.glm.qlll`
#' 
#' @param eset eset (or most likely eset subclass) object
#' @param form.alt character formulation of alternative model
#' @param form.nul character formulation of NULL model
#' @param facs data frame with the factors in its columns.
#'          If NULL, then pData(eset) will be used.
#' @param norm.coef vector with sample-to-sample normalization 
#'          coefficients in log2 scale.
#' @return data.frame
#'      \describe{
#'          \item{\code{effect}}{either max contrast in case of factor 
#'                               or slope in case of continuous variable}
#'          \item{\code{F.stat}}{F statistic}
#'          \item{\code{p.value}}{p-value}
#'      }
#' @export eset_lm
#' 
eset_lm <- function (eset, form.alt, form.nul, facs = NULL, norm.coef = NULL) 
{
    data <- exprs(eset)
    if (is.null(facs)) 
        facs <- pData(eset)
    if (is.null(norm.coef)) 
        norm.coef <- rep(0, ncol(data)) # log2 scale
    Np <- nrow(data)
    result <- data.frame(effect = numeric(Np), 
                         F.stat = numeric(Np), 
                         p.value = numeric(Np))
    # make parallel
    
    result <- parallel::mclapply(seq_len(nrow(eset)), function(i){
        lm_for_one(data[i,], form.alt, form.nul, 
                   facs, norm.coef)}, 
        mc.cores=parallel::detectCores())
    result <- Reduce(rbind, result)
    
    #     for (i in seq_len(Np)) 
    #         result[i, ] <- lm.for.one(data[i,], form.alt, form.nul, 
    #                                   facs, norm.coef)
    
    rownames(result) <- rownames(data)
    return(as.data.frame(result)) # coerce to specific type in case of mat/df ambig
}





lm_for_one <- function(ints, form.alt, form.nul, facs, off) {
    data <- data.frame(y = ints, facs, off)
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



