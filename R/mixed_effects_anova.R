#' Mixed Effects ANOVA
#'
#' A wrapper around \code{\link[lme4]{lmer}} functionality.
#'
#' @param m ExpressionSet or MSnSet object
#' @param alt_model character formulation of alternative model
#' @param test_factor name of the tested factor
#' @param order_by "none" or "p-value". Default is "none".
#'
#' @return data.frame
#'      \describe{
#'          \item{\code{feature}}{protein or feature names in general.
#'          Should be as in sampleNames(m).}
#'          \item{\code{p.val}}{ANOVA p-value}
#'          \item{\code{p.adj}}{Benjamini-Hochberg adjusted p-value}
#'          \item{\code{q.val}}{q-value.
#'          Parameters for qvalue() are: lambda=0.05, pi0.method="bootstrap"}
#'      }
#' @importFrom lme4 lmer
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr inner_join do select summarize mutate
#'                   arrange matches group_by %>%
#' @importFrom tidyr gather
#' @importFrom Biobase exprs pData featureNames
#' @importFrom qvalue qvalue
#' @importFrom stats as.formula anova p.adjust
#'
#' @export mixed_effects_anova
#'
#' @examples
#' data(srm_msnset)
#' # -- make sure the factors are factors (not characters or integers)
#' # This is intentionally left for user's control.
#' msnset$subject.type <- as.factor(msnset$subject.type)
#' msnset$match.group <- as.factor(msnset$match.group)
#' msnset$PrepBatch <- as.factor(msnset$PrepBatch)
#' # -- model description
#' # subject.type - fixed effect
#' # match.group - random effect affecting the intercept
#' # PlateCol - another random effect influencing intercept
#' out <- mixed_effects_anova(msnset,
#'                            alt_model="~ subject.type + (1|match.group) + (1|PrepBatch)",
#'                            test_factor="subject.type",
#'                            order_by="p-value")
#' head(out)


mixed_effects_anova <- function(m, alt_model, test_factor,
                                      order_by=c("none","p-value")){

    # 1 check that factor is part of model
    factor_is_in_model <- grepl(test_factor, alt_model)
    stopifnot(factor_is_in_model)

    # 2 check that there are no characters
    vrs <- all.vars(as.formula(alt_model))
    no_characters_in_vars <- all(!sapply(pData(m)[,vrs],is.character))
    stopifnot(no_characters_in_vars)

    # 3 check ordering argument
    order_by <- match.arg(order_by, c("none","p-value"))

    xe <- exprs(m) %>%
        as.data.frame() %>%
        rownames_to_column("feature") %>%
        gather(sample_name, abundance, -feature)
    xp <- pData(m) %>%
        select(-matches("^sample_name$")) %>%
        rownames_to_column("sample_name")
    x <- inner_join(xe, xp, by="sample_name")

    null_model <- sub(test_factor, "1", alt_model)
    alt_model <- paste("abundance", alt_model)
    null_model <- paste("abundance", null_model)

    safe_qvalue <- function(x){
        tryCatch(qvalue(x, lambda=0.05, pi0.method="bootstrap")$qvalue,
                          error=function(e) NA)
    }

    res <- x %>%
        group_by(feature) %>%
        do(
            alt = lmer(eval(parse(text=alt_model)), data=., REML=F),
            nul = lmer(eval(paste(text=null_model)), data=., REML=F)) %>%
        summarize(p.val = anova(alt,nul)[2,'Pr(>Chisq)'],
                  feature = feature) %>%
        mutate(p.adj = p.adjust(p.val, method="BH"),
               q.val = safe_qvalue(p.val)) %>%
        select(feature, p.val, p.adj, q.val)
    if(order_by == "p-value")
        res <- res %>% arrange(p.val)
    return(res)
}

utils::globalVariables(c("sample_name", "abundance", "feature", ".",
                         "alt", "nul", "p.val", "p.adj", "q.val"))

