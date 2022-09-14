#' This is a modified version of the PeCorA function.
#' For the original paper see here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8592057/.
#'
#' The modification involves aggregating the intensities of "all other peptides" within
#' the test described in the PeCorA paper. We aggregate "all other peptides" by sample using the median, this
#' way each peptide mapping to a particular protein is not treated as an 'independent' observation.
#' The effect is that p-values are much higher (less significant) than in the original version of PeCorA.
#'
#' This median modification is turned OFF by default, and only affects the p-values reported by PeCorA.
#'
#'
#' @param t PeCorA input table. See the original PeCorA function for details.
#' @param median_mod logical for whether to use the median version of PeCorA (median_mod = TRUE), or the original (median_mod = FALSE)
#'
#' @return PeCorA results table.
#'
#' @importFrom stats lm
#'
#' @export PeCorA_mod
#'


PeCorA_mod <- function (t, median_mod = FALSE) {

  print("checking which proteins still have at least 2 peptides")
  pgs <- levels(as.factor(t$Protein))
  pgs_morethan2 <- c()
  # Extremely slow
  # for (x in pgs) {
  #   if (length(unique(t[t$Protein %in% x, "modpep_z"])) >
  #       1) {
  #     pgs_morethan2 <- c(pgs_morethan2, x)
  #   }
  # }
  pgs_morethan2 <- t %>%
    select(Protein, modpep_z) %>%
    unique() %>%
    group_by(Protein) %>%
    summarize(n_peptides = n()) %>%
    filter(n_peptides > 1) %>%
    pull(Protein)

  allp <- list()
  j = 1
  t0 <- Sys.time()
  print("computing the interaction p-values")
  pb <- txtProgressBar(min = 0, max = length(pgs_morethan2),
                       style = 3)
  for (x in pgs_morethan2) {
    tmpdf <- t[which(t$Protein == x), ]
    tmpdf["allothers"] <- rep("allothers", times = nrow(tmpdf))
    pvalues <- c(rep(0, length(unique(tmpdf$modpep_z))))
    i = 1
    for (y in unique(tmpdf$modpep_z)) {
      subtmpdf <- tmpdf
      subtmpdf[which(tmpdf$modpep_z == y), "allothers"] <- y
      if (median_mod){
        subtmpdf <- subtmpdf %>%
          group_by(Sample, allothers) %>%
          mutate(ms1adj = median(ms1adj, na.rm = T)) %>%
          ungroup() %>%
          select(Sample, ms1adj, Condition, allothers) %>%
          unique()
      }
      tmplm <- lm(subtmpdf$ms1adj ~ subtmpdf$Condition *
                    subtmpdf$allothers)
      tmpanova <- car::Anova(tmplm)
      pvalues[i] <- tmpanova$`Pr(>F)`[3]
      i = i + 1
    }
    allp[[x]] <- pvalues
    setTxtProgressBar(pb, j)
    j = j + 1
  }
  print(" ")
  print(paste("PeCorA finished in ", round(Sys.time() - t0,
                                           2), " minutes", sep = ""))
  print(paste("number of proteins tested =", length(allp),
              sep = " "))
  print(paste("number of peptides tested =", length(unlist(allp)),
              sep = " "))
  print("started making data table")
  alldf = data.frame()
  x <- names(allp)[1]
  for (x in names(allp)) {
    tmpdf <- t[which(t$Protein == x), ]
    tmp_peps <- unique(tmpdf$modpep_z)
    if (length(tmp_peps) > 0) {
      tmp_pval <- allp[[x]]
      tmpout = cbind.data.frame(protein = rep(x, length(allp[[x]])),
                                tmp_peps, tmp_pval = as.numeric(tmp_pval))
      alldf = rbind(alldf, tmpout)
    }
  }
  print("correcting p-values")
  alldf$adj_pval <- p.adjust(alldf$tmp_pval, method = "BH")
  alldf_ordered <- alldf[order(alldf$adj_pval), ]
  print(paste("number of uncorrelated peptides =", nrow(alldf[alldf$adj_pval <=
                                                                0.01, ]), sep = " "))
  print(paste("number of proteins with uncorrelated peptides =",
              length(unique(alldf[alldf$adj_pval <= 0.01, ]$protein)),
              sep = " "))
  colnames(alldf_ordered)[2] <- "peptide"
  colnames(alldf_ordered)[3] <- "pvalue"
  alldf_ordered
}
