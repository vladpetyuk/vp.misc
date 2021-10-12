
#' Adjust for batch effects using an empirical Bayes framework
#'
#' This is a modified version of ComBat, which can handle when all values in a feature + batch
#' combination are missing. It does this by ignoring missing values when calculating
#' prior distributions for the batch effects, and adjusting the values that ARE there, while doing
#' nothing to values that are missing. Particularly useful in data with many batches, as previously
#' a feature with only one batch fully missing would have been thrown out before correction.
#'
#' Can also retain specific variable effects (use mod parameter for this), and estimate the final batch effects
#' using parametric or non-parametric adjustments (paramtric is default).
#'
#' ComBat allows users to adjust for batch effects in datasets where the batch covariate is known, using methodology
#' described in Johnson et al. 2007. It uses either parametric or non-parametric empirical Bayes frameworks for adjusting data for
#' batch effects.  Users are returned an expression matrix that has been corrected for batch effects. The input
#' data are assumed to be cleaned and normalized before batch effect removal.
#'
#' @param dat Genomic measure matrix (dimensions probe x sample) - for example, expression matrix
#' @param batch {Batch covariate (only one batch allowed)}
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' @param par.prior (Optional) TRUE indicates parametric adjustments will be used, FALSE indicates non-parametric adjustments will be used
#' @param prior.plots (Optional) TRUE give prior plots with black as a kernel estimate of the empirical batch effect density and red as the parametric
#' @param mean.only (Optional) FALSE If TRUE ComBat only corrects the mean of the batch effect (no scale adjustment)
#' @param ref.batch (Optional) NULL If given, will use the selected batch as a reference for batch adjustment.
#'
#' @return data A probe x sample genomic measure matrix, adjusted for batch effects.
#'
#' @importFrom limma lmFit
#' @importFrom graphics lines par
#' @importFrom genefilter rowVars
#' @importFrom stats cor density dnorm model.matrix pf ppoints prcomp predict qgamma qnorm qqline qqnorm qqplot smooth.spline var
#' @importFrom utils read.delim
#'
#' @examples
#' library(bladderbatch)
#' data(bladderdata)
#' dat <- bladderEset[1:50,]
#'
#' pheno = pData(dat)
#' edata = exprs(dat)
#' batch = pheno$batch
#' mod = model.matrix(~as.factor(cancer), data=pheno)
#'
#' # parametric adjustment
#' combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
#'
#' # non-parametric adjustment, mean-only version
#' combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)
#'
#' # reference-batch version, with covariates
#' combat_edata3 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, ref.batch=3)
#'
#' @export ComBat.NA
#'

ComBat.NA <- function(dat, batch, mod = NULL, par.prior = TRUE, mean.only = FALSE,
                      prior.plots = FALSE, ref.batch = NULL) {

  if(length(dim(batch))>1){
    stop("This version of ComBat only allows one batch variable")
  }  ## to be updated soon!

  ## coerce dat into a matrix
  dat <- as.matrix(dat)

  ## find genes with zero variance in any of the batches
  batch <- as.factor(batch)
  # zero.rows.lst <- lapply(levels(batch), function(batch_level){
  #   if(sum(batch==batch_level)>1){
  #     return(which(apply(dat[, batch==batch_level], 1, function(x){var(x)==0})))
  #   }else{
  #     return(which(rep(1,3)==2))
  #   }
  # })
  zero.rows.lst <- lapply(levels(batch), function(batch_level){
    if(sum(batch==batch_level)>1){
      return(which(apply(dat[, batch==batch_level], 1, function(x){
        if (all(is.na(x))){
          F
        } else {
          var(x, na.rm = T) %in% c(0, NA)
        }
      })))
    }else{
      return(which(rep(1,3)==2))
    }
  })
  zero.rows <- Reduce(union, zero.rows.lst)
  keep.rows <- setdiff(1:nrow(dat), zero.rows)

  if (length(zero.rows) > 0) {
    cat(sprintf("Found %d genes with uniform expression within a single batch, leading to a variance of zero;
                these will not be adjusted for batch.\n", length(zero.rows)))
    # keep a copy of the original data matrix and remove zero var rows
    dat.orig <- dat
    dat <- dat[keep.rows, ]
  }

  if (length(keep.rows) <= 2) {
    cat("Within batch variance is zero for all but (at most) 2 features!")
    stop("Within batch variance is zero for all but (at most) 2 features!")
  }

  ## make batch a factor and make a set of indicators for batch
  if(any(table(batch)==1)){mean.only=TRUE}
  if(mean.only==TRUE){
    message("Using the 'mean only' version of ComBat")
  }

  batchmod <- model.matrix(~-1+batch)
  if (!is.null(ref.batch)){
    ## check for reference batch, check value, and make appropriate changes
    if (!(ref.batch%in%levels(batch))) {
      stop("reference level ref.batch is not one of the levels of the batch variable")
    }
    message("Using batch =",ref.batch, "as a reference batch (this batch won't change)")
    ref <- which(levels(as.factor(batch))==ref.batch) # find the reference
    batchmod[,ref] <- 1
  } else {
    ref <- NULL
  }
  message("Found", nlevels(batch), "batches")

  ## A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  } # list of samples in each batch
  n.batches <- sapply(batches, length)
  if(any(n.batches==1)){
    mean.only=TRUE
    message("Note: one batch has only one sample, setting mean.only=TRUE")
  }
  n.array <- sum(n.batches)
  ## combine batch variable and covariates
  design <- cbind(batchmod,mod)

  ## check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  if(!is.null(ref)){
    check[ref] <- FALSE
  } ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])

  ## Number of covariates or covariate levels
  message("Adjusting for", ncol(design)-ncol(batchmod), 'covariate(s) or covariate level(s)')

  ## Check if the design is confounded
  if(qr(design)$rank < ncol(design)) {
    ## if(ncol(design)<=(n.batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n.batch+1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
    }
    if(ncol(design)>(n.batch+1)) {
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
        stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      } else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
      }
    }
  }

  ## Check for missing values
  NAs <- any(is.na(dat))
  if(NAs){
    message(c('Found',sum(is.na(dat)),'Missing Data Values'), sep=' ')}
  ## print(dat[1:2,])

  ##Standardize Data across genes
  message('Standardizing Data across genes')
  if (!NAs){
    B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(dat)))
  } else {
    # B.hat <- apply(dat, 1, Beta.NA, design) # FIXME

    ### Computing coefficients using lmFit, as this can deal with NA coefficients.
    B.hat <- lmFit(dat, design)$coefficients
    B.hat <- t(as.matrix(B.hat))
  }

  B.hat.0 <- B.hat
  B.hat.0[is.na(B.hat.0)] <- 0

  ## change grand.mean for ref batch
  if(!is.null(ref.batch)){
    grand.mean <- t(B.hat[ref, ])
  } else {
    ## Computes the mean of the batch means. This is similar (but not EXACTLY equal) to the overall mean
    # grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])

    # Adapted for NA values
    grand.mean <- t(B.hat[1:n.batch, ]) %>%
      as.data.frame() %>%
      apply(1, mean, na.rm = T) %>%
      as.matrix() %>%
      t()
  }


  ## change var.pooled for ref batch
  if (!NAs){
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- ((ref.dat-t(design[batches[[ref]], ] %*% B.hat))^2) %*% rep(1/n.batches[ref],n.batches[ref]) # FIXME
    } else {
      var.pooled <- ((dat-t(design %*% B.hat))^2) %*% rep(1/n.array,n.array) # FIXME
    }
  } else {
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- rowVars(ref.dat-t(design[batches[[ref]], ]%*%B.hat), na.rm=TRUE)
    } else {
      # var.pooled <- rowVars(dat-t(design %*% B.hat), na.rm=TRUE)

      # Adapted for NA values
      var.pooled <- rowVars(dat - t(design %*% B.hat.0), na.rm = TRUE)
    }
  }


  ### mean over each gene plus retained covariate effects
  stand.mean <- t(grand.mean) %*% t(rep(1,n.array)) # FIXME
  if(!is.null(design)){
    tmp <- design
    tmp[,c(1:n.batch)] <- 0
    stand.mean <- stand.mean+t(tmp %*% B.hat.0) #FIXME
  }

  ### This is the Z_{ijg} in the original paper.
  ### Removing the mean and retained covariates when standardizing
  s.data <- (dat-stand.mean)/(sqrt(var.pooled) %*% t(rep(1,n.array))) # FIXME

  ##Get regression batch effect parameters
  message("Fitting L/S model and finding priors")
  batch.design <- design[, 1:n.batch]
  if (!NAs){
    gamma.hat <- solve(crossprod(batch.design), tcrossprod(t(batch.design),
                                                           as.matrix(s.data)))
  } else{
    # gamma.hat <- apply(s.data, 1, Beta.NA, batch.design) # FIXME

    # Adapted for NA values
    gamma.hat <- lmFit(s.data, batch.design)$coefficients
    gamma.hat <- t(as.matrix(gamma.hat))
  }
  delta.hat <- NULL
  for (i in batches){
    if(mean.only==TRUE) {
      delta.hat <- rbind(delta.hat,rep(1,nrow(s.data)))
    } else {
      delta.hat <- rbind(delta.hat, rowVars(s.data[,i], na.rm=TRUE))
    }
  }

  ##Find Priors
  # gamma.bar <- rowMeans(gamma.hat)
  gamma.bar <- rowMeans(gamma.hat, na.rm = TRUE)
  # t2 <- rowVars(gamma.hat)
  t2 <- rowVars(gamma.hat, na.rm = TRUE)

  aprior.na <- function(delta.hat) {
    m <- mean(delta.hat, na.rm = T)
    s2 <- var(delta.hat, na.rm= T)
    (2 * s2 + m^2)/s2
  }

  bprior.na <- function(delta.hat) {
    m <- mean(delta.hat, na.rm = T)
    s2 <- var(delta.hat, na.rm= T)
    (m * s2 + m^3)/s2
  }

  # a.prior <- apply(delta.hat, 1, aprior) # FIXME
  # b.prior <- apply(delta.hat, 1, bprior) # FIXME

  a.prior <- apply(delta.hat, 1, aprior.na)
  b.prior <- apply(delta.hat, 1, bprior.na)

  ## Plot empirical and parametric priors.
  if (prior.plots && par.prior) {
    old_pars <- par(no.readonly = TRUE)
    on.exit(par(old_pars))
    par(mfrow=c(2,2))

    ## Top left
    tmp <- density(gamma.hat[1,], na.rm =  T)
    plot(tmp,  type='l', main=expression(paste("Density Plot of First Batch ",  hat(gamma))))
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)

    ## Top Right
    qqnorm(gamma.hat[1,], main=expression(paste("Normal Q-Q Plot of First Batch ", hat(gamma))))
    qqline(gamma.hat[1,], col=2)

    ## Bottom Left
    tmp <- density(delta.hat[1,], na.rm = T)
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    tmp1 <- list(x=xx, y=dinvgamma(xx, a.prior[1], b.prior[1]))
    plot(tmp, typ="l", ylim=c(0, max(tmp$y, tmp1$y)),
         main=expression(paste("Density Plot of First Batch ", hat(delta))))
    lines(tmp1, col=2)

    ## Bottom Right
    invgam <- 1/qgamma(ppoints(ncol(delta.hat)), a.prior[1], b.prior[1])
    qqplot(invgam, delta.hat[1,],
           main=expression(paste("Inverse Gamma Q-Q Plot of First Batch ", hat(delta))),
           ylab="Sample Quantiles", xlab="Theoretical Quantiles")
    lines(c(0, max(invgam)), c(0, max(invgam)), col=2)
  }


  ## Find EB batch adjustments
  ## Dummy values, needed so that the vector operations don't fail
  na.terms <- is.na(gamma.hat)
  gamma.hat[na.terms] <- 0
  delta.hat[na.terms] <- 1

  gamma.star <- delta.star <- matrix(NA, nrow=n.batch, ncol=nrow(s.data))
  if (par.prior) {
    message("Finding parametric adjustments")
    results <- lapply(1:n.batch, function(i) {
      if (mean.only) {
        gamma.star <- sva:::postmean(gamma.hat[i,], gamma.bar[i], 1, 1, t2[i])
        delta.star <- rep(1, nrow(s.data))
      }
      else {
        temp <- sva:::it.sol(s.data[, batches[[i]]], gamma.hat[i, ],
                             delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i],
                             b.prior[i])
        gamma.star <- temp[1, ]
        delta.star <- temp[2, ]
      }
      list(gamma.star=gamma.star, delta.star=delta.star)
    })
    for (i in 1:n.batch) {
      gamma.star[i,] <- results[[i]]$gamma.star
      delta.star[i,] <- results[[i]]$delta.star
    }
  } else {
    message("Finding nonparametric adjustments")
    results <- lapply(1:n.batch, function(i) {
      if (mean.only) {
        delta.hat[i, ] = 1
      }
      # temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
      #                    gamma.hat[i, ], delta.hat[i, ])

      sdat <- as.matrix(s.data[, batches[[i]]])
      g.cle <- gamma.hat
      g.cle[na.terms] <- NaN
      g.hat <- g.cle[i, ]

      d.cle <- delta.hat
      d.cle[na.terms] <- NaN
      d.hat <- d.cle[i, ]

      g.star <- d.star <- NULL
      r <- nrow(sdat)
      for (j in 1:r) {
        g <- g.hat[-j]
        d <- d.hat[-j]
        x <- sdat[j, !is.na(sdat[j, ])]
        n <- length(x)
        q <- numeric(n) + 1
        dat <- matrix(as.numeric(x), length(g), n, byrow = TRUE)
        resid2 <- (dat - g)^2
        sum2 <- resid2 %*% q
        LH <- 1/(2 * pi * d)^(n/2) * exp(-sum2/(2 * d))
        LH[LH == "NaN"] <- 0
        g[g == "NaN"] <- 0
        d[d == "NaN"] <- 0
        g.star <- c(g.star, sum(g * LH)/sum(LH))
        d.star <- c(d.star, sum(d * LH)/sum(LH))
      }

      temp <- rbind(g.star, d.star)
      rownames(temp) <- c("g.star", "d.star")

      list(gamma.star=temp[1,], delta.star=temp[2,])
    })
    for (i in 1:n.batch) {
      gamma.star[i,] <- results[[i]]$gamma.star
      delta.star[i,] <- results[[i]]$delta.star
    }
  }

  if(!is.null(ref.batch)){
    gamma.star[ref,] <- 0  ## set reference batch mean equal to 0
    delta.star[ref,] <- 1  ## set reference batch variance equal to 1
  }

  ## Normalize the Data ###
  message("Adjusting the Data\n")

  bayesdata <- s.data
  j <- 1
  for (i in batches){
    bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j]))) # FIXME
    j <- j+1
  }

  bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean # FIXME

  ## Do not change ref batch at all in reference version
  if(!is.null(ref.batch)){
    bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]
  }

  ## put genes with 0 variance in any batch back in data
  if (length(zero.rows) > 0) {
    dat.orig[keep.rows, ] <- bayesdata
    bayesdata <- dat.orig
  }

  out <- list("corrected data" = bayesdata, "priors" = list("gamma mean" = gamma.bar,
                                                            "gamma variance" = t2,
                                                            "delta shape" = a.prior,
                                                            "delta rate" = b.prior))
  return(out)
}
