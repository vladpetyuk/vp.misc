



select_features_lasso <- function(x, y, ...){
    cv.glmmod <- cv.glmnet(x=as.matrix(x), 
                           y=y, 
                           alpha=1, 
                           nfolds=length(y),
                           family="binomial")
    lambda.sel <- cv.glmmod$lambda.min
    if(lambda.sel == cv.glmmod$lambda[1])
        lambda.sel <- cv.glmmod$lambda[2]
    features <- names(which(coef(cv.glmmod, s=lambda.sel)[-1,1] != 0))
    return(features)
}


train_model_lr <- function(x, y, ...){
    mod.valid <- glm(subject.type ~ .,
                     data=data.frame(x,subject.type=y), 
                     family=binomial(link="logit")) 
    return(mod.valid)
}




lr_cv_round <- function( i, dSet, features, response, select){
    # feature selection first.
    if(select){
        features.sel <- select_features_lasso(x=dSet[-i,features],
                                              y=dSet[-i,response])
    }else{
        features.sel <- features
    }
    # train model
    mdl <- train_model_lr(x=dSet[-i,features.sel,drop=F], 
                          y=dSet[-i,response])
    # predict
    newdata <- dSet[i,features.sel,drop=F]
    colnames(newdata) <- make.names(colnames(newdata))
    #
    list(predict(mdl, newdata=newdata, type='response'), 
         features.sel)
}



cv_round <- function(FUN, ...){
    FUN(...)
}



loocv_mclapply <- function(response, dSet, features, select){
    res <- mclapply(1:nrow(dSet), function(i){
                        cv_round(lr_cv_round, i, dSet=dSet, features=features, 
                                 response=response, select=select)
                        }, 
                  mc.cores=detectCores())
    predProb <- sapply(res, '[[', 1)
    names(predProb) <- NULL # for compatibility
    selected <- sapply(res, '[[', 2)
    list(predProb, selected)
}


loocv_foreach <- function(response, dSet, features, select){
    #
    # nCores <- ifelse(is.null(nCores), detectCores(), nCores)
    nCores <- detectCores()
    cl <- makeCluster(nCores) #, outfile=ifelse(verbose, '', NULL))
    on.exit(stopCluster(cl))
    registerDoParallel(cl)
    #
    n <- nrow(dSet)
    f <- function(){
        pb <- txtProgressBar(min=1, max=n-1,style=3)
        count <- 0
        function(...) {
            count <<- count + length(list(...)) - 1
            setTxtProgressBar(pb,count)
            Sys.sleep(0.01)
            flush.console()
            c(...)
        }
    }
    res <- foreach( i=icount(n), .packages=NULL, .combine=f()) %dopar%
        cv_round(lr_cv_round, i, dSet, features, response, select)
    #
    predProb <- as.numeric(res[(1:nrow(dSet))*2-1])
    selected <- res[(1:nrow(dSet))*2]
    list(predProb, selected)
}


loocv_onethread <- function(response, dSet, features, select){
    predProb <- numeric(nrow(dSet))
    selected <- vector("list",nrow(dSet))
    #---
    pb <- txtProgressBar(min=1, max=nrow(dSet)-1, style=3)
    for( i in 1:nrow(dSet)){
        res <- cv_round(lr_cv_round, i, dSet, features, response, select)
        predProb[i] <- res[[1]]
        selected[[i]] <- res[[2]]
        setTxtProgressBar(pb, i)
    }
    close(pb)
    list(predProb, selected)
}


#' Logistic Regression Predictive Models
#' 
#' The objective is - given the number of features,
#' select the most informative ones and evaluate
#' the predictive logistic regression model. 
#' The feature and model selection performed independently
#' for each round of LOOCV. Feature selection performed using
#' LASSO approach.
#' 
#' @param msnset MSnSet object. Note - can it be generalized to eset?
#' @param features character vector features to select from for 
#'          building prediction model. The features can be either in
#'          featureNames(msnset) or in pData(msnset).
#' @param response factor to classify along. Must be only 2 levels.
#' @param pred.cls character, class to predict
#' @param sel.feat logical to select features using LASSO
#'           or use the entire set?
#' @param par.backend type of backend to support parallelizattion. 
#'          'mc' uses mclapply from parallel,
#'          'foreach' is based on 'foreach',
#'          'none' - just a single thread.
#' @return list 
#'      \describe{
#'          \item{\code{prob}}{is the probabilities (response) from LOOCV
#'          that the sample is "case".  That is how well model trained 
#'          on other samples, predicts this particular one.}
#'          \item{\code{features}}{list of selected features
#'          for each iteration of LOOCV}
#'          \item{\code{top}}{top features over all iterations}
#'          \item{\code{auc}}{AUC}
#'          \item{\code{pred}}{prediction perfomance obtained by
#'                  \code{ROCR::prediction}}
#'      }
#' @import glmnet
#' @importFrom ROCR prediction performance
#' @importFrom parallel makeCluster stopCluster detectCores mclapply
#' @importFrom iterators icount
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' 
#' @export lr_modeling
#' 
#' @examples
#' data(srm_msnset)
#' head(varLabels(msnset))
#' head(msnset$subject.type)
#' # reduce to two classes
#' msnset <- msnset[,msnset$subject.type != "control.1"]
#' msnset$subject.type <- as.factor(msnset$subject.type)
#' # Note, par.backend="none" is for the example only.
#' out <- lr_modeling(msnset, 
#'                    features=featureNames(msnset), 
#'                    response="subject.type", 
#'                    pred.cls="case", par.backend="none")
#' plotAUC(out)
#' # top features consistently recurring in the models during LOOCV
#' print(out$top)
#' # the AUC
#' print(out$auc)
#' # probabilities of classifying the sample right, if the feature selection
#' # and model training was performed on other samples
#' plot(sort(out$prob))
#' abline(h=0.5, col='red')
#' 
#' 
#' 
lr_modeling <- function( msnset, features, response, pred.cls, sel.feat=T,
                         par.backend=c('mc','foreach','none')){
    # features - name vector. What if I want to add stuff from pData? G'head.
    # msnset - MSnSet (eSet) object
    # sel.feat - to select features or use the entire set?
    # verbose - 0,1,2 from silent to printing rounds and features
    #---
    par.backend <- match.arg(par.backend)
    if(.Platform$OS.type != 'unix' & par.backend == 'mc')
        par.backend <- 'foreach'
    # prepare data
    dSet <- cbind(pData(msnset), t(exprs(msnset)))
    #
    stopifnot(length(unique(dSet[,response])) == 2)
    stopifnot(pred.cls %in% dSet[,response] )
    if(unique(dSet[,response])[1] == pred.cls){
        lvlz <- rev(unique(dSet[,response]))
    }else{
        lvlz <- unique(dSet[,response])
    }
    dSet[,response] <- factor(dSet[,response], levels=lvlz)
    #
    res <- switch(par.backend,
                  mc = loocv_mclapply(response, dSet, features, sel.feat),
                  foreach = loocv_foreach(response, dSet, features, sel.feat),
                  none = loocv_onethread(response, dSet, features, sel.feat))
    pred <- prediction(res[[1]], dSet[,response] == pred.cls)
    return(list(prob = res[[1]], 
                features = res[[2]], 
                top = rev(sort(table(unlist(res[[2]])))),
                auc = performance(pred,"auc")@y.values[[1]],
                pred = pred))
}

