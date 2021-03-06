# 
# 
# 
# 
# select_features_lasso <- function(x, y, ...){
#     cv.glmmod <- cv.glmnet(x=as.matrix(x), 
#                            y=y, 
#                            alpha=1, 
#                            nfolds=length(y),
#                            family="binomial")
#     lambda.sel <- cv.glmmod$lambda.min
#     if(lambda.sel == cv.glmmod$lambda[1])
#         lambda.sel <- cv.glmmod$lambda[2]
#     features <- names(which(coef(cv.glmmod, s=lambda.sel)[-1,1] != 0))
#     return(features)
# }
# 
# 
# train_model_lr <- function(x, y, ...){
#     mod.valid <- glm(subject.type ~ .,
#                      data=data.frame(x,subject.type=y), 
#                      family=binomial(link="logit")) 
#     return(mod.valid)
# }
# 
# 
# 
# 
# lr_cv_round <- function( i, dSet, features, response, select){
#     # feature selection first.
#     if(select){
#         features.sel <- select_features_lasso(x=dSet[-i,features],
#                                               y=dSet[-i,response])
#     }else{
#         features.sel <- features
#     }
#     # train model
#     mdl <- train_model_lr(x=dSet[-i,features.sel,drop=F], 
#                           y=dSet[-i,response])
#     # predict
#     newdata <- dSet[i,features.sel,drop=F]
#     colnames(newdata) <- make.names(colnames(newdata))
#     #
#     list(predict(mdl, newdata=newdata, type='response'), 
#          features.sel)
# }
# 
# 
# 
# cv_round <- function(FUN, ...){
#     FUN(...)
# }
# 
# 
# 
# loocv_mclapply <- function(response, dSet, features, select){
#     res <- mclapply(1:nrow(dSet), function(i){
#                         cv_round(lr_cv_round, i, dSet=dSet, features=features, 
#                                  response=response, select=select)
#                         }, 
#                   mc.cores=detectCores())
#     predProb <- sapply(res, '[[', 1)
#     names(predProb) <- NULL # for compatibility
#     selected <- sapply(res, '[[', 2)
#     list(predProb, selected)
# }
# 
# 
# loocv_foreach <- function(response, dSet, features, select){
#     #
#     # nCores <- ifelse(is.null(nCores), detectCores(), nCores)
#     nCores <- detectCores()
#     cl <- makeCluster(nCores) #, outfile=ifelse(verbose, '', NULL))
#     on.exit(stopCluster(cl))
#     registerDoParallel(cl)
#     #
#     n <- nrow(dSet)
#     f <- function(){
#         pb <- txtProgressBar(min=1, max=n-1,style=3)
#         count <- 0
#         function(...) {
#             count <<- count + length(list(...)) - 1
#             setTxtProgressBar(pb,count)
#             Sys.sleep(0.01)
#             flush.console()
#             c(...)
#         }
#     }
#     res <- foreach( i=icount(n), .packages=NULL, .combine=f()) %dopar%
#         cv_round(lr_cv_round, i, dSet, features, response, select)
#     #
#     predProb <- as.numeric(res[(1:nrow(dSet))*2-1])
#     selected <- res[(1:nrow(dSet))*2]
#     list(predProb, selected)
# }
# 
# 
# loocv_onethread <- function(response, dSet, features, select){
#     predProb <- numeric(nrow(dSet))
#     selected <- vector("list",nrow(dSet))
#     #---
#     pb <- txtProgressBar(min=1, max=nrow(dSet)-1, style=3)
#     for( i in 1:nrow(dSet)){
#         res <- cv_round(lr_cv_round, i, dSet, features, response, select)
#         predProb[i] <- res[[1]]
#         selected[[i]] <- res[[2]]
#         setTxtProgressBar(pb, i)
#     }
#     close(pb)
#     list(predProb, selected)
# }
# 
# 
# #' Logistic Regression Predictive Models
# #' 
# #' The objective is - given the number of features,
# #' select the most informative ones and evaluate
# #' the predictive logistic regression model. 
# #' The feature and model selection performed independently
# #' for each round of LOOCV. Feature selection performed using
# #' LASSO approach.
# #' 
# # ' @param response character name of the predicted variable
# # '          it has to have two values and one of them should be "case"
# #' @param msnset MSnSet object
# #' @param features character vector features to select from for 
# #'          building prediction model. The features can be either in
# #'          featureNames(msnset) or in pData(msnset).
# #' @param predCls class to predict
# #' @param select logical to select features using LASSO
# #'           or use the entire set?
# #' @param par.backend type of backend to support parallelizattion. 
# #'          'mc' uses mclapply from parallel,
# #'          'foreach' is based on 'foreach',
# #'          'none' - just a single thread.
# # ' @param verbose 0 - silent, 1 - count round, 
# # '          2 - print selected features at each round
# # ' @param nCores number of cores for parallel processing
# #' @return list 
# #'      \describe{
# #'          \item{\code{prob}}{is the probabilities (response) from LOOCV
# #'          that the sample is "case"}
# #'          \item{\code{features}}{list of selected features
# #'          for each iteration of LOOCV}
# #'          \item{\code{top}}{top features over all iterations}
# #'          \item{\code{auc}}{AUC}
# #'          \item{\code{pred}}{prediction perfomance obtained by
# #'                  \code{ROCR::prediction}}
# #'      }
# #' @import glmnet
# #' @importFrom ROCR prediction performance
# #' @importFrom parallel makeCluster stopCluster detectCores mclapply
# #' @importFrom iterators icount
# #' @importFrom doParallel registerDoParallel
# #' @importFrom foreach foreach %dopar%
# #' @export lr_modeling
# #' 
# lr_modeling <- function( msnset, features, predCls="case", select=T,
#                          par.backend=c('mc','foreach','none')){
#     # features - name vector. What if I want to add stuff from pData? G'head.
#     # msnset - MSnSet (eSet) object
#     # select - to select features or use the entire set?
#     # verbose - 0,1,2 from silent to printing rounds and features
#     #---
#     par.backend <- match.arg(par.backend)
#     # prepare data
#     response <- "subject.type" # in this project it is always subject.type
#     dSet <- cbind(pData(msnset), t(exprs(msnset)))
#     #
#     stopifnot(length(unique(dSet[,response])) == 2)
#     stopifnot(predCls %in% dSet[,response] )
#     if(unique(dSet[,response])[1] == predCls){
#         lvlz <- rev(unique(dSet[,response]))
#     }else{
#         lvlz <- unique(dSet[,response])
#     }
#     dSet[,response] <- factor(dSet[,response], levels=lvlz)
#     #
#     res <- switch(par.backend,
#                   mc = loocv_mclapply(response, dSet, features, select),
#                   foreach = loocv_foreach(response, dSet, features, select),
#                   none = loocv_onethread(response, dSet, features, select))
#     pred <- prediction(res[[1]], dSet[,response] == predCls)
#     return(list(prob = res[[1]], 
#                 features = res[[2]], 
#                 top = rev(sort(table(unlist(res[[2]])))),
#                 auc = performance(pred,"auc")@y.values[[1]],
#                 pred = pred))
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # # ' @export lr_modeling.2
# # # ' 
# # lr_modeling.2 <- function( msnset, features, predCls="case", 
# #                          select=T, verbose=1, nCores=NULL){
# #     # features - name vector. What if I want to add stuff from pData? G'head.
# #     # msnset - MSnSet (eSet) object
# #     # select - to select features or use the entire set?
# #     # verbose - 0,1,2 from silent to printing rounds and features
# #     #---
# #     # prepare data
# #     response <- "subject.type" # in this project it is always subject.type
# #     dSet <- cbind(pData(msnset), t(exprs(msnset)))
# #     #
# #     stopifnot(length(unique(dSet[,response])) == 2)
# #     stopifnot(predCls %in% dSet[,response] )
# #     if(unique(dSet[,response])[1] == predCls){
# #         lvlz <- rev(unique(dSet[,response]))
# #     }else{
# #         lvlz <- unique(dSet[,response])
# #     }
# #     dSet[,response] <- factor(dSet[,response], levels=lvlz)
# #     #############################
# #     nCores <- ifelse(is.null(nCores), detectCores(), nCores)
# #     cl <- makeCluster(nCores) #, outfile=ifelse(verbose, '', NULL))
# #     on.exit(stopCluster(cl))
# #     registerDoParallel(cl)
# #     #Progress combine function
# #     n <- nrow(dSet)
# #     f <- function(){
# #         pb <- txtProgressBar(min=1, max=n-1,style=3)
# #         count <- 0
# #         function(...) {
# #             count <<- count + length(list(...)) - 1
# #             setTxtProgressBar(pb,count)
# #             Sys.sleep(0.01)
# #             flush.console()
# #             c(...)
# #         }
# #     }
# #     res <- foreach( i=icount(n), .packages=NULL, .combine=f()) %dopar%
# #     {
# #         # feature selection first.
# #         if(select){
# #             features.sel <- suppressWarnings(
# #                 select_features_lasso(x=dSet[-i,features],
# #                                       y=dSet[-i,response])
# #             )
# #         }else{
# #             features.sel <- features
# #         }
# #         # train model
# #         mdl <- train_model_lr(x=dSet[-i,features.sel,drop=F], 
# #                               y=dSet[-i,response])
# #         # predict
# #         newdata <- dSet[i,features.sel,drop=F]
# #         colnames(newdata) <- make.names(colnames(newdata))
# #         #
# #         list(predict(mdl, newdata=newdata, type='response'), 
# #              features.sel)
# #     }
# #     #
# #     predProb <- as.numeric(res[(1:nrow(dSet))*2-1])
# #     selected <- res[(1:nrow(dSet))*2]
# #     #############################
# #     pred <- prediction(predProb, dSet[,response] == predCls)
# #     auc <- performance(pred,"auc")@y.values[[1]]
# #     top <- rev(sort(table(unlist(selected))))
# #     #
# #     return(list(prob=predProb, 
# #                 features=selected, 
# #                 top=top,
# #                 auc=auc,
# #                 pred=pred))
# # }
# 
# 
# 
# #' @export lr_modeling.v1
# #' 
# lr_modeling.v1 <- function( msnset, features, predCls="case", select=T, verbose=1){
#     # features - name vector. What if I want to add stuff from pData? G'head.
#     # msnset - MSnSet (eSet) object
#     # select - to select features or use the entire set?
#     # verbose - 0,1,2 from silent to printing rounds and features
#     #---
#     # prepare data
#     response <- "subject.type" # in this project it is always subject.type
#     dSet <- cbind(pData(msnset), t(exprs(msnset)))
#     #
#     stopifnot(length(unique(dSet[,response])) == 2)
#     stopifnot(predCls %in% dSet[,response] )
#     if(unique(dSet[,response])[1] == predCls){
#         lvlz <- rev(unique(dSet[,response]))
#     }else{
#         lvlz <- unique(dSet[,response])
#     }
#     dSet[,response] <- factor(dSet[,response], levels=lvlz)
#     #############################
#     predProb <- numeric(nrow(dSet))
#     selected <- vector("list",nrow(dSet))
#     #---
#     for( i in 1:nrow(dSet)){
#         # feature selection first.
#         if(select){
#             features.sel <- select_features_lasso(x=dSet[-i,features],
#                                                   y=dSet[-i,response])
#         }else{
#             features.sel <- features
#         }
#         selected[[i]] <- features.sel
#         # train model
#         mdl <- train_model_lr(x=dSet[-i,features.sel,drop=F], 
#                               y=dSet[-i,response])
#         # predict
#         newdata <- dSet[i,features.sel,drop=F]
#         colnames(newdata) <- make.names(colnames(newdata))
#         predProb[i] <- predict(mdl, newdata=newdata, type='response')
#         #
#         if(verbose > 0) print(i)
#         if(verbose > 1) print(features.sel)
#     }
#     #
#     pred <- prediction(predProb, dSet[,response] == predCls)
#     auc <- performance(pred,"auc")@y.values[[1]]
#     #############################
#     #
#     top <- rev(sort(table(unlist(selected))))
#     #
#     return(list(prob=predProb, 
#                 features=selected, 
#                 top=top,
#                 auc=auc,
#                 pred=pred))
# }
# 
# 
