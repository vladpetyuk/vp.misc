# 
# 
# 
# #' @import varSelRF
# select_features_varSelRF <- function(x, y, ...){
#     rfvs <- varSelRF(xdata = x, Class = y,
#                      verbose=FALSE, 
#                      recompute.var.imp=FALSE, # tweak
#                      whole.range = TRUE,
#                      ntree=10000, # tweak
#                      ntreeIterat=10000, # tweak
#                      vars.drop.num=NULL, 
#                      vars.drop.frac=0.2, 
#                      ...)
#     hi <- rfvs$selec.history
#     idx <- which.min(hi$OOB)
#     idx <- min(idx + 5, nrow(hi)) # start from 5 steps back
#     vrz <- hi[idx, "Vars.in.Forest"] 
#     vrz <- strsplit(as.character(vrz)," \\+ ")[[1]]
#     #
#     rfvs <- varSelRF(xdata = x[,vrz], Class = y,
#                      verbose=FALSE, 
#                      recompute.var.imp=FALSE, # tweak
#                      whole.range = TRUE,
#                      ntree=10000, # tweak
#                      ntreeIterat=10000, # tweak
#                      vars.drop.num=1, 
#                      vars.drop.frac=NULL, 
#                      ...)
#     return(rfvs$selected.vars)
#     
# }
# 
# 
# #' @import Boruta
# select_features_Boruta <- function(x, y, ...){
#     
# #     # .. add on ..
# #     # .. preselection of features using plain random forest
# #     # .. perhaps I can skimp on Boruta settings
# #     #
# #     rf <- randomForest(x, y, data=dSet, ntree = 10000)
# #     f2 <- names(head(importance(rf)[order(importance(rf)[,1],decreasing=T),],20))
# #     x <- x[,f2]
# #     #
#     
#     rfbo <- Boruta(x, y, 
#                    #maxRuns = 100, # 1000
#                    doTrace=0, 
#                    #ntree=5000, # 10000
#                    ...) 
#     # 5h with 10000
#     # 10m with 1000
#     # ~ 5 min with 500
#     # features_i <- getSelectedAttributes(TentativeRoughFix(rfbo))
#     features_i <- getSelectedAttributes(rfbo)
#     if(identical(character(),features_i))
#         features_i <- getSelectedAttributes(rfbo, withTentative = T)
#     if(identical(character(),features_i)){
#         # return just the best one
#         meds <- apply(rfbo$ImpHistory,2,median)
#         meds <- meds[names(rfbo$finalDecision)]
#         features_i <- names(rev(sort(meds)))[1]
#     }
#     return(features_i)
# }
# 
# 
# 
# select_features_top <- function(x, y, ...){
#     r <- randomForest(x, y, ntree=1000, ...)
#     ordr <- order(importance(r), decreasing = T)
#     imp <- importance(r)[ordr,,drop=F]
#     return(rownames(head(imp, 5)))
# }
# 
# 
# 
# train_model_rf <- function(x, y, ...){
#     r <- randomForest(x, y, ntree=1000, ...)
#     return(r)
# }
# 
# 
# 
# #' Random Forest Predictive Models
# #' 
# #' The objective is - given the number of features,
# #' select the most informative ones and evaluate
# #' the predictive random forest model. 
# #' The feature and model selection performed independently
# #' for each round of LOOCV.
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
# #' @param selAgl character.
# #'      \itemize{
# #'          \item varSelRF
# #'          \item Boruta
# #'          \item top - just selecting top 3 features
# #'      }
# #' @param verbose 0 - silent, 1 - count round, 
# #'          2 - print selected features at each round
# #'          
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
# #' @import randomForest
# #' @importFrom ROCR prediction performance
# #' @importFrom parallel mclapply detectCores
# #' @export rf_modeling
# #'
# rf_modeling <- function( msnset, features, predCls="case", select=T,
#                          selAlg=c("varSelRF","Boruta","top"), ...){
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
#     #---
#     FUN <- switch(selAlg, 
#                   varSelRF = select_features_varSelRF,
#                   Boruta = function(x,y) select_features_Boruta(x,y,...),
#                   top = select_features_top)
#     res <- mclapply(1:nrow(dSet), function(i){
#         if(select){
#             features.sel <- FUN(x=dSet[-i,features],
#                                y=dSet[-i,response])
#         }else{
#             features.sel <- features
#         }
#         # train model
#         x=dSet[-i,features.sel,drop=F]
#         colnames(x) <- make.names(colnames(x))
#         mdl <- train_model_rf(x=x, y=dSet[-i,response])
#         # predict
#         newdata <- dSet[i,features.sel,drop=F]
#         colnames(newdata) <- make.names(colnames(newdata))
#         predProb <- predict(mdl, newdata=newdata, type='prob')[1,predCls]
#         print(i)
#         list(predProb, features.sel)
#     }, mc.cores=detectCores(), mc.set.seed = FALSE)
#     #
#     predProb <- sapply(res, '[[', 1)
#     names(predProb) <- NULL # for compatibility
#     selected <- sapply(res, '[[', 2)
#     #---
#     pred <- prediction(predProb, dSet[,response] == predCls)
#     return(list(prob = predProb, 
#                 features = selected, 
#                 top = rev(sort(table(unlist(selected)))),
#                 auc = performance(pred,"auc")@y.values[[1]],
#                 pred = pred))
# }
# 
# 
# 
# 
# #' @export rf_modeling.v1
# rf_modeling.v1 <- function( msnset, features, predCls="case", select=T, 
#                          selAlg=c("varSelRF","Boruta","top"), 
#                          verbose=1){
#     # msnset - MSnSet (eSet) object
#     # features - name vector. What if I want to add stuff from pData? G'head.
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
#     #
#     predProb <- numeric(nrow(dSet))
#     selected <- vector("list",nrow(dSet))
#     #---
#     for( i in 1:nrow(dSet)){
#         # feature selection first.
#         if(select){
#             features.sel <- switch(selAlg,
#                 varSelRF = select_features_varSelRF(x=dSet[-i,features],
#                                                     y=dSet[-i,response]),
#                 Boruta = select_features_Boruta(x=dSet[-i,features],
#                                                   y=dSet[-i,response]),
#                 top = select_features_top(x=dSet[-i,features],
#                                           y=dSet[-i,response]))
#         }else{
#             features.sel <- features
#         }
#         selected[[i]] <- features.sel
#         # train model
#         x=dSet[-i,features.sel,drop=F]
#         colnames(x) <- make.names(colnames(x))
#         mdl <- train_model_rf(x=x, y=dSet[-i,response])
#         # predict
#         newdata <- dSet[i,features.sel,drop=F]
#         colnames(newdata) <- make.names(colnames(newdata))
#         predProb[i] <- predict(mdl, newdata=newdata, type='prob')[1,predCls]
#         #
#         if(verbose > 0) print(i)
#         if(verbose > 1) print(features.sel)
#     }
#     #
#     pred <- prediction(predProb, dSet[,response] == predCls)
#     auc <- performance(pred,"auc")@y.values[[1]]
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
