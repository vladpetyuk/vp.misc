#' Plot AUC
#' 
#' Plot AUC after LOOCV Model Evaluation
#' 
#' @param modelingResult output from either 
#'          \code{lr_modeling} or \code{rf_modeling}
#' 
#' @importFrom ROCR performance
#' 
#' @export plotAUC
#' 
plotAUC <- function(modelingResult){
    perf <- performance(modelingResult$pred,"tpr","fpr")
    # x=1-spec, y=sens
    plot(perf, 
         main=sprintf("AUC: %s", round(modelingResult$auc,2)), 
         col=2, lwd=2) 
    abline(a=0,b=1,lwd=2,lty=2,col="gray")
}









