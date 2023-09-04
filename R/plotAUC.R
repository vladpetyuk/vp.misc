#' Plot AUC
#'
#' Plot AUC after LOOCV Model Evaluation
#'
#' @param modelingResult output from either
#'          \code{lr_modeling} or \code{rf_modeling}
#' @param CI (logical) Plot confidence interval or not. Default is False.
#' @param ... further arguments passes to the \code{plot} method of the
#'             \code{\link[pROC]{ci.se}} object of the \code{pROC} package.
#'
#'
#' @importFrom ROCR performance
#' @importFrom graphics abline
#' @importFrom pROC roc ci.se
#'
#' @export plotAUC
#'
plotAUC <- function(modelingResult, CI = FALSE, ...){

   if(!CI){
      perf <- performance(modelingResult$pred,"tpr","fpr")
       # x=1-spec, y=sens
       plot(perf,
            main=sprintf("AUC: %s", round(modelingResult$auc,2)),
            col=2, lwd=2)
       abline(a=0,b=1,lwd=2,lty=2,col="gray")
   }else{
      old_par <- par()
      par(pty="s")
      pROC_obj <- roc(modelingResult$pred@labels[[1]],
                      modelingResult$pred@predictions[[1]],
                      direction = "<",
                      smoothed = T,
                      # arguments for ci
                      ci=TRUE,
                      ci.alpha=0.95,
                      stratified=FALSE,
                      # arguments for plot
                      plot=TRUE,
                      auc.polygon=F,
                      max.auc.polygon=F,
                      grid=F,
                      print.auc=TRUE,
                      print.auc.y=0.1,
                      print.auc.x=0.8,
                      show.thres=TRUE)
      sens.ci <- ci.se(pROC_obj)
      plot(sens.ci, type="shape", col="#FF8888", conf=95, ...)
      abline(a=1,b=-1,lty=2,lwd=2,col="grey50")
      par(old_par)
   }

}









