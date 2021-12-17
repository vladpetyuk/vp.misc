#' Normalization of LC-MS/MS Data
#'
#' Fits a nonparametric trend into relative
#' abundance data. Any significant (at this point just any)
#' trend is corrected.
#' Converts one ExpressionSet/MSnSet to another ExpressionSet/MSnSet.
#'
#' @param eset ExpressionSet/MSnSet object
#' @param property character the column in the fData that the
#'              relative intensities regress against.
#' @param method character. either "lowess" or "loess" at this point.
#' @param partition_by character. Column name in \code{fData}, according to
#'                                which features will be partitioned. This is
#'                                when, for example, normalizing by LC
#'                                and there are multiple fractions per sample.
#' @param ... passed to \code{lowess} or \code{loess}. We suggest to use span = 0.5 with loess.
#'
#' @note So far the only property I have in mind is elution time in
#'      label-free LC-MS/MS data.
#'
#' @importFrom Biobase exprs fData sampleNames exprs<-
#' @importFrom stats lowess loess predict loess.control
#' @importFrom graphics points grid
#'
#' @export normalize_by_feature_property
#' @export normalize_by_feature_property_partitioned
#' @export evaluate_parameter_effect
#'


normalize_by_feature_property <- function(eset,
                                          property,
                                          method=c("lowess","loess"),
                                          ...){
    # elution time is a typical problem in a label-free LC-MS/MS
    method <- match.arg(method)
    o <- order(fData(eset)[[property]])
    if(method == "lowess"){
        for(sample in sampleNames(eset)){
            m <- lowess(fData(eset)[[property]], exprs(eset)[,sample], ...)
            exprs(eset)[,sample] <- exprs(eset)[,sample] - m$y[order(o)]
        }
    }else if(method == "loess"){
        x <- cbind(fData(eset)[,property,drop=F], exprs(eset))
        for(sample in sampleNames(eset)){
            xx <- x[,c(property,sample)]
            colnames(xx) <- c('x','y')
            m <- loess(y ~ x, data=xx, ...)
            exprs(eset)[,sample] <- exprs(eset)[,sample] - predict(m, xx)
        }
    }else{
        warning("method is not defined")
    }
    return(eset)
}


#' @rdname normalize_by_feature_property
normalize_by_feature_property_partitioned <- function(eset,
                                                      property,
                                                      partition_by,
                                                      method=c("lowess","loess"),
                                                      ...){
  partition_vec <- fData(eset)[[partition_by]]
  for(part_i in unique(partition_vec)){
    # subset the object
    idx <- partition_vec == part_i
    eset_part_i <- eset[idx,]
    # correct the ESI effect
    eset_part_i <- normalize_by_feature_property(eset_part_i,
                                                 property,
                                                 method,
                                                 ...)
    # update the main object
    exprs(eset)[idx,] <- exprs(eset_part_i)
  }
  return(eset)
}



#' @rdname normalize_by_feature_property
evaluate_parameter_effect <- function(eset,
                                      sample,
                                      property,
                                      partition_by,
                                      partition_value,
                                      ...){
  #
  eset <- eset[fData(eset)[[partition_by]] == partition_value,]
  model_i <- loess(exprs(eset)[,sample] ~ fData(eset)[[property]],
                   # span = 0.5, # hardcoded so far
                   ...,
                   control=loess.control(surface="direct"))
  bias_i <- predict(model_i, fData(eset)[[property]])
  plot(fData(eset)[[property]], exprs(eset)[,sample])
  grid()
  points(fData(eset)[[property]], bias_i, col="red")
}
























