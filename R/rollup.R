
#' Rolling-Up Features Abundance Measurements
#' 
#' A primary purpose of this function is to roll-up peptide
#' measurements to protein level. Ultimately it can be used
#' for rolling-up or combining any features according to a grouping
#' factor.
#' 
#' @param msnset msnset (or most likely eset subclass) object
#' @param rollBy character. A column name in pData(msnset)
#' @param rollFun function for roll-up. "-" for log-transformed data, 
#'                      "/" for normal scale
#' @param verbose logical. Controls real-time output.
#' 
#' @return MSnSet object
#' @importFrom Biobase exprs fData
#' @importFrom outliers grubbs.test outlier
#' @export rrollup
#' 
#' @examples
#' data(srm_msnset)
#' dim(msnset)
#' msnset2 <- rrollup(msnset, rollBy="Gene", rollFun="-", verbose=TRUE)
#' dim(msnset2)


rrollup <- function(msnset, rollBy, rollFun, verbose=TRUE){
    summarisedFeatures <- by(exprs(msnset), fData(msnset)[[rollBy]], 
                             rrollup_a_feature_set, rollFun, verbose)
    exprs.new <- do.call(rbind, as.list(summarisedFeatures))
    if(nrow(unique(fData(msnset))) == nrow(exprs.new)){
        fData.new <- unique(fData(msnset))
    }else{
        fData.new <- data.frame(rollBy=rownames(exprs.new), 
                                stringsAsFactors = FALSE)
        colnames(fData.new) <- rollBy
    }
    fData(msnset) <- fData.new
    exprs(msnset) <- exprs.new
    featureNames(msnset) <- rownames(exprs.new)
    return(msnset)
}


rrollup_a_feature_set <- function(mat, rollFun, verbose){
    if(verbose)
        print(rownames(mat))
    refIdx <- which.max(rowSums(!is.na(mat)))
    mat.ratios <- sweep(mat, 2, as.numeric(mat[refIdx,]), FUN=rollFun)
    scale.factors <- apply(mat.ratios, 1, median, na.rm=TRUE)
    mat <- sweep(mat, 1, scale.factors, FUN=rollFun)
    #
    # for each sample, check if any of the feature (peptide) 
    # abundance measurements
    # representing the set (protein) are not in line with other.
    maxVals <- apply(mat, 1, max, na.rm=TRUE)
    if(nrow(mat) > 2){
        while(grubbs.test(maxVals)$p.value < 0.05 &
              sum(!is.na(maxVals)) > 2){
            i <- which(maxVals == outlier(maxVals))
            maxVals[i] <- NA
        }
    }
    # protein abundance estimates
    mat <- mat[!is.na(maxVals),]
    protProfile <- apply(mat, 2, median, na.rm=TRUE)
    return(protProfile)
}







