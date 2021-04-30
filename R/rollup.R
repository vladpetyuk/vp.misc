
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
                             rrollup_a_feature_set, rollFun, verbose,
                             simplify = FALSE)
    exprs.new <- do.call(rbind, as.list(summarisedFeatures))
    if(nrow(unique(fData(msnset))) == nrow(exprs.new)){
        fData.new <- unique(fData(msnset))
    }else{
        fData.new <- data.frame(rollBy=rownames(exprs.new), 
                                stringsAsFactors = FALSE)
        colnames(fData.new) <- rollBy
        rownames(fData.new) <- fData.new[,rollBy]
    }
    msnset.new <- MSnSet(exprs = as.matrix(exprs.new), 
                         fData = fData.new, 
                         pData = pData(msnset))
    return(msnset.new)
}



rrollup_a_feature_set <- function(mat, rollFun, verbose){
    if (nrow(mat) > 1) {
        # Calculate reference index
        refIdx <- which.max(rowSums(!is.na(mat)))
        # Count number of times each row overlaps with reference
        overlap_count <- 
            apply(mat, 1,
                  function(x) 
                      sum(!is.na(as.numeric(x) == as.numeric(mat[refIdx, ]))))
        overlap_count <- as.numeric(overlap_count)
        # Subset to rows that overlap more than once with reference
        mat <- mat[overlap_count != 1, , drop = F]
        
        if (any(overlap_count == 1) & verbose) {
            warning("Some peptides associated with this protein will be removed.")
        }
        
        # Recalculate index using subset mat
        refIdx <- which.max(rowSums(!is.na(mat)))
        
        mat.ratios <- sweep(mat, 2, as.numeric(mat[refIdx,]), FUN = rollFun)
        scale.factors <- apply(mat.ratios, 1, median, na.rm = TRUE)
        mat <- sweep(mat, 1, scale.factors, FUN = rollFun)
        #
        # for each sample, check if any of the feature (peptide) 
        # abundance measurements
        # representing the set (protein) are not in line with other.
        maxVals <- apply(mat, 1, max, na.rm = TRUE)
        # If all values in a row are NA, the max is -Inf
        # Change infinite values to NA
        maxVals[is.infinite(maxVals)] <- NA
        
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
    } else {
        protProfile <- mat
    }
    
    return(protProfile)
}







