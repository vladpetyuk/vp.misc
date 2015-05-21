
#' @name msnset
#' @title SRM proteomics data from human CSF
#' @description Targeted SRM analysis of 33 proteins in CSF
#' @docType data
# ' @usage CpG.human.GRCh37
#' @format a \code{MSnSet} instance
# ' @source UCSC Table Browser
#' @author Vlad Petyuk, 2014-03-24
# ' @export
NULL


# #' Prices of 50,000 round cut diamonds.
# #'
# #' A dataset containing the prices and other attributes of almost 54,000
# #' diamonds. The variables are as follows:
# #'
# #' \itemize{
# #'   \item price. price in US dollars (\$326--\$18,823)
# #'   \item carat. weight of the diamond (0.2--5.01)
# #'   ...
# #' }
# #'
# #' @format A data frame with 53940 rows and 10 variables
# #' @source \url{http://www.diamondse.info/}
# #' @name diamonds
# NULL


#' CPTAC Ovarian Cancer Proteomics Dataset
#' 
#' Processed iTRAQ data from both PNNL and JHU sites.
#' 
#' @format an \code{MSnID} instance.  8103 proteins by 73 features.
#' @source CPTAC study
#' @name oca.set
#' @examples
#' data(cptac_oca)
#' head(pData(oca.set))
#' head(fData(oca.set))
#' head(exprs(oca.set))
NULL