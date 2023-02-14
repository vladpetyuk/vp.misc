#'
#' @importMethodsFrom MSnID map_peptide_position
#' @importFrom MSnID MSnID `psms<-`
#' @importFrom Biostrings readAAStringSet
#'
NULL

#' @export
#'
#' @docType methods


setMethod("map_peptide_position", "MSnSet",
          definition=function(object,
                              fasta,
                              accession_col = "accession",
                              peptide_col = "peptide",
                              ...)
          {
            object <- .map_peptide_position(object, fasta, accession_col, peptide_col, ...)
            return(object)
          }
)


.map_peptide_position <- function(m, fasta, accession_col, peptide_col, ...){

   ms <- MSnID(".")
   psms(ms) <- fData(m) %>%
      select(!!accession_col, !!peptide_col) %>%
      dplyr::rename(peptide = !!peptide_col) %>%
      dplyr::rename(accession = !!accession_col) %>%
      mutate(isDecoy = NA,
             calculatedMassToCharge = NA,
             experimentalMassToCharge = NA,
             chargeState = NA,
             spectrumFile = NA,
             spectrumID = NA)

    fa <- readAAStringSet(fasta)
    names(fa) <- sub("([^ ]+) .*", "\\1", names(fa))
    x <- MSnID::map_peptide_position(ms, fa)
    xsel <- psms(x) %>%
       select(accession, cleanSeq, First_AA_First, Last_AA_First, ProtLen) %>%
       dplyr::rename(First_AA = First_AA_First,
                  Last_AA = Last_AA_First) %>%
   as.data.frame()
   rownames(xsel) <- xsel$cleanSeq

   safe_features <- intersect(rownames(xsel), featureNames(m)) # for safety
   xsel <- xsel[safe_features,]
   m <- m[safe_features,]
   fData(m) <- cbind(fData(m), select(xsel, First_AA, Last_AA, ProtLen))

   return(m)

  return(NULL)

}






