#' Map positions of peptides on protein sequences
#'
#' @param object \code{MSnSet} object with \code{fData}.
#' @param fasta character; path to FASTA file.
#' @param accession_col character; name of column in \code{fData(object)}
#'   containing accessions (usually proteins).
#' @param peptide_col character; name of column in \code{fData(object)}
#'   containing peptides.
#' @param ... additional arguments
#'
#' @importMethodsFrom MSnID map_peptide_position
#' @importFrom MSnID MSnID psms `psms<-`
#' @importFrom MSnbase `fData<-`
#' @importFrom Biostrings readAAStringSet
#' @importFrom dplyr %>% select rename mutate
#' @importFrom rlang sym !!
#' @importFrom methods setMethod
#'
#' @aliases map_peptide_position,MSnSet-character-method
#'
#' @export
#'
#' @docType methods

setMethod(f = "map_peptide_position",
          signature = c(object = "MSnSet", fasta = "character"),
          definition = function(object,
                                fasta,
                                accession_col = "accession",
                                peptide_col = "peptide",
                                ...)
          {
            object <- .map_peptide_position(object,
                                            fasta,
                                            accession_col,
                                            peptide_col,
                                            ...)
            return(object)
          }
)


.map_peptide_position <- function(m,
                                  fasta,
                                  accession_col,
                                  peptide_col,
                                  ...)
{
  fa <- readAAStringSet(fasta)
  names(fa) <- sub(" .*", "", names(fa))

  ms <- MSnID(".")

  psms(ms) <- fData(m) %>%
    select(!!accession_col, !!peptide_col) %>%
    dplyr::rename(peptide = !!peptide_col,
                  accession = !!accession_col) %>%
    mutate(isDecoy = NA,
           calculatedMassToCharge = NA,
           experimentalMassToCharge = NA,
           chargeState = NA,
           spectrumFile = NA,
           spectrumID = NA)

  ms <- MSnID::map_peptide_position(object = ms, fasta = fa, ...)

  xsel <- psms(ms) %>%
    select(accession, cleanSeq,
           First_AA_First, Last_AA_First,
           ProtLen) %>%
    dplyr::rename(First_AA = First_AA_First,
                  Last_AA = Last_AA_First) %>%
    as.data.frame()
  rownames(xsel) <- xsel$cleanSeq

  safe_features <- intersect(rownames(xsel), featureNames(m)) # for safety
  xsel <- xsel[safe_features, ]
  m <- m[safe_features, ]

  fData(m) <- cbind(fData(m), select(xsel, First_AA, Last_AA, ProtLen))

  return(m)
}



