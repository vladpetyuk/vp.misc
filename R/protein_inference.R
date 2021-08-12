#' Protein Inference
#'
#' A function for inferring parsimonious peptide to protein mapping.
#'
#'
#' @param x table (data.frame, data.frame or tbl_df) with one column correspoding to peptides
#'          and the other to proteins or in general case accessions
#' @param Accession character corresponding to column name with protein accesssions
#' @param PepSeq character corresponding to column name with peptide sequences
#' @return data.frame with two additional columns
#'      \describe{
#'          \item{\code{UniqueAndRazor}}{TRUE if raw correspond to unique or razor peptide to protein mapping}
#'          \item{\code{JustifiedAccession}}{TRUE if protein has any unique peptide that justify its presence}
#'      }
#' @importFrom data.table setDT rbindlist
#' @importFrom dplyr rename full_join mutate case_when
#' @importFrom magrittr %>%
#' @importFrom rlang !! :=
#'
#' @export protein_inference
#'
#' @examples
#' library(tibble)
#' x <- tribble(
#'     ~pep, ~prot,
#'     "AAAAK", "P1",
#'     "CCCCK", "P1",
#'     "LLLLK", "P1",
#'     "AAAAK", "P2",
#'     "DDDDK", "P2",
#'     "AAAAK", "P3")
#' protein_inference(x, Accession = "prot", PepSeq = "pep")


protein_inference <- function(x, Accession = "Accession", PepSeq = "PepSeq"){

    infer_parsimonious_set <- function(x){
        res <- list()
        setDT(x)
        while(nrow(x) > 0){
            top_prot <- x[, .N, by=Accession][which.max(N),,]$Accession
            top_peps <- subset(x, Accession == top_prot)
            res <- c(res, list(top_peps))
            x <- subset(x, !(PepSeq %in% top_peps$PepSeq))
        }
        return(rbindlist(res, use.names=F, fill=FALSE, idcol=NULL))
    }


    x <- x %>%
        rename(PepSeq = !!PepSeq) %>%
        rename(Accession = !!Accession)

    # main execution
    res <- infer_parsimonious_set(x)

    res$UniqueAndRazor <- TRUE
    x <- x %>%
        full_join(res, by = c("PepSeq", "Accession")) %>%
        mutate(UniqueAndRazor = case_when(is.na(UniqueAndRazor) ~ FALSE,
                                          TRUE ~ UniqueAndRazor))

    # justified protein
    x <- x %>%
        mutate(JustifiedAccession = Accession %in% unique(res$Accession))

    # switch back peptide and protein column names
    x <- x %>%
        rename(!!PepSeq := PepSeq) %>%
        rename(!!Accession := Accession)
    return(x)
}

utils::globalVariables(c(".N", "N"))

