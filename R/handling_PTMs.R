#' PTM Mapping Utilities
#'
#' Map PTM locations from peptide to protein. Creates new columns that describe
#' which amino acids are modified and where the modifications occur in the
#' protein sequences.
#'
#' @param ids \code{data.frame} object. Must have column containing protein IDs
#'            that match names in the provided FASTA file. The other column is
#'            peptide IDs containing flanking AAs in the form \code{X.X...X.X}.
#' @param fasta \code{\link[Biostrings]{AAStringSet}} object.
#'        \code{names(fasta)} must match protein IDs in the corresponding column
#'        of the `ids` object.
#' @param prot_id_col character. Name of the column with protein IDs in the
#'        `ids` object.
#' @param peptide_col character. Name of the column with peptides in the `ids`
#'        object.
#' @param mod_char character. Character denoting mapped PTM in peptides.
#'        Typically \code{"*"}.
#'
#' @return
#' \code{data.frame}
#'
#' @importFrom dplyr mutate rename inner_join filter select row_number
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_detect str_replace_all str_locate_all str_length
#'             str_sub
#' @importFrom rlang !! parse_expr
#' @importFrom purrr map map_chr map2
#'
#' @export map_PTM_sites
#'
#' @examples
#' # Path to FASTA file
#' fasta_file_name <- system.file("extdata/FASTAs",
#'                                "rattus_norvegics_uniprot_2018_09.fasta.gz",
#'                                package = "MSnSet.utils")
#' library(Biostrings)
#' # Create AAStringSet object from FASTA file
#' fst <- readAAStringSet(fasta_file_name, format="fasta",
#'                        nrec=-1L, skip=0L, use.names=TRUE)
#'
#' # Extract UniProt Accessions
#' names(fst) <- sub("^.*\\|(.*)\\|.*$","\\1",names(fst))
#'
#' # PTM data
#' data(phospho_identifications_rat)
#'
#' # Get protein PTM locations
#' ids_with_sites <- map_PTM_sites(ids = ids, fasta = fst,
#'                                 prot_id_col = "UniProtAccFull",
#'                                 peptide_col = "Peptide", mod_char = "*")
#'


map_PTM_sites <- function(ids, fasta, prot_id_col, peptide_col, mod_char){

    # check if fasta entry names are unique
    if(any(duplicated(names(fasta)))){
        message("FASTA entry names are not unique!")
        stop()
    }

    # check if there is at least some agreement in IDs
    if(length(intersect(ids[[prot_id_col]], names(fasta))) == 0){
        message("There is zero overlap in protein IDs and FASTA entry names!")
        stop()
    }

    # check the characters
    # there should be nothing except the AAs and modification character
    # ids <- ids %>%
    #    mutate_(TrimmedPeptide := sub(".\\.(.*)\\..", "\\1", peptide_col))
    ids <- ids %>%
        mutate(TrimmedPeptide = sub(".\\.(.*)\\..", "\\1", .[[peptide_col]]))
    present_chars <- paste0(ids$TrimmedPeptide, collapse = '') %>%
        strsplit(split='') %>%
        `[[`(1) %>%
        unique()
    other_chars <- setdiff(present_chars, c(AA_STANDARD, mod_char))
    if(length(other_chars) > 0){
        message("Detected extra chararacters in the peptide sequences!")
        # erase other chars in TrimmedPeptide
        other_chars_pttrn <- other_chars %>%
            map_chr(~paste0("\\",.x)) %>%
            paste0(collapse='') %>%
            paste0("[",.,"]")
        ids <- ids %>%
            mutate(TrimmedPeptide = str_replace_all(TrimmedPeptide, other_chars_pttrn, ""))
    }

    # check for missing peptide sequences
    if(any(is.na(ids$TrimmedPeptide))){
        message("Entries with missing peptide sequences were detected.
                The correspoding rows will be removed.")
        ids <- ids %>% filter(!is.na(TrimmedPeptide))
    }

    # check if there are peptides without mods
    mod_char_pttrn <- paste0("[\\", mod_char, "]")
    if(any(!str_detect(ids$TrimmedPeptide, mod_char_pttrn))){
        message("Peptide with no PTMs were detected. The correspoding rows will be removed.")
        ids <- ids %>%
            filter(str_detect(TrimmedPeptide, mod_char_pttrn))
    }

    # extract clean sequence
    mod_char_pttrn <- paste0("[\\", mod_char, "]")
    ids <- ids %>%
        mutate(CleanSeq = str_replace_all(TrimmedPeptide, mod_char_pttrn, ""))

    # merger of identifications and FASTA
    res <- fasta %>%
        as.data.frame() %>%
        rownames_to_column(prot_id_col) %>%
        rename(ProtSeq = x) %>%
        inner_join(ids, ., by = prot_id_col)

    # the core part
    res <- res %>%
        # locating peptide within protein
        mutate(PepLoc = map2(ProtSeq, CleanSeq, ~ as.numeric(str_locate_all(.x, .y)[[1]][,1])),
               PepLocFirst = map(PepLoc, ~ .[1])) %>%
        # adding protein length
        mutate(ProtLength = str_length(ProtSeq)) %>%
        # drop protein sequence
        select(-ProtSeq) %>%
        # locations of PTM within peptide
        mutate(ModShift = map(TrimmedPeptide, ~ as.numeric(str_locate_all(.,mod_char_pttrn)[[1]][,1])),
               ModShift = map(ModShift, ~ . - seq_along(.) - 1)) %>%
        # extract AAs
        mutate(ModAAs = map2(CleanSeq, ModShift, ~ str_sub(.x, .y+1, .y+1))) %>%
        # calculate site positions: peptide location + mod shift
        mutate(SiteLoc = map2(PepLoc, ModShift, ~ lapply(.x, `+`, .y))) %>%
        # create site notation: aa + location
        mutate(Site = map2(SiteLoc, ModAAs, ~ lapply(.x, function(..) paste0(.y, ..)))) %>%
        # collapsing site notation
        mutate(SiteCollapsed = map(Site, ~ lapply(.x, paste0, collapse =','))) %>%
        # picking the first one
        mutate(SiteCollapsedFirst = map(SiteCollapsed, 1))

    # clean-up
    res <- res %>%
        select(-c(TrimmedPeptide,CleanSeq))

    return(res)

}

utils::globalVariables(c(".", "AA_STANDARD", "TrimmedPeptide", "x",
                         "ProtSeq", "CleanSeq", "PepLoc", "ModShift",
                         "SiteLoc", "ModAAs", "Site", "SiteCollapsed"))



#' Resolving ambiguity in PTM mapping.
#'
#' Typically, the downstream use of data is based on gene symbols. One gene
#' symbol may map to multiple UniProt or RefSeq IDs corresponding to different
#' isoforms.
#' This utility simply retains the longest isoform per gene.
#'
#' @param ids \code{data.frame} object. Must contain 3 columns described below.
#' @param gene_id_col character. Name of the column with gene IDs in the `ids`
#'        object.
#' @param isoform_id_col character. Name of the column with protein isoform IDs
#'        in the `ids` object.
#' @param isoform_len_col character. Name of the column with protein isoform
#'        lengths in the `ids` object.
#'
#' @return
#' \code{data.frame}
#'
#' @importFrom dplyr mutate rename inner_join filter select row_number
#'             semi_join distinct one_of
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_detect str_replace_all str_locate_all str_length
#'             str_sub
#' @importFrom rlang !! parse_expr
#' @importFrom purrr map map2
#'
#' @export keep_longest_isoform_per_gene
#'
#' @examples
#' # Path to FASTA file
#' fasta_file_name <- system.file("extdata/FASTAs",
#'                                "rattus_norvegics_uniprot_2018_09.fasta.gz",
#'                                package = "MSnSet.utils")
#'
#' # Create AAStringSet object from FASTA file
#' library(Biostrings)
#' fst <- readAAStringSet(fasta_file_name, format="fasta",
#'                        nrec=-1L, skip=0L, use.names=TRUE)
#'
#' # Extract UniProt Accessions
#' names(fst) <- sub("^.*\\|(.*)\\|.*$", "\\1", names(fst))
#'
#' # PTM data
#' data(phospho_identifications_rat)
#'
#' # Get protein PTM locations
#' ids_with_sites <- map_PTM_sites(ids = ids, fasta = fst,
#'                                 prot_id_col = "UniProtAccFull",
#'                                 peptide_col = "Peptide", mod_char = "*")
#'
#' # Adding gene annotation. Note: this is rat data searched against UniProt.
#' # 10116 is rat taxonomy ID
#' URL <- "http://www.uniprot.org/uniprot/?query=organism:10116&columns=id,genes(PREFERRED)&format=tab"
#'
#' library(dplyr)
#'
#' # Add GeneMain column
#' ids_with_sites <- read.delim(URL, check.names = FALSE, stringsAsFactors = FALSE) %>%
#'    rename(GeneMain = `Gene names  (primary )`,
#'           UniProtAcc = `Entry`) %>%
#'    inner_join(ids_with_sites, ., by = "UniProtAcc")
#' nrow(ids_with_sites) # before
#'
#' # Resolve PTM mapping abiguity
#' ids_with_sites <- keep_longest_isoform_per_gene(
#'   ids = ids_with_sites, gene_id_col = "GeneMain",
#'   isoform_id_col = "UniProtAccFull", isoform_len_col = "ProtLength"
#' )
#' nrow(ids_with_sites) # after
#'
#' @md

keep_longest_isoform_per_gene <- function(ids, gene_id_col, isoform_id_col, isoform_len_col){

    # check if there are NA or "" values for genes
    if(any(is.na(ids[[gene_id_col]]))){
        message("NA values in gene column detected. The corresponding rows will be removed.")
        ids <- ids[!is.na(ids[[gene_id_col]]),]
    }
    if(any(is.na(ids[[gene_id_col]]))){
        message("\"\" values in gene column detected. The corresponding rows will be removed.")
        ids <- ids[ids[[gene_id_col]] != "",]
    }

    ids2 <- ids %>%
        select(one_of(gene_id_col, isoform_id_col, isoform_len_col)) %>%
        distinct()

    gene_id_col_expr <- parse_expr(gene_id_col)
    isoform_len_col_expr <- parse_expr(isoform_len_col)
    top_iso <- ids2 %>%
        group_by(!!gene_id_col_expr) %>%
        filter(row_number() == which.max(!!isoform_len_col_expr))

    # filter join
    res <- semi_join(ids, top_iso,
                     by=c(gene_id_col, isoform_id_col, isoform_len_col))
    return(res)

}


