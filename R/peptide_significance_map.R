#' Peptide significance map
#'
#' Mapping peptides onto the protein sequence with showing statisical significance.
#'
#' @param x \code{MSnSet} object.
#'
#' @param accession \code{character}. Name of the protein
#'        (or accession in general terms) to plot peptides for.
#'
#' @param p_val_from character. Name of the column in `fData` with p-values.
#'
#' @param type character. Different ways of representing p-values.
#'        Types are: "continuous" - continuous color scheme,
#'        "discrete" - discrete color scheme,
#'        "y" - p-values shown on y-axis akin Manhattan plot.
#'
#' @param name_from character. Column in `fData` that `accession` is taken from.
#'        Typically it is "accession".
#'
#' @param ... to be documented
#'
#' @return
#' \code{ggplot} object
#'
#' @importFrom dplyr filter mutate arrange distinct_at slice
#'
#' @importFrom ggplot2 ggplot geom_rect geom_hline
#' theme theme_classic scale_x_continuous
#' scale_y_continuous scale_fill_manual scale_fill_gradientn
#' xlab ylab element_text element_blank
#' element_rect element_text scale_fill_gradient
#'
#' @importFrom scales trans_new pretty_breaks scientific_format
#'
#'
#' @importFrom stringr str_detect str_replace_all str_locate_all str_length
#'             str_sub
#' @importFrom rlang !! parse_expr
#' @importFrom purrr map map_chr map2
#' @importFrom rlang sym
#'
#' @export peptide_significance_map
#'
#' @examples
#'
#' data("peptide_count_data")
#'
#' library(msmsTests)
#'
#'
#' alt.f <- "y ~ group + 1"
#' null.f <- "y ~ 1"
#' div <- colSums(exprs(peptide_count_data)) # normalization factor
#' res <- msms.glm.qlll(peptide_count_data, alt.f, null.f, div=div)
#' res$p.val.adj <- p.adjust(res$p.value, "BH")
#'
#' # adding significance metrics
#' f_data <- fData(peptide_count_data)
#' f_data$P.Value <- res[featureNames(peptide_count_data),"p.value"]
#' f_data$adj.P.Val <- res[featureNames(peptide_count_data),"p.val.adj"]
#' fData(peptide_count_data) <- f_data
#'
#' # Path to FASTA file
#' fasta <- system.file("extdata/FASTAs",
#'                      "H_sapiens_UniProt_SPROT_2021-06-20.fasta.gz",
#'                      package = "MSnSet.utils")
#'
#' peptide_count_data <- map_peptide_position(peptide_count_data,
#'                                            fasta,
#'                                            accession_col = "protein",
#'                                            peptide_col = "Sequence")
#'
#' peptide_significance_map(fData(peptide_count_data),
#'                          "sp|O15240|VGF_HUMAN",
#'                          name_from = "protein",
#'                          "adj.P.Val",
#'                          type = "y")
#'
#'


peptide_significance_map <- function(x, accession, p_val_from,
                                     type=c("discrete","continuous","y"),
                                     name_from = c("accession"),
                                     ...){

  x <- x %>%
    filter(!!rlang::sym(name_from) == !!accession) %>%
    mutate(Length = Last_AA - First_AA + 1) %>%
    arrange(First_AA, -Length)

  # === plot type selector ===
  type <- match.arg(type, several.ok = FALSE)
  fun_switch <- list(peptide_significance_map_discrete,
                     peptide_significance_map_continuous,
                     peptide_significance_map_Y)
  names(fun_switch) <- c("discrete","continuous","y")
  FUN <- fun_switch[[type]]
  p <- FUN(x, p_val_from, ...)

  # === plot title ===
  prot_name_cols <- distinct_at(x, vars(all_of(name_from)))
  if(nrow(prot_name_cols) == 1){
    prot_name <- paste0(prot_name_cols, collapse = ", ")
  }else{
    stop("name_from  isn't unique!")
  }

  p <- p + ggtitle(prot_name)

  return(p)
}



#===== helper functions below ==================================================


manhattan_arrangement_of_the_peptides <- function(prot, p_val_from){

  prot$significance <- -log10(prot[[p_val_from]])
  sig_scaler <- max(c(ceiling(max(prot$significance)), 5))
  min_y <- 0 # placeholder
  width_y <- 0.025
  prot$ymin <- min_y
  prot$ymin <- prot$ymin + prot$significance/sig_scaler
  prot$ymin <- prot$ymin - width_y/2
  prot$ymax <- prot$ymin + width_y
  return(prot)
}

compact_stacking_of_the_peptides <- function(prot){

  # setting staggered ymin
  min_y <- 0 # 0.05
  step_y <- 0.033
  width_y <- 0.025

  prot$ymin <- 0 # min_y
  if(nrow(prot) > 1){
    for(i in 2:nrow(prot)){
      current_y <- min_y
      while(TRUE){
        # is there a conflict
        max_last_residue <- prot %>%
          dplyr::slice(1:(i-1)) %>%
          dplyr::filter(ymin == current_y) %>%
          pull(Last_AA) %>%
          max()
        if(max_last_residue + 0 >= prot[i,"First_AA",drop=TRUE]){
          current_y <- current_y + step_y
        }else{
          break()
        }
      }
      prot[i,"ymin"] <- current_y
    }
  }

  prot$ymax <- prot$ymin + width_y

  return(prot)

}



plot_peptide_on_protein_map <- function(prot,
                                        fill_by,
                                        border_color = "white",
                                        border_size = NULL,
                                        aa_step = 20){


  prot_len <- prot %>%
    distinct(ProtLen) %>%
    as.numeric()

  if(is.null(border_size)){
    border_size <- 1.5/log10(prot_len)
  }


  width_y <- prot$ymax[[1]] - prot$ymin[[1]]

  p <-
    ggplot(data = prot) +
    geom_rect(aes(xmin = 1 - 0.5,
                  xmax = prot_len + 0.5,
                  ymin = -width_y*2,
                  ymax = -width_y
    )) +
    geom_rect(aes(xmin=First_AA - 0.5,
                  xmax=Last_AA + 0.5,
                  ymin=ymin,
                  ymax=ymax,
                  fill=!!rlang::sym(fill_by)),
              color=border_color,
              size=border_size) +
    theme_classic() +
    xlab("residue") +
    ylab(NULL) +
    scale_x_continuous(breaks = c(1, seq(aa_step,prot_len,aa_step)),
                       expand = c(0.02, 0.02)) +
    theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
          axis.line.x = element_blank(),
          legend.key = element_rect(color = "black"),
          plot.title = element_text(hjust = 0.5, size=16))

  return(p)

}



plot_colored_peptide_map <- function(prot, fill_by, ...){

  # arranging peptide positions
  prot <- compact_stacking_of_the_peptides(prot) # returns with ymin and ymax

  p <- plot_peptide_on_protein_map(prot, fill_by, ...)

  width_y <- prot$ymax[[1]] - prot$ymin[[1]]
  p <- p + scale_y_continuous(
    expand =c(0, 0),
    limits = if (max(prot$ymax) < 0.45) c(-width_y*2, 0.45) else NULL) +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank())

  return(p)
}



peptide_significance_map_discrete <- function(x, p_val_from, ...){

  cuts <- c(0,0.0001, 0.001, 0.01, 0.05, 1)
  x[[p_val_from]] <- cut(x[[p_val_from]], cuts, ordered_result=T)

  fillings <- scale_fill_manual(
    breaks = levels(x[[p_val_from]]),
    values = c("#E53935", "#FDD835", "#43A047", "#1E88E5", "#BDBDBD"),
    limits = levels(x[[p_val_from]]),
    labels = c("< 1e-04", "(1e-04, 1e-03]", "(1e-03, 0.01]",
               "(0.01, 0.05]", "> 0.05"),
  )

  p <- plot_colored_peptide_map(x, p_val_from, ...)
  p <- p + fillings
  return(p)
}



peptide_significance_map_continuous <- function(x, p_val_from, ...){

  p <- plot_colored_peptide_map(x, p_val_from, ...)

  log10_rev_trans <- trans_new("log10_rev",
                               transform = function(x) -log10(x),
                               inverse = function(x) 10 ^ (-x),
                               breaks = function(x) rev(log_breaks(10)(x)))

  breaks <- min(p$data[[p_val_from]], na.rm = TRUE) %>%
    {signif(10 ^ (pretty_breaks()(0 : floor(log10(.) * 2)) / 2), 1)}

  labels <- ifelse(breaks > 1e-3,
                   format(breaks, scientific = FALSE, drop0trailing = TRUE),
                   format(breaks))

  p <- p + scale_fill_gradientn(trans=log10_rev_trans, colours=MSnSet.utils::jet2.colors(21),
                                breaks=breaks, labels=breaks, limits=c(1, min(breaks)))
  return(p)
}



peptide_significance_map_Y <- function(prot, p_val_from, fill_by = NULL, ...){


  # arranging peptide positions
  prot <- manhattan_arrangement_of_the_peptides(prot, p_val_from)

  if(is.null(fill_by))
    fill_by <- p_val_from

  p <- plot_peptide_on_protein_map(prot, fill_by = fill_by, ...)

  prot$significance <- -log10(prot[[p_val_from]])
  sig_scaler <- max(c(ceiling(max(prot$significance)), 5))
  num_labels <- sig_scaler
  if(num_labels > 10)
    num_labels <- 10
  if(num_labels < 5)
    num_labels <- 5
  # some integer number between 5 and 10, inclusive

  lblz <- pretty_breaks(num_labels)(0:sig_scaler)
  width_y <- prot$ymax[[1]] - prot$ymin[[1]]
  p <- p + scale_y_continuous(breaks = lblz/sig_scaler,
                              labels = scientific_format(digits = 1)(10^-(lblz)),
                              limits = c(-2*width_y, 1 + width_y),
                              expand = c(0, 0)) +
    geom_hline(yintercept = -log10(0.05)/sig_scaler,
               color="grey",
               linetype="dashed") +
    ylab(p_val_from)

  # handling fill
  if(fill_by == p_val_from)
    p <- p + scale_fill_gradient(low="#212121", high="#212121", guide=NULL)

  return(p)
}






