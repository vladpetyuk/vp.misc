#' Plotting Individual Protein (Feature)
#' 
#' The simplest way to plot an abundace profile.
#' 
#' @param m MSnSet or ExpressionSet object
#' @param feature Name of the feature to plot. Should be in featureNames(m).
#' @param feature_name_col Name of the column in fData(m) that 
#'                         should contain name of the feature. 
#'                         The default value is NULL and matches feature
#'                         against featureNames(m).
#' @param color_by One of the varLabels(m) that will be used to color points.
#'                 Default is no coloring.
#' @param order_by One of the varLabels(m) for ordering points.
#'                 Defaul is `color_by`. If no ordering desired use NULL.
#' @export plot_feature
#' 
#' @examples 
#' data(cptac_oca)
#' plot_feature(oca.set, "NP_001077422.1")
#' plot_feature(oca.set, "NP_001077422.1", color_by = "Batch")
#' plot_feature(oca.set, "NP_001077422.1", color_by = "Batch", order_by = "iTRAQ_ID")
#' plot_feature(oca.set, feature = "NP_001077422.1", feature_name_col = "RefSeq", color_by = "Batch")

plot_feature <- function(m, feature, feature_name_col=NULL, 
                         color_by=NULL, order_by=color_by){
    
    p_data <- pData(m) %>%
        select(-matches("sample name")) %>% # in case sample.name already exists
        rownames_to_column("sample name")
    
    if(is.null(feature_name_col))
        feature_names <- featureNames(m)
    else
        feature_names <- fData(m)[[feature_name_col]]
    
    idx <- feature == feature_names
    if(sum(idx) == 0)
        stop("No such feature found in the object.")
    if(sum(idx) > 1)
        stop("Multiple features found satisfying the criteria. 
        Handing of multiple features is not yet implemented.")
    
    x <- exprs(m[idx,]) %>%
        as.data.frame() %>%
        gather(`sample name`, abundance) %>%
        inner_join(p_data, by = "sample name") %>%
        `if`(is.null(order_by), ., arrange_at(.,order_by)) %>%
        mutate(`sample name` = ordered(`sample name`, levels=`sample name`))
    p <- x %>%
        ggplot() +
        aes(x=`sample name`, y=abundance) +
        geom_point(size=3) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(feature)
    if(!is.null(color_by))
        p <- p + aes_string(color = color_by)
    p
}
