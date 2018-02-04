#' @include RNAmodR.R
NULL

#' @rdname heatmapModifications
#' @aliases heatmapModifications saveHeatmapModifications
#'
#' @title Comparative visualization of modifications
#' 
#' @description 
#' compares different experiments as a heatmap for selected genes and/or 
#' modifications
#'
#' @param .Object a RNAmod object.
#' @param number an experiment number
#' @param ses a named list of SummarizedExperiment objects containg the 
#' experimental data to be used for comparison
#' @param grl a named GRangesList containg the experimental data.
#' @param genes gene names to be included in the analysis
#' @param modifications name of modification to be used for analysis. These
#' should be defined as the format gene-name_nucleotide_position.
#'
#' @return a heatmap plot for given experiments and modifications positions
#' @export
#'
#' @import ggplot2
#'
#' @examples
#' \donttest{
#' heatmapModifications(grl,"RDN18-1_G_1575")
#' }
setMethod(
  f = "heatmapModifications",
  signature = signature(.Object = "RNAmodR",
                        number = "numeric",
                        genes = "character"),
  definition = function(.Object,
                        number,
                        genes){
    NULL  
  }
)
#' @rdname heatmapModifications
#' 
#' @export
setMethod(
  f = "heatmapModifications",
  signature = signature(ses = "list",
                        modifications = "character",
                        grl = "missing",
                        .Object = "missing",
                        number = "missing",
                        genes = "missing"),
  definition = function(ses,
                        modifications){
    RNAmodR::assert_all_are_SummarizedExperiment(ses)
    assertive::assert_all_are_non_missing_nor_empty_character(modifications)
    
    # get data as list of DataFrames
    data <- .extract_modification_info_from_se(ses)
    
    
  }
)
#' @rdname heatmapModifications
#' 
#' @export
setMethod(
  f = "heatmapModifications",
  signature = signature(grl = "GRangesList",
                        modifications = "character",
                        ses = "missing",
                        .Object = "missing",
                        number = "missing",
                        genes = "missing"),
  definition = function(grl,
                        modifications){
    assertive::assert_all_are_non_missing_nor_empty_character(modifications)
    
    # get data as list of DataFrames
    data <- .extract_modification_info_from_grl(grl)
    # get melted data for plotting
    data <- .assemble_modification_info(data,
                                        modifications)
    
    
    
  }
)

.assemble_modification_info <- function(data,
                                        modifications){
  res <- lapply(seq_along(data), function(i){
    x <- data[[i]]
    x <- x[rownames(x) %in% modifications,]
    if(nrow(x) == 0) return(NULL)
    x$sample <- names(data[i])
    x[,c("sample","ID","Parent","RNAmodR_type","RNAmodR_signal","RNAmodR_signal_sd",
         "RNAmodR_z","RNAmodR_nbReplicates")]
  })
  res <- .clean_output(res)
  if(sum(vapply(res,nrow,numeric(1))) == 0) 
    stop("None of the given modifications found: '",
         paste(modifications, collapse = "', '"),
         "'.",
         call. = FALSE)
  res <- do.call(rbind, res)
  res
}

.get_modification_heatmap_plot <- function(data){
  requireNamespace("ggplot2", quietly = FALSE)
  
  
  plot <- ggplot2::ggplot(as.data.frame(data), 
                          ggplot2::aes_(x = ~sample, 
                                        y = ~ID)) +
    ggplot2::geom_raster(ggplot2::aes_(fill = ~RNAmodR_signal)) +
    ggplot2::scale_x_discrete(name = "experiment name",
                              expand = c(0,0)) +
    ggplot2::scale_y_discrete(name = "modification\nname",
                              expand = c(0,0)) +
    ggplot2::scale_fill_gradientn(colours = c("tan1","cyan","lightblue","blue")) +
    ggplot2::theme_minimal() +
    ggplot2::theme_update(axis.text.x = ggplot2::element_text(angle = 90,
                                                              hjust = 1,
                                                              vjust = 0.5),
                          panel.background = ggplot2::element_rect(fill = NULL),
                          panel.grid = ggplot2::element_blank(),
                          panel.border = ggplot2::element_blank())
}


#' @rdname heatmapModifications
#'
#' @export
setMethod(
  f = "saveHeatmapModifications",
  signature = signature(.Object = "RNAmodR",
                        number = "numeric",
                        genes = "character"),
  definition = function(.Object,
                        number,
                        genes){
    NULL
  }
)
#' @rdname heatmapModifications
#'
#' @export
setMethod(
  f = "saveHeatmapModifications",
  signature = signature(ses = "list",
                        modifications = "character",
                        grl = "missing",
                        .Object = "missing",
                        number = "missing",
                        genes = "missing"),
  definition = function(ses,
                        modifications){
    NULL
  }
)
#' @rdname heatmapModifications
#'
#' @export
setMethod(
  f = "saveHeatmapModifications",
  signature = signature(grl = "GRangesList",
                        modifications = "character",
                        ses = "missing",
                        .Object = "missing",
                        number = "missing",
                        genes = "missing"),
  definition = function(grl,
                        modifications){
    NULL
  }
)