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
                        sampleNames = "character",
                        ses = "missing",
                        grl = "missing",
                        modifications = "character",
                        genes = "character"),
  definition = function(.Object,
                        number,
                        sampleNames,
                        modifications,
                        genes){
    # Input checks
    assertive::assert_all_are_whole_numbers(numbers)
    assertive::assert_all_are_non_missing_nor_empty_character(modifications)
    assertive::assert_all_are_non_empty_character(sampleNames)
    if(length(numbers) != length(sampleNames)){
      stop("The same number of experiment numbers and sample names must given.",
           call. = FALSE)
    }
    browser()
    l <- lapply(number, function(i){getGffResult(.Object,i,modifications)})
    NULL  
  }
)
#' @rdname heatmapModifications
#' 
#' @export
setMethod(
  f = "heatmapModifications",
  signature = signature(.Object = "missing",
                        number = "missing",
                        sampleNames = "missing",
                        ses = "list",
                        grl = "missing",
                        modifications = "character",
                        genes = "character"),
  definition = function(ses,
                        modifications,
                        genes){
    # Input check
    RNAmodR::assert_all_are_SummarizedExperiment(ses)
    assertive::assert_all_are_non_missing_nor_empty_character(modifications)
    # get data as list of DataFrames
    modData <- .extract_modification_info_from_se(ses)
    posData <- .extract_position_info_from_se(ses)
    # get melted data for plotting
    modData <- .assemble_modification_info(modData,
                                           genes)
    posData <- .patch_modification_info_with_posdata(modData,
                                                     posData)
    modPlot <- .get_modification_heatmap_plot(modData)
    posPlot <- .get_modification_heatmap_plot(data)
    # to inhibit plotting by grid.arrange
    grDevices::pdf(file = NULL)
    grid <-  gridExtra::grid.arrange(posPlot,
                                     modPlot,
                                     ncol = 1, 
                                     nrow = 1,
                                     width = 420,
                                     height = 210)
    # to inhibit plotting by grid.arrange
    grDevices::dev.off()
    return(grid)
  }
)

.patch_modification_info_with_posdata <- function(modData,
                                                  posData){
  
}

#' @rdname heatmapModifications
#' 
#' @export
setMethod(
  f = "heatmapModifications",
  signature = signature(.Object = "missing",
                        number = "missing",
                        sampleNames = "missing",
                        ses = "missing",
                        grl = "GRangesList",
                        modifications = "character",
                        genes = "character"),
  definition = function(grl,
                        modifications,
                        genes){
    # Input check
    assertive::assert_all_are_non_missing_nor_empty_character(modifications)
    # get data as list of DataFrames
    modData <- .extract_modification_info_from_grl(grl)
    # get melted data for plotting
    modData <- .assemble_modification_info(modData,
                                           genes)
    modPlot <- .get_modification_heatmap_plot(modData)
    return(modPlot)
  }
)

.assemble_modification_info <- function(data,
                                        genes){
  res <- lapply(seq_along(data), function(i){
    x <- data[[i]]
    x <- x[as.character(x$Parent) %in% genes,]
    if(nrow(x) == 0) return(NULL)
    x$sample <- names(data[i])
    x[,c("sample","ID","Parent","RNAmodR_type","RNAmodR_signal","RNAmodR_signal_sd",
         "RNAmodR_z","RNAmodR_nbReplicates")]
  })
  res <- .clean_output(res)
  if(sum(vapply(res,nrow,numeric(1))) == 0) 
    stop("No modifications for given gene names: '",
         paste(genes, collapse = "', '"),
         "'.",
         call. = FALSE)
  res <- do.call(rbind, res)
  res
}

.get_modification_heatmap_plot <- function(data){
  requireNamespace("ggplot2", quietly = FALSE)
  data <- data[order(data$RNAmodR_type),]
  data$ID <- factor(data$ID, levels = unique(data$ID))
  
  plot <- ggplot2::ggplot(as.data.frame(data), 
                          ggplot2::aes_(x = ~sample, 
                                        y = ~ID)) +
    ggplot2::geom_raster(ggplot2::aes_(fill = ~RNAmodR_signal)) +
    ggplot2::scale_x_discrete(name = "experiment name",
                              expand = c(0,0)) +
    ggplot2::scale_y_discrete(name = "modification\nname",
                              expand = c(0,0)) +
    ggplot2::scale_fill_gradientn(colours = c("tan1","cyan","blue")) +
    ggplot2::theme_minimal() +
    ggplot2::theme_update(axis.text.x = ggplot2::element_text(angle = 90,
                                                              hjust = 1,
                                                              vjust = 0.5),
                          panel.background = ggplot2::element_rect(fill = NULL),
                          panel.grid = ggplot2::element_blank(),
                          panel.border = ggplot2::element_blank())
  return(plot)
}


#' @rdname heatmapModifications
#'
#' @export
setMethod(
  f = "saveHeatmapModifications",
  signature = signature(.Object = "RNAmodR",
                        number = "numeric",
                        sampleNames = "character",
                        ses = "missing",
                        grl = "missing",
                        modifications = "character",
                        genes = "character"),
  definition = function(.Object,
                        number,
                        sampleNames,
                        modifications,
                        genes){
    # Input checks
    assertive::assert_all_are_whole_numbers(numbers)
    assertive::assert_all_are_non_missing_nor_empty_character(modifications)
    assertive::assert_all_are_non_empty_character(sampleNames)
    
    NULL
  }
)
#' @rdname heatmapModifications
#'
#' @export
setMethod(
  f = "saveHeatmapModifications",
  signature = signature(.Object = "missing",
                        number = "missing",
                        sampleNames = "missing",
                        ses = "list",
                        grl = "missing",
                        modifications = "character",
                        genes = "character"),
  definition = function(ses,
                        modifications,
                        genes){
    # Input check
    RNAmodR::assert_all_are_SummarizedExperiment(ses)
    assertive::assert_all_are_non_missing_nor_empty_character(modifications)
    
    NULL
  }
)
#' @rdname heatmapModifications
#'
#' @export
setMethod(
  f = "saveHeatmapModifications",
  signature = signature(.Object = "missing",
                        number = "missing",
                        sampleNames = "missing",
                        ses = "missing",
                        grl = "GRangesList",
                        modifications = "character",
                        genes = "character"),
  definition = function(grl,
                        modifications,
                        genes){
    # Input check
    assertive::assert_all_are_non_missing_nor_empty_character(modifications)
    
    NULL
  }
)