#' @include RNAmodR.R
NULL

RNAMODR_HEATMAP_WIDTH_PER_SAMPLE <- 210
RNAMODR_HEATMAP_HEIGHT_PER_GENE <- 5

#' @rdname heatmapModifications
#' @aliases heatmapModifications saveHeatmapModifications
#'
#' @title Comparative visualization of modifications
#' 
#' @description 
#' compares different experiments as a heatmap for selected genes and/or 
#' modifications
#'
#' @param x a RNAmod object or a named list of SummarizedExperiment objects containg the 
#' experimental data to be used for comparison or a named GRangesList containg the experimental data.
#' @param numbers two or more experiment numbers
#' @param sampleNames two or more names. This length of \code{numbers} and 
#' \code{sampleNames} must match
#' @param modifications name of modification to be used for analysis. These
#' should be defined as the format gene-name_nucleotide_position.
#' @param genes gene names to be included in the analysis
#'
#' @return a heatmap plot for given experiments and modifications positions
#' @export
#'
#' @import ggplot2
#' @importFrom scales rescale
#' @importFrom gridExtra arrangeGrob
#' @importFrom grDevices pdf
#'
#' @examples
#' \donttest{
#' heatmapModifications(grl,"RDN18-1_G_1575")
#' }
setMethod(
  f = "heatmapModifications",
  signature = signature(x = "RNAmodR",
                        numbers = "numeric",
                        sampleNames = "character",
                        modifications = "character",
                        genes = "character"),
  definition = function(x,
                        numbers,
                        sampleNames,
                        modifications,
                        genes){
    browser()
    # Input checks
    assertive::assert_all_are_whole_numbers(numbers)
    assertive::assert_all_are_non_empty_character(sampleNames)
    assertive::assert_all_are_non_missing_nor_empty_character(modifications)
    if(length(numbers) != length(sampleNames)){
      stop("The same number of experiment numbers and sample names must given.",
           call. = FALSE)
    }
    # get SummarizedExperiments
    ses <- lapply(numbers, function(i){getSummarizedExperiment(x,
                                                               i,
                                                               modifications)})
    return(heatmapModifications(ses,
                                modifications = modifications,
                                genes = genes)) 
  }
)
#' @rdname heatmapModifications
#' @export
setMethod(
  f = "heatmapModifications",
  signature = signature(x = "GRangesList",
                        numbers = "missing",
                        sampleNames = "missing",
                        modifications = "character",
                        genes = "character"),
  definition = function(x,
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
#' @rdname heatmapModifications
#' @export
setMethod(
  f = "heatmapModifications",
  signature = signature(x = "list",
                        numbers = "missing",
                        sampleNames = "missing",
                        modifications = "character",
                        genes = "character"),
  definition = function(x,
                        modifications,
                        genes){
    # Input check
    lapply(ses, RNAmodR::assert_is_SummarizedExperiment)
    assertive::assert_all_are_non_missing_nor_empty_character(modifications)
    # get data as list of DataFrames
    modData <- .extract_modification_info_from_ses(ses,
                                                   genes,
                                                   modifications)
    posData <- .extract_position_info_from_ses(ses,
                                               genes)
    # get melted data for plotting
    modData <- .assemble_modification_info(modData,
                                           genes)
    modData2 <- .rescale_modification_info(modData)
    # create plots
    modPlot <- .get_modification_heatmap_plot(modData,
                                              "RNAmodR_signal",
                                              "signal strength")
    modPlotZ <- .get_modification_heatmap_plot(modData,
                                               "RNAmodR_z",
                                               "Z score")
    modPlot2 <- .get_modification_heatmap_plot(modData2,
                                               "RNAmodR_signal",
                                               "signal strength\n(normalized)")
    # create grid
    grid <- gridExtra::arrangeGrob(modPlot,
                                   modPlot2,
                                   modPlotZ,
                                   nrow = 1,
                                   ncol = 3)
    return(grid)
  }
)

# splits ID of modifications into data.frame
.convert_mod_id_to_dataframe <- function(IDs){
  df <- as.data.frame(matrix(unlist(strsplit(IDs, "_")),
                             ncol = 3, 
                             byrow = TRUE),
                      stringsAsFactors = FALSE)
  colnames(df) <- c("gene","nucleotide","pos")
  df
}
# joins ID information in data.frame to character
.convert_id_dataframe_to_id <- function(df){
  paste0(df$gene,"_",df$nucleotide,"_",df$pos)
}

# rescales the intensity to 0%-100% (normalization)
.rescale_modification_info <- function(modData){
  # get information on modifications: plottype, gene, nucleotide, pos, modtype
  plotTypes <- .get_plot_types_for_modifications(unique(modData$RNAmodR_type))
  ids <- .convert_mod_id_to_dataframe(modData$ID)
  ids$origID <- modData$ID
  ids$type <- modData$RNAmodR_type
  # get information on all positions to look at
  ids <- ids[!duplicated(ids$origID),]
  rownames(ids) <- 1:nrow(ids)
  # patch data
  res <- lapply(seq_len(nrow(ids)), function(i){
    origID <- ids[i,]$origID
    x <- modData[modData$ID == origID,]
    x$RNAmodR_signal <- scales::rescale(x$RNAmodR_signal, from = c(0,max(x$RNAmodR_signal)), to = c(0,100))
    x
  })
  modData <- do.call(rbind,res[!vapply(res, is.null, logical(1))])
  return(modData)
}

# exchanges the signal data of modifications with position data
# does not produce usable output
.patch_modification_info_with_posdata <- function(modData,
                                                  posData){
  browser()
  # get information on modifications: plottype, gene, nucleotide, pos, modtype
  plotTypes <- .get_plot_types_for_modifications(unique(modData$RNAmodR_type))
  positionOffsets <- .get_position_offset(unique(modData$RNAmodR_type))
  ids <- .convert_mod_id_to_dataframe(modData$ID)
  ids$origID <- modData$ID
  ids$type <- modData$RNAmodR_type
  # expand information for all samples, even if the modification was not called
  ids <- ids[!duplicated(ids$origID),]
  rownames(ids) <- 1:nrow(ids)
  ids <- do.call(rbind,lapply(unique(modData$sample), function(x){
    sampleIDs <- ids
    sampleIDs$sample <- x
    sampleIDs
  }))
  # get data
  l <- lapply(seq_len(nrow(ids)), function(i){
    pos <- ids[i,]$pos
    gene <- ids[i,]$gene
    sample <- ids[i,]$sample
    type <- ids[i,]$type
    signalPos <- as.numeric(pos) + positionOffsets[[type]]
    if( gene %in% rownames(posData[[sample]]) ){
      return(list(
        sig = posData[[sample]][gene,plotTypes[[type]]][[1]][signalPos,]$mean,
        sig_sd = posData[[sample]][gene,plotTypes[[type]]][[1]][signalPos,]$sd))
    }
    return(list(sig = NA,
                sig_sd = NA))
  })
  ids$RNAmodR_signal <- unlist(vapply(l,"[[",numeric(1),"sig")) * 100
  ids$RNAmodR_signal_sd <- unlist(vapply(l,"[[",numeric(1),"sig_sd")) * 100
  # patch data
  l <- list(id = .convert_id_dataframe_to_id(ids),
            sample = ids$sample,
            sig = ids$RNAmodR_signal,
            sd = ids$RNAmodR_signal_sd)
  res <- lapply(seq_along(l$id), function(i){
    id <- l$id[i]
    sample <- l$sample[i]
    if( nrow(modData[modData$ID %in% id &
                     modData$sample == sample,]) != 1 ){
      res <- modData[modData$ID %in% id,]
      res <- res[1,]
      res$sample <- sample
      res$RNAmodR_signal <- l$sig[i]
      res$RNAmodR_signal_sd <- l$sd[i]
      return(res)
    }
    return(modData[modData$ID %in% id &
                     modData$sample == sample,])
  })
  modData <- do.call(rbind,res[!vapply(res, is.null, logical(1))])
  return(modData)
}

# extracts modifications information from input data (DataFrame from mcols of
# GRanges or from rowData of SummarizedExperiment)
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

.get_modification_heatmap_plot <- function(data,
                                           column,
                                           title){
  requireNamespace("ggplot2", quietly = FALSE)
  # check that column exists
  if(!(column %in% colnames(data))){
    stop("Column name is invalid.")
  }
  # optional title
  if(missing(title)){
    title <- "signal"
  }
  title <- as.character(title)
  # sort data before plotting
  data <- data[order(data$RNAmodR_type),]
  data$ID <- factor(data$ID, levels = unique(data$ID))
  # plot data
  plot <- ggplot2::ggplot(as.data.frame(data), 
                          ggplot2::aes_(x = ~sample, 
                                        y = ~ID)) +
    ggplot2::geom_raster(ggplot2::aes_string(fill = column)) +
    ggplot2::scale_x_discrete(name = "experiment name",
                              expand = c(0,0)) +
    ggplot2::scale_y_discrete(name = "modification\nname",
                              expand = c(0,0)) +
    ggplot2::scale_fill_gradientn(name = title,
                                  limits = c(0,100),
                                  colours = c("tan1","cyan","blue"),
                                  na.value = "white") +
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
  signature = signature(x = "RNAmodR",
                        numbers = "numeric",
                        sampleNames = "character",
                        modifications = "character",
                        genes = "character",
                        folder = "missing"),
  definition = function(x,
                        numbers,
                        sampleNames,
                        modifications,
                        genes){
    # Input checks
    assertive::assert_all_are_whole_numbers(numbers)
    assertive::assert_all_are_non_missing_nor_empty_character(modifications)
    assertive::assert_all_are_non_empty_character(sampleNames)
    
    # set up some parameters
    width <- length(numbers) * RNAMODR_HEATMAP_WIDTH_PER_SAMPLE
    height <- 80 + length(genes) * RNAMODR_HEATMAP_HEIGHT_PER_GENE
    
    # create folder if it does not exist
    folder  <- paste0(getOutputFolder(x),
                      "heatmap/")
    if(!assertive::is_dir(folder)){
      dir.create(folder, 
                 recursive = TRUE)
    }
    fileName <- paste0(folder,
                       "/RNAMODR_",
                       paste(sampleNames,collapse = "_"),
                       ".pdf")
    # get plot
    grid <- heatmapModifications(x,
                                 numbers = numbers,
                                 sampleNames = sampleNames,
                                 modifications = modifications,
                                 genes = genes)
    ggplot2::ggsave(fileName,
                    plot = grid,
                    units = "mm",
                    width = width,
                    height = height)
    return(invisible(TRUE))
  }
)
#' @rdname heatmapModifications
#'
#' @export
setMethod(
  f = "saveHeatmapModifications",
  signature = signature(x = "GRangesList",
                        numbers = "missing",
                        sampleNames = "missing",
                        modifications = "character",
                        genes = "character",
                        folder = "character"),
  definition = function(x,
                        modifications,
                        genes,
                        folder){
    requireNamespace("ggplot2")
    # Input check
    assertive::assert_all_are_non_missing_nor_empty_character(modifications)
    # set up some parameters
    width <- length(x) * RNAMODR_HEATMAP_WIDTH_PER_SAMPLE
    
    height <- 80 + length(genes) * RNAMODR_HEATMAP_HEIGHT_PER_GENE
    # stop if folder does not exist
    if( !assertive::is_dir(folder)){
      stop("Given output folder does not exists. Please create it first.",
           call. = FALSE)
    }
    fileName <- paste0(folder,
                       "/RNAMODR_",
                       paste(names(x),collapse = "_"),
                       ".pdf")
    # get plot
    grid <- heatmapModifications(x,
                                 modifications = modifications,
                                 genes = genes)
    ggplot2::ggsave(fileName,
                    plot = grid,
                    units = "mm",
                    width = width,
                    height = height)
    return(invisible(TRUE))
  }
)
#' @rdname heatmapModifications
#'
#' @export
setMethod(
  f = "saveHeatmapModifications",
  signature = signature(x = "list",
                        numbers = "missing",
                        sampleNames = "missing",
                        modifications = "character",
                        genes = "character",
                        folder = "character"),
  definition = function(x,
                        modifications,
                        genes,
                        folder){
    requireNamespace("ggplot2")
    # Input check
    lapply(ses, RNAmodR::assert_is_SummarizedExperiment)
    assertive::assert_all_are_non_missing_nor_empty_character(modifications)
    # set up some parameters
    width <- length(x) * RNAMODR_HEATMAP_WIDTH_PER_SAMPLE
    
    height <- 80 + length(genes) * RNAMODR_HEATMAP_HEIGHT_PER_GENE
    # stop if folder does not exist
    if( !assertive::is_dir(folder)){
      stop("Given output folder does not exists. Please create it first.",
           call. = FALSE)
    }
    fileName <- paste0(folder,
                       "/RNAMODR_",
                       paste(names(x),collapse = "_"),
                       ".pdf")
    # get plot
    grid <- heatmapModifications(x,
                                 modifications = modifications,
                                 genes = genes)
    ggplot2::ggsave(fileName,
                    plot = grid,
                    units = "mm",
                    width = width,
                    height = height)
    return(invisible(TRUE))
  }
)