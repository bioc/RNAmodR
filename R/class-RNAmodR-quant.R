#' @include RNAmodR.R
#' @include class-RNAmodR.R
NULL

#' @name RNAmodR-quant-class
#' @title RNAmodR quant class
#' @description
#' The class is virtual and has to be extended from, eg. \code{RNAmodRquant_5end}
#' The data class is the main function for detecting modifications. It works
#' usually in a three step manner.
#' \itemize{
#'   \item{"1."}{\code{convertReadsToPositions}: This function aggregates the
#'   information of reads per bam file, eg. per replicate and stores them for
#'   subsequent analysis}
#'   \item{"2."}{\code{parseMod}: This function searches for a defined pattern. 
#'   For this \code{mod} class objects will be created on the fly based on the 
#'   names of the modification one wants to detect and which can understand the
#'   data aggregated by \code{convertReadsToPositions}.}
#'   \item{"3."}{\code{mergePositionsOfReplicates}: This function aggregates the
#'   base data further for storage. This involves usually generating a per
#'   position mean and sd or equivalent.}
#' }
#' 
#' @slot plotType 
#' @slot data 
#' @slot modifications 
#' @return a RNAmodR data class object
#' @import methods
#' @export
setClass("RNAmodRquant",
         contains = "VIRTUAL",
         slots = c(dataLabel = "character",
                   dataFormat = "function",
                   plotType = "character",
                   data = "list",
                   conditions = "character",
                   modifications = "list"),
         prototype = list(data = list(),
                          conditions = c(),
                          modifications = list())
)
#' @rdname RNAmodR-data-class
#' @param object a RNAmodR data object. 
#' @export
setMethod(
  f = "show", 
  signature = signature(object = "RNAmodRquant"),
  definition = function(object) {
    NULL
  }
)

#' @rdname RNAmodRquant-accessors
#' @aliases getPlotType getPositions getModifications
#' @title Accessor for \code{RNAmodRquant} class objects
#' 
#' @description
#' The accessor function to \code{RNAmodRquant} class objects can be used to access
#' the RNAmodRquant saved in slots of the object. See examples for available functions.
#' 
#' @param object an RNAmodRquant object 
#' @return character defining the plot type for this RNAmodRquant class
#' @export
#' @examples
#' \donttest{
#' getPlotType(quant)
#' getPositions(quant)
#' getModifications(quant)
#' }
setMethod(
  f = "getPlotType", 
  signature = signature(x = "RNAmodRquant"),
  definition = function(x) {
    return(x@plotType)
  }
)

#' @rdname RNAmodRquant-accessors
#' @return a list of data as a lists of list(replicate) of DataFrame(transcript)
#' @export
setMethod(
  f = "getPositions", 
  signature = signature(x = "RNAmodRquant"),
  definition = function(x) {
    return(x@data)
  }
)

#' @rdname RNAmodRquant-accessors
#' @return a list of data as a lists of list(replicate) of DataFrame(transcript)
#' @export
setMethod(
  f = "getModifications", 
  signature = signature(x = "RNAmodRquant"),
  definition = function(x) {
    return(x@modifications)
  }
)

#' @rdname RNAmodRquant-accessors
#' @return the description of the type of position data, which can be used eg. 
#' as a label for printing
#' @export
setMethod(
  f = "getDataLabel", 
  signature = signature(x = "RNAmodRquant"),
  definition = function(x) {
    return(x@dataLabel)
  }
)

#' @rdname RNAmodRquant-accessors
#' @return a function for formating the label of the position data
#' @export
setMethod(
  f = "getDataFormat", 
  signature = signature(x = "RNAmodRquant"),
  definition = function(x) {
    return(x@dataFormat)
  }
)

.get_data_format_for_modifications <- function(modifications){
  assertive::assert_all_are_non_missing_nor_empty_character(modifications)
  modifications <- unique(modifications)
  l <- lapply(.load_data_classes(.get_data_type(modifications)),
              getDataFormat)
  if(length(l) != length(modifications))
    stop("Incompatibble modification and data classes.")
  names(l) <- modifications
  l
}

.get_plot_types_for_modifications <- function(modifications){
  assertive::assert_all_are_non_missing_nor_empty_character(modifications)
  modifications <- unique(modifications)
  l <- lapply(.load_data_classes(.get_data_type(modifications)),
              getPlotType)
  if(length(l) != length(modifications))
    stop("Incompatibble modification and data classes.")
  names(l) <- modifications
  l
}

.get_data_label_for_modifications <- function(modifications){
  assertive::assert_all_are_non_missing_nor_empty_character(modifications)
  modifications <- unique(modifications)
  l <- lapply(.load_data_classes(.get_data_type(modifications)),
              getDataLabel)
  if(length(l) != length(modifications))
    stop("Incompatibble modification and data classes.")
  names(l) <- modifications
  l
}

#' @rdname quantifyReadData
#' 
#' @title Convert reads to data for search for modifications
#'
#' @description
#' load the input files, eg. bam files andpre-processes the information to be 
#' used in the modification detection. 
#' 
#' @param x RNAmodRquant object 
#' @param files files which to analyze
#' @param param ScanBamParam to use for loading the bam files
#' @param gff a GRanges object for the genome
#' 
#' @return the read data quantified in the flavour of each specific incarnation
#' 
#' @importFrom IRanges findOverlaps extractList
#' @importFrom S4Vectors split from to
setMethod(
  f = "quantifyReadData",
  signature = signature(x = "RNAmodRquant",
                        args = "RNAmodRargs",
                        gff = "GRanges",
                        param = "ScanBamParam"),
  definition = function(x,
                        args,
                        gff,
                        param) {
    # detect modifications in each file
    data <- lapply(getInputFiles(args),
                   FUN = .quantify_read_data_per_file,
                   x,
                   args,
                   gff,
                   param)
    names(data) <- getConditions(args)
    data <- data[!is.null(data)]
    if(length(data) == 0){
      stop("No reads detected in any bam file :\n",
           paste(files, collapse = "\n"),
           call. = FALSE)
    }
    # Process only genes found in all treated datasets
    IDs <- lapply(data,names)
    IDs <- IDs[names(IDs) == "Treated"]
    IDs <- Reduce(intersect, IDs)
    res <- lapply(IDs,
                  FUN = .subset_to_intersecting_transcripts,
                  data = data,
                  conditions = names(data))
    names(res) <- IDs
    res <- res[!vapply(res,is.null,logical(1))]
    # attach some attributes?
    
    return(res)
  }
)

# detect modifications in each file
.quantify_read_data_per_file <- function(bamFile,
                                         quant,
                                         args,
                                         gff,
                                         param){
  # load the bamfile
  bamData <- GenomicAlignments::readGAlignments(bamFile,
                                                param = param)
  # Total counts
  totalCounts <- Rsamtools::idxstatsBam(bamFile,
                                        param = param)
  totalCounts <- sum(totalCounts$mapped)
  # process result and split into chunks based on gff
  IDs <- .get_IDs_from_scanBamParam(param)
  gff_subset <- gff[.get_unique_identifiers(gff) %in% IDs,]
  hits <- IRanges::findOverlaps(bamData,
                                gff_subset)
  bamData <- IRanges::extractList(bamData, 
                                  S4Vectors::split(
                                    S4Vectors::from(hits),
                                    as.factor(S4Vectors::to(hits))))
  names(bamData) <- .get_unique_identifiers(gff_subset)[
    as.numeric(names(bamData))]
  # check for data
  if(length(bamData) == 0){
    warning("No reads detected in bam file '",
            bamFile,
            "'")
    return(NULL)
  }
  # get data for transcript. this is analysis specific. each list entry must 
  # have aggregated read end position data and coverage data: 
  # list(data = ..., coverage = ...)
  transcripts <- quantifiyReadDataPerTranscript(quant,
                                                bamData,
                                                args,
                                                totalCounts,
                                                gff)
  # remove entries for transcript for which position data is sufficient
  transcripts <- transcripts[!vapply(transcripts, is.null, logical(1))]
  if(length(transcripts) == 0) return(NULL)
  return(transcripts)
}

# returns the position data for analysis as a list of data per replicate for
# individual transcript
.subset_to_intersecting_transcripts <- function(ID,
                                                data,
                                                conditions){
  res <- lapply(data,"[[",ID)
  res <- res[!vapply(res,is.null,logical(1))]
  return(res)
}

#' @rdname quantifiyReadDataPerTranscript
#'
#' @param x a RNAmodRquant object. 
#' @param bamData a list of GAlignments objects. 
#' @param counts total read count in bam file. 
#' @param gff GRanges annotation data. 
#'
#' @return a named list
NULL



#' @rdname mergeDataOfReplicates
#' 
#' @title Merge data of replicates
#'
#' @description
#' merges data from replicates for each data class.
#' 
#' @param x a mod object 
#' 
#' @return a modified \code{data} class object now containing the aggregated
#' pre-processed data of the input. This data is used for saving.
#' 
#' @importFrom matrixStats rowSds
#' @importFrom stats complete.cases
setMethod(
  f = "mergeDataOfReplicates",
  signature = signature(x = "RNAmodRargs"),
  definition = function(x) {
    res <- mapply(
      # res <- BiocParallel::bpmapply(
      FUN = .merge_positions,
      names(x@data),
      x@data)
    names(res) <- names(x@data)
    res <- res[!is.null(res)]
    # If not results are present return NA instead of NULL
    if(is.null(res)){
      res <- NA
    }
    x@data <- res
    return(x)
  }
)

# merge positions in one transcript
.merge_positions <- function(ID,data){
  # debug
  if( getOption("RNAmodR_debug") ){
    message(ID)
  }
  # get mean of coverage
  coverage <- lapply(data,"[[","coverage")
  coverage <- split(coverage,names(coverage))
  coverage <- .merge_coverage_data(coverage)
  # get mean of position stop data
  data <- lapply(data,"[[","data")
  positions <- as.numeric(unique(unlist(lapply(data,names))))
  data <- split(data,names(data))
  data <- .merge_position_data(data)
  df <- data.frame(pos = positions,
                   mean = data$mean,
                   sd = data$sd,
                   coverage = coverage$mean,
                   stringsAsFactors = FALSE)
  return(df)
}

# iterates on every position and calculates the difference of the means and sd
.merge_position_data <- function(data){
  treated <- data$Treated
  control <- data$Control
  if(is.null(control)){
    df <- S4Vectors::DataFrame(treated)
    colnames(df) <- names(treated)
    # remove incomplete 
    df[!stats::complete.cases(as.data.frame(df)),] <- 0
    mean <- rowMeans(as.matrix(df))
    sd <- matrixStats::rowSds(as.matrix(df))
  } else {
    df <- S4Vectors::DataFrame(append(treated,control))
    colnames(df) <- c(names(treated),names(control))
    treated <- df[,colnames(df) == "Treated",drop = FALSE]
    treated[!stats::complete.cases(as.data.frame(treated)),] <- 0
    control <- df[,colnames(df) == "Control",drop = FALSE]
    control[!stats::complete.cases(as.data.frame(control)),] <- 0
    mean <- rowMeans(as.matrix(treated)) - 
      rowMeans(as.matrix(control))
    mean[mean < 0] <- 0
    sd <- matrixStats::rowSds(as.matrix(treated)) + 
      matrixStats::rowSds(as.matrix(control))
  }
  return(list(mean = mean,
              sd = sd))
}

# iterates on every position and calculates a mean coverage
.merge_coverage_data <- function(coverage){
  treated <- coverage$Treated
  df <- data.frame(treated)
  colnames(df) <- names(treated)
  df[!complete.cases(df),] <- rep(0, ncol(df))
  mean <- rowMeans(df)
  return(list(mean = mean))
}
