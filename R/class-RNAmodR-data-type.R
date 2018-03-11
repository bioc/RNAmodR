#' @include RNAmodR.R
#' @include class-RNAmodR.R
NULL

#' @name RNAmodR-data-class
#' @title RNAmodR data class
#' @description
#' The class is virtual and has to be extended from, eg. \code{data_5end}
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
setClass("data",
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
  signature = signature(object = "data"),
  definition = function(object) {
    NULL
  }
)

#' @rdname data-accessors
#' @aliases getPlotType getPositions getModifications
#' @title Accessor for \code{data} class objects
#' 
#' @description
#' The accessor function to \code{data} class objects can be used to access
#' the data saved in slots of the object. See examples for available functions.
#' 
#' @param object an data object 
#' @return character defining the plot type for this data class
#' @export
#' @examples
#' \donttest{
#' getPlotType(data)
#' getPositions(data)
#' getModifications(data)
#' }
setMethod(
  f = "getPlotType", 
  signature = signature(object = "data"),
  definition = function(object) {
    return(object@plotType)
  }
)

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

#' @rdname data-accessors
#' @return a list of data as a lists of list(replicate) of DataFrame(transcript)
#' @export
setMethod(
  f = "getPositions", 
  signature = signature(object = "data"),
  definition = function(object) {
    return(object@data)
  }
)

#' @rdname data-accessors
#' @return a list of data as a lists of list(replicate) of DataFrame(transcript)
#' @export
setMethod(
  f = "getModifications", 
  signature = signature(object = "data"),
  definition = function(object) {
    return(object@modifications)
  }
)

#' @rdname data-accessors
#' @return the description of the type of position data, which can be used eg. 
#' as a label for printing
#' @export
setMethod(
  f = "getDataLabel", 
  signature = signature(object = "data"),
  definition = function(object) {
    return(object@dataLabel)
  }
)

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

#' @rdname data-accessors
#' @return a function for formating the label of the position data
#' @export
setMethod(
  f = "getDataFormat", 
  signature = signature(object = "data"),
  definition = function(object) {
    return(object@dataFormat)
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


# data class handling ------------------------------------------------------

# load classes for modification data
.load_data_classes <- function(data){
  dataClasses <- vector(mode = "list", length = length(data))
  names(dataClasses) <- data
  for(i in seq_along(data)){
    className <- paste0("data_",data[[i]])
    # try to create modification detection classes
    tryCatch(
      class <- new(className),
      error = function(e) stop("Class for '",
                               data[[i]],
                               "' data does not exist (",className,").",
                               call. = FALSE)
    )
    # if( !existsMethod("convertReadsToPositions",signature(class(class),
    #                                                       "numeric",
    #                                                       "GRanges",
    #                                                       "DataFrame") ) )
    #   stop("Function convertReadsToPositions() not defined for ",class(class))
    # if( !existsMethod("parseMod",signature(class(class),
    #                                        "GRanges",
    #                                        "FaFile",
    #                                        "list") ) )
    #   stop("Function parseMod() not defined for ",class(class))
    # if( !existsMethod("mergePositionsOfReplicates",signature(class(class),
    #                                                          "GRanges",
    #                                                          "FaFile",
    #                                                          "list") ) )
    #   stop("Function mergePositionsOfReplicates() not defined for ",
    #        class(class))
    dataClasses[[i]] <- class
  }
  return(dataClasses)
}


# modification class handling --------------------------------------------------

# load classes for modification analysis
.load_mod_classes <- function(modifications){
  modClasses <- vector(mode = "list", length = length(modifications))
  for(i in seq_along(modifications)){
    className <- paste0("mod_",modifications[[i]])
    # try to create modification detection classes
    tryCatch(
      class <- new(className),
      error = function(e) stop("Class for detecting ",
                               modifications[[i]],
                               " does not exist (",className,").",
                               call. = FALSE)
    )
    # if( !existsMethod("convertReadsToPositions",signature(class(class),
    #                                                       "numeric",
    #                                                       "GRanges",
    #                                                       "DataFrame") ) )
    #   stop("Function convertReadsToPositions() not defined for ",class(class))
    # if( !existsMethod("parseMod",signature(class(class),
    #                                        "GRanges",
    #                                        "FaFile",
    #                                        "list") ) )
    #   stop("Function parseMod() not defined for ",class(class))
    # if( !existsMethod("mergePositionsOfReplicates",signature(class(class),
    #                                                          "GRanges",
    #                                                          "FaFile",
    #                                                          "list") ) )
    #   stop("Function mergePositionsOfReplicates() not defined for ",
    #        class(class))
    modClasses[[i]] <- class
  }
  names(modClasses) <- modifications
  return(modClasses)
}


#' @rdname convertReadsToData
#' 
#' @title Convert reads to data for search for modifications
#'
#' @description
#' load the input files, eg. bam files andpre-processes the information to be 
#' used in the modification detection. 
#' 
#' @param object data object 
#' @param files files which to analyze
#' @param param ScanBamParam to use for loading the bam files
#' @param gff a GRanges object for the genome
#' 
#' @return a modifified \code{data} class object containing now the read 
#' data required for modification detection.
#' 
#' @importFrom IRanges findOverlaps extractList
#' @importFrom S4Vectors split from to
setMethod(
  f = "convertReadsToData",
  signature = signature(x = "data",
                        files = "character",
                        conditions = "character",
                        gff = "GRanges",
                        param = "ScanBamParam"),
  definition = function(x,
                        files,
                        conditions,
                        gff,
                        param) {
    x@conditions <- conditions
    # detect modifications in each file
    data <- lapply(files,
                   FUN = .get_transcript_data,
                   gff,
                   param,
                   x)
    names(data) <- conditions
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
                  FUN = .get_data,
                  data = data,
                  conditions = names(data))
    names(res) <- IDs
    res <- res[!vapply(res,is.null,logical(1))]
    x@data <- res
    return(x)
  }
)

# returns the position data for analysis as a list of data per replicate for
# individual transcript
.get_data <- function(ID,
                      data,
                      conditions){
  res <- lapply(data,"[[",ID)
  res <- res[!vapply(res,is.null,logical(1))]
  return(res)
}

# detect modifications in each file
.get_transcript_data <- function(bamFile,
                                 gff,
                                 param,
                                 x){
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
  transcripts <- .getDataOfTranscript(x,
                                      bamData,
                                      totalCounts,
                                      gff)
  # remove entries for transcript for which position data is sufficient
  transcripts <- transcripts[!vapply(transcripts, is.null, logical(1))]
  if(length(transcripts) == 0) return(NULL)
  return(transcripts)
}


#' @rdname parseMod
#' 
#' @description 
#' \code{data_5end}: parses the pre-processed data in the data slot
#' to detect modifications
#' 
#' @return a modified \code{data} class object containting pre-processed
#' data and found modifications
#' 
#' @importFrom stringr str_locate_all
setMethod(
  f = "parseMod",
  signature = signature(object = "data",
                        gff = "GRanges",
                        fafile = "FaFile",
                        modClasses = "list"),
  definition = function(object,
                        gff,
                        fafile,
                        modClasses) {
    # detect modification per transcript
    # res <- mapply(
    res <- BiocParallel::bpmapply(
      FUN = .analyze_transcript_prep,
      names(object@data),
      object@data,
      MoreArgs = list(gff = gff,
                      fafile = fafile,
                      modClasses = modClasses),
      SIMPLIFY = FALSE)
    names(res) <- names(object@data)
    res <- res[!vapply(res, is.null, logical(1))]
    # If not results are present return NA instead of NULL
    if(is.null(res)){
      res <- NA
    }
    object@modifications <- res
    return(object)
  }
)

# detect and merge modification positions in one transcript
.analyze_transcript_prep <- function(ID,
                                     data,
                                     gff,
                                     fafile,
                                     modClasses){
  # debug
  if( getOption("RNAmodR_debug") ){
    message(ID)
  }
  # get sequence of transcript and subset gff for single transcript data
  gr <- .subset_gff_for_unique_transcript(gff, ID)
  # get the genomic sequences
  seq <- .get_seq_for_unique_transcript(gr,fafile,ID)
  # get transcript sequence by removing intron sequences
  seq <- .get_transcript_sequence(gff,gr$ID,seq)
  # generate a location vector
  locations <- 1:length(seq)
  names(locations) <- seq
  # analyze the transcript
  # Retrieve modifications positions
  modifications <- lapply(locations,
                          .check_for_modification,
                          modClasses,
                          data,
                          locations)
  # name the locations based on sequence position
  names(modifications) <- paste0(ID,
                                 "_",
                                 names(locations),
                                 "_",
                                 locations)
  modifications <- modifications[!vapply(modifications,is.null,logical(1))]
  if(is.null(modifications) || length(modifications) == 0) return(NULL)
  modifications <- modifications[order(as.numeric(unlist(lapply(modifications, 
                                                                "[[", 
                                                                "location"))))]
  return(modifications)
}

# check for modifications at given position
.check_for_modification <- function(location, 
                                    modClasses,
                                    data,
                                    locations){
  res <- lapply(modClasses, function(class){
    return(checkForModification(class,
                                location = location,
                                locations = locations,
                                data = data))
  })
  res <- res[!vapply(res, is.null, logical(1))]
  if(length(res) != 1) return(NULL)
  # Return data
  return(list(location = location,
              type = res[[1]]$type,
              signal = res[[1]]$signal,
              signal.sd = res[[1]]$signal.sd,
              z = res[[1]]$z,
              nbsamples = res[[1]]$nbsamples))
}


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
  signature = signature(x = "data"),
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
