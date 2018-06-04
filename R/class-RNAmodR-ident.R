#' @include RNAmodR.R
#' @include class-RNAmodR.R
NULL

#' @name RNAmodR-ident-class
#' 
#' @title RNAmodR identification class
#' 
#' @description 
#' Virtual class for modification detection.
#' 
#' @slot analysisType 
#' @slot modType
#' 
#' @import methods
#' @export
setClass("RNAmodRident",
         contains = "VIRTUAL",
         slots = c(dataType = "character",
                   modType = "character",
                   param = "list"),
         prototype = list(param = list()
                          )
)
setMethod(
  f = "initialize", 
  signature = signature(.Object = "RNAmodRident"),
  definition = function(.Object, 
                        param) {
    param <- as.list(param)
    # subset to valid params, which do not overlap with build in params
    param <- param[!(names(param) %in% names(.Object@param))]
    param <- c(.Object@param, param)
    # force numeric values
    param[vapply(param, assertive::is_numeric_string, logical(1))] <-
      as.numeric(param[vapply(param, assertive::is_numeric_string, logical(1))])
    .Object@param <- param
    return(.Object)
  }
)
#' @rdname RNAmodR-ident-class
#'
#' @param object a RNAmodRident class object 
#' 
#' @export
setMethod(
  f = "show", 
  signature = signature(object = "RNAmodRident"),
  definition = function(object) {
    cat(paste0("# Data type: ", object@dataType, "\n"))
    cat(paste0("# Modification type: ", object@modType, "\n"))
    cat(paste0("# Arguments: ", object@args, "\n"))
  }
)

#' @rdname RNAmodRident-accessors
#' @aliases getDataType getModType
#'
#' @param x a RNAmodR ident class object 
#' 
#' @return character defining the modification type for this modification
#' class
#' @export
#'
#' @examples
#' \donttest{
#' getDataType(mod)
#' getModType(mod)
#' }
setMethod(
  f = "getModType", 
  signature = signature(x = "RNAmodRident"),
  definition = function(x) {
    return(x@modType)
  }
)

#' @rdname RNAmodRident-accessors
#'
#' @return character defining the plot type for this modification class
#' 
#' @export
setMethod(
  f = "getDataType", 
  signature = signature(x = "RNAmodRident"),
  definition = function(x) {
    return(x@dataType)
  }
)

.check_modification_type <- function(modifications){
  assertive::assert_all_are_non_missing_nor_empty_character(modifications)
  modifications <- unique(modifications)
  l <- lapply(.load_mod_classes(modifications),
              getModType)
  if(unlist(l) != modifications){
    stop("Modification types do not match.")
  }
  names(l) <- modifications
  l
}

.get_analysis_type <- function(modifications){
  assertive::assert_all_are_non_missing_nor_empty_character(modifications)
  modifications <- unique(modifications)
  l <- lapply(.load_mod_classes(modifications),
              getDataType)
  if(length(l) != length(modifications)){
    stop("Modification types do not match.")
  }
  names(l) <- modifications
  l
}

.get_position_offset <- function(modifications){
  assertive::assert_all_are_non_missing_nor_empty_character(modifications)
  modifications <- unique(modifications)
  l <- lapply(.load_mod_classes(modifications),
              getPositionOffset)
  if(length(l) != length(modifications)){
    stop("Modification types do not match.")
  }
  names(l) <- modifications
  l
}

#' @rdname RNAmodRident-accessors
#' 
#' @param data the data to be subset and is assumed to be a named list of
#' GRangesList. The names have to be a valid data type of a quantifier.
#'
#' @return the subset GRangesList of data per transcript of a single type
#' 
#' @export
setMethod(
  f = "subsetData", 
  signature = signature(x = "RNAmodRident"),
  definition = function(x,
                        data) {
    data <- data[[getDataType(x)]]
    # subset to nucleotide needed
    
    if(is.null(data)){
      stop("Something went wrong.")
    }
    return(data)
  }
)

#' @rdname RNAmodRident-accessors
#' 
#' @param data the data to be subset and is assumed to be a named list of
#' GRangesList. The names have to be a valid data type of a quantifier.
#'
#' @return the subset GRangesList of scores or modifications per transcript 
#' of a single type
#' 
#' @export
setMethod(
  f = "subsetResults", 
  signature = signature(x = "RNAmodRident"),
  definition = function(x,
                        data) {
    data <- data[[getModType(x)]]
    # subset to nucleotide needed
    
    if(is.null(data)){
      stop("Something went wrong.")
    }
    return(data)
  }
)

#' @rdname scoreModifications
#' 
#' @description 
#' title
#' 
#' @param x a RNAmodRident object
#' @param data the aggregated read data per file and per transcript
#' @param args a RNAmodRargs object 
#' 
#' @return the result
setMethod(
  f = "scoreModifications",
  signature = signature(x = "RNAmodRident",
                        data = "list",
                        args = "RNAmodRargs"),
  definition = function(x,
                        data,
                        args) {
    # Input check
    nucleotide <- x@param$nucleotide
    if(!(nucleotide %in% names(IUPAC_CODE_MAP))){
      stop("Invalid nucleotide '",nucleotide,"'. Valid codes are IUPAC ",
           "conform: ",
           paste(names(IUPAC_CODE_MAP), collapse = ", "),
           ".",
           call. = FALSE)
    }
    # detect modification per transcript
    res <- mapply(
    # res <- BiocParallel::bpmapply(
      FUN = .score_mod_in_transcript,
      names(data),
      data,
      MoreArgs = list(ident = x,
                      args = args,
                      nucleotide = nucleotide),
      SIMPLIFY = FALSE)
    names(res) <- names(data)
    res <- res[vapply(res, function(z){as.logical(length(z) > 0)}, logical(1))]
    res <- GRangesList(res)
    return(res)
  }
)

# detect and merge modification positions in one transcript
.score_mod_in_transcript <- function(ID,
                                     data,
                                     ident,
                                     args,
                                     nucleotide){
  # debug
  if( getOption("RNAmodR_debug") ){
    message(ID)
  }
  # get locations which match the nucleotide preference
  nucleotide <- IUPAC_CODE_MAP[names(IUPAC_CODE_MAP) == nucleotide]
  nucleotide <- substring(nucleotide, 1:nchar(nucleotide), 1:nchar(nucleotide))
  locations <- unique(unlist(lapply(data, 
                             function(d){
    pos(d[mcols(d)$nucleotide %in% nucleotide])
  })))
  locations <- locations[order(locations)]
  endLocation <- max(locations)
  # analyze the transcript
  # Retrieve scores for positions
  scores <- scoreModificationsPerTranscript(ident,
                                            locations = locations,
                                            data = data,
                                            args = args,
                                            endLocation = endLocation)
  # name the locations based on sequence position
  mcols(scores)$ID <- ID 
  mcols(scores)$ModID <- paste0(as.character(mcols(scores)$nucleotide),
                                "m_",
                                pos(scores))
  return(scores)
}

#' @rdname scoreModificationsPerTranscript
#'
#' @param x a RNAmodRident object. 
#' @param location current location to be analyzed
#' @param data a GRangesList of GPos objects
#' @param args a RNAmodRargs object
#' @param endLocation end position of the transcript
#'
#' @return the data input list with calculated scores
NULL



#' @rdname identifyModifications
#' 
#' @description 
#' title
#' 
#' @param x a RNAmodRident object
#' @param data the aggregated read data per file and per transcript
#' @param args a RNAmodRargs object 
#' 
#' @return the result
setMethod(
  f = "identifyModifications",
  signature = signature(x = "RNAmodRident",
                        data = "GRangesList",
                        args = "RNAmodRargs"),
  definition = function(x,
                        data,
                        args) {
    browser()
    # detect modification per transcript
    res <- mapply(
      # res <- BiocParallel::bpmapply(
      FUN = .identify_mod_in_transcript,
      names(data),
      data,
      MoreArgs = list(ident = x,
                      args = args),
      SIMPLIFY = FALSE)
    names(res) <- names(object@data)
    res <- res[!vapply(res, is.null, logical(1))]
    # If not results are present return NA instead of NULL
    if(is.null(res)){
      res <- NA
    }
    return(res)
  }
)

# detect and merge modification positions in one transcript
.identify_mod_in_transcript <- function(ID,
                                        data,
                                        ident,
                                        args){
  # debug
  if( getOption("RNAmodR_debug") ){
    message(ID)
  }
  # identify modifications based on scores
  modifications <- identifyModificationsPerTranscript(ident,
                                                      data = data,
                                                      args = args)
  return(modifications)
}

#' @rdname identifyModificationsPerTranscript
#'
#' @param x a RNAmodRident object. 
#' @param data a GRangesList of GPos objects
#' @param args a RNAmodRargs object
#'
#' @return a subset of input data, which pass the thresholds set
NULL
