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
                   args = "numeric")
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
    NULL
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
    
    return(data)
  }
)

#' @rdname identifyModifications
#' 
#' @description 
#' title
#' 
#' @param x a RNAmodRident object
#'
#' @param args a RNAmodRargs object 
#' @param data the qaggregated read data per file and per transcript
#' @param gff the gene annotation
#' @param fafile the genomic sequence file 
#' 
#' @return the result
#' 
#' @importFrom stringr str_locate_all
setMethod(
  f = "identifyModifications",
  signature = signature(x = "RNAmodRident",
                        args = "RNAmodRargs",
                        data = "GRangesList",
                        gff = "GRanges",
                        fafile = "FaFile"),
  definition = function(x,
                        args,
                        data,
                        gff,
                        fafile) {
    
    # get sequence of transcript and subset gff for single transcript data
    grl <- .subset_gff_for_unique_transcript(gff, ID)
    # get the genomic sequences
    seqs <- .get_seq_for_unique_transcript(gr,fafile,ID)
    # get transcript sequence by removing intron sequences
    seqs <- .get_transcript_sequence(gff,gr$ID,seq)
    
    # detect modification per transcript
    # res <- mapply(
    res <- BiocParallel::bpmapply(
      FUN = .identify_mod_in_transcript,
      names(object@data),
      object@data,
      grl,
      seqs,
      MoreArgs = list(args = args),
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
                                        gr,
                                        seq,
                                        args){
  # debug
  if( getOption("RNAmodR_debug") ){
    message(ID)
  }
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


#' @rdname identifyModificationPerTranscript
#'
#' @param x a RNAmodRident object. 
#' @param bamData a list of GAlignments objects. 
#' @param counts total read count in bam file. 
#' @param gff GRanges annotation data. 
#'
#' @return a named list
NULL
