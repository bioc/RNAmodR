#' @include class-RNAmodR-mod-type.R
#' @include class-RNAmodR-data-type.R
NULL

RNAMODR_ALLENDS_COVERAGE_MIN <- 500
RNAMODR_ALLENDS_AVR_COVERAGE_MIN <- 20

#' @rdname RNAmodR-data-class
#'
#' @description 
#' \code{data_allends}: the class can be used for analyzing both the 5'-end 
#' and 3'-ends of reads to detect specific post-transcriptional modifications.
#' Methylation of the ribose is detected with this data type.
#' 
#' @export
#' @importFrom scales scientific
#'
#' @examples
#' \donttest{
#' ad <- new("data_allends")
#' }
setClass("data_allends",
         contains = "data",
         prototype = list(plotType = "allends",
                          dataLabel = "mean(relative arrest rate)",
                          dataFormat = scales::percent))


#' @rdname .getDataOfTranscript
#'
#' @param x an object for data class. 
#' @param bamData a list of GAlignments objects. 
#' @param counts total read count in bam file. 
#' @param gff GRanges annotation data. 
#'
#' @return a named list
setMethod(
  f = ".getDataOfTranscript",
  signature = signature(x = "data_allends",
                        bamData = "GAlignments",
                        counts = "numeric",
                        gff = "GRanges"),
  definition = function(x,
                        bamData,
                        counts,
                        gff) {
    transcripts <- BiocParallel::bpmapply(
      FUN = .get_position_data_of_transcript_allends,
      bamData,
      names(bamData),
      MoreArgs = list(totalCounts,
                      gff),
      SIMPLIFY = FALSE)
    names(transcripts) <- names(bamData)
    return(transcripts)
  }
)
# For each transcript get positional data
# This can be individually done for different modification types
.get_position_data_of_transcript_allends <- function(data,
                                             id,
                                             counts,
                                             gff){
  # debug
  if( getOption("RNAmodR_debug") ){
    message(id)
  }
  # skip if transcript does not have data
  if(length(data) == 0) return(NULL)
  # get ID and GRanges
  gr <- .subset_gff_for_unique_transcript(gff, 
                                          id)
  # get a list of introns and the position which are void
  posToBeRemoved <- .get_intron_positions(gff,
                                          gr$ID)
  # move position based on strand
  data <- data[.is_on_correct_strand(data,.get_unique_strand(gr))]
  # discard reads out of boundaries
  data <- data[BiocGenerics::end(data) <= BiocGenerics::end(gr),]
  data <- data[BiocGenerics::start(data) >= BiocGenerics::start(gr),]
  # get coverage per position
  coverage <- .convert_global_coverage_to_local_coverage(gff,
                                                         gr,
                                                         data,
                                                         posToBeRemoved)
  # If average coverage is to low
  if((sum(coverage)/length(coverage)) < RNAMODR_ALLENDS_AVR_COVERAGE_MIN) {
    return(NULL)
  }
  # do position conversion to translate genomic position to local transcript
  # position. take care of introns, etc
  endData <- .convert_global_ends_to_local_ends(gff,
                                                gr,
                                                data,
                                                posToBeRemoved)
  return(list(data = endData,
              coverage = coverage))
}

# converts global to local positions and modifies data accoringly
.convert_global_ends_to_local_ends <- function(gff,
                                               gr,
                                               data,
                                               posToBeRemoved){
  # interest in read's 5' position
  # reset to relative positions to gene start
  if(.is_on_minus_strand(gr)){
    start <- BiocGenerics::end(data)
    stops <- BiocGenerics::start(data)
  } else {
    start <- BiocGenerics::start(data)
    stops <- BiocGenerics::end(data)
  }
  ends <- c(start,stops)
  # offset stops based on how many stops the read has passed from
  # transcription start
  ends <- .move_positions(ends, 
                          posToBeRemoved, 
                          .get_unique_strand(gr))
  # reset to relative stops to gene start
  if(.is_on_minus_strand(gr)){
    ends <- BiocGenerics::end(gr) - ends + 1
  } else {
    ends <- ends - BiocGenerics::start(gr) + 1
  }
  # get stops per position
  ends <- table(ends)
  # get the transcript length
  length <- width(gr) - length(unlist(posToBeRemoved))
  # spread table with zero values to the length of transcript
  ends <- stats::setNames(as.double(unlist(
    lapply(1:length, 
           function(i){
             if(length(ends[names(ends) == i]) == 0) return(0)
             ends[names(ends) == i]
           }))),1:length)
  return(ends)
}
