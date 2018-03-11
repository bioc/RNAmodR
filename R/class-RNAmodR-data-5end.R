#' @include class-RNAmodR-mod-type.R
#' @include class-RNAmodR-data-type.R
NULL

RNAMODR_5END_COVERAGE_MIN <- 50
RNAMODR_5END_AVR_COVERAGE_MIN <- 10

#' @rdname RNAmodR-data-class
#'
#' @description 
#' \code{data_5end}: the class can be used for analyzing the 5'-end of
#' reads to detect specific post-transcriptional modifications. This includes
#' m7G, m3C and Dihydrouridine modifications.
#' 
#' @export
#' @importFrom scales scientific
#'
#' @examples
#' \donttest{
#' ad <- new("data_5end")
#' }
setClass("data_5end",
         contains = "data",
         prototype = list(plotType = "5end",
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
  signature = signature(x = "data_5end",
                        bamData = "GAlignments",
                        counts = "numeric",
                        gff = "GRanges"),
  definition = function(x,
                        bamData,
                        counts,
                        gff) {
    # transcripts <- mapply(
    transcripts <- BiocParallel::bpmapply(
      FUN = .get_position_data_of_transcript_5end,
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
.get_position_data_of_transcript_5end <- function(data,
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
  if((sum(coverage)/length(coverage)) < RNAMODR_5END_AVR_COVERAGE_MIN) {
    return(NULL)
  }
  # do position conversion to translate genomic position to local transcript
  # position. take care of introns, etc
  stopsData <- .convert_global_stops_to_local_stops(gff,
                                                    gr,
                                                    data,
                                                    posToBeRemoved)
  # calculate relative amount of stops per coverage as percent
  stopsData <- .calc_relative_stop_data(stopsData,
                                        coverage)
  return(list(data = stopsData,
              coverage = coverage))
}

.calc_relative_stop_data <- function(stopsData,
                                     coverage){
  # remove low coverage positions by setting pos data to zero
  toLowCoverage <- as.numeric(names(coverage[coverage < RNAMODR_5END_COVERAGE_MIN]))
  # avoid deviding by zero
  coverage[toLowCoverage] <- 1
  # calc relative stop data
  res <- stopsData / coverage
  # this quenches unsepecific results
  # offset to not cover the end of coverage
  res[c(1,toLowCoverage)] <- NA
  res
}

# converts global to local positions and modifies data accoringly
.convert_global_stops_to_local_stops <- function(gff,
                                                 gr,
                                                 data,
                                                 posToBeRemoved){
  # interest in read's 5' position
  # reset to relative positions to gene start
  if(.is_on_minus_strand(gr)){
    stops <- BiocGenerics::end(data)
  } else {
    stops <- BiocGenerics::start(data)
  }
  # offset stops based on how many stops the read has passed from
  # transcription start
  stops <- .move_positions(stops, 
                           posToBeRemoved, 
                           .get_unique_strand(gr))
  # reset to relative stops to gene start
  if(.is_on_minus_strand(gr)){
    stops <- BiocGenerics::end(gr) - stops + 1
  } else {
    stops <- stops - BiocGenerics::start(gr) + 1
  }
  # get stops per position
  stops <- table(stops)
  # get the transcript length
  length <- width(gr) - length(unlist(posToBeRemoved))
  # spread table with zero values to the length of transcript
  stops <- stats::setNames(as.double(unlist(
    lapply(1:length, 
           function(i){
             if(length(stops[names(stops) == i]) == 0) return(0)
             stops[names(stops) == i]
           }))),1:length)
  return(stops)
}
