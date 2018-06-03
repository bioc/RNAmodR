#' @include class-RNAmodR-ident.R
#' @include class-RNAmodR-quant.R
NULL

#' @rdname RNAmodR-quant-class
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
#' ad <- new("RNAmodRquant_5end")
#' }
setClass("RNAmodRquant_5end",
         contains = "RNAmodRquant",
         prototype = list(dataType = "5end",
                          dataLabel = "mean(relative arrest rate)",
                          dataFormat = scales::percent))

#' @rdname quantifiyReadDataPerTranscript
#' @export
setMethod(
  f = "quantifiyReadDataPerTranscript",
  signature = signature(x = "RNAmodRquant_5end",
                        bamData = "GAlignmentsList",
                        args = "RNAmodRargs",
                        counts = "numeric",
                        gff = "GRanges",
                        fafile = "FaFile"),
  definition = function(x,
                        bamData,
                        args,
                        counts,
                        gff,
                        fafile) {
    transcriptPos <- mapply(
    # transcriptPos <- BiocParallel::bpmapply(
      FUN = .get_position_data_of_transcript_5end,
      bamData,
      names(bamData),
      MoreArgs = list(args,
                      counts,
                      gff,
                      fafile),
      SIMPLIFY = FALSE)
    names(transcriptPos) <- names(bamData)
    return(transcriptPos)
  }
)
# For each transcript get positional data
# This can be individually done for different modification types
.get_position_data_of_transcript_5end <- function(data,
                                                  id,
                                                  args,
                                                  counts,
                                                  gff,
                                                  fafile){
  # debug
  if( getOption("RNAmodR_debug") ){
    message(id)
  }
  # skip if transcript does not have data
  if(length(data) == 0) return(NULL)
  # get ID and GRanges
  gr <- .subset_gff_for_unique_transcript(gff, 
                                          id)
  # get the genomic sequences
  seq <- .get_seq_for_unique_transcript(gr,fafile)
  # move position based on strand
  data <- data[.is_on_correct_strand(data,.get_unique_strand(gr))]
  # discard reads out of boundaries
  data <- data[BiocGenerics::end(data) <= BiocGenerics::end(gr),]
  data <- data[BiocGenerics::start(data) >= BiocGenerics::start(gr),]
  # create result GRanges object
  respos <- .create_per_position_grange(gr,
                                        seq)
  respos <- .save_5end_counts_to_gpos(respos,
                                       data)
  respos <- .keep_non_intron_positions(respos,
                                       gff,
                                       gr$ID)
  return(respos)
}

.save_5end_counts_to_gpos <- function(gpos,
                                      data){
  t <- table(start(data))
  t <- as.data.frame(t)
  colnames(t) <- c("pos","value")
  gpos[pos(gpos) %in% t$pos]$counts <- t$value
  return(gpos)
}

# .calc_relative_stop_data <- function(stopsData,
#                                      coverage){
#   # remove low coverage positions by setting pos data to zero
#   toLowCoverage <- as.numeric(names(coverage[coverage < RNAMODR_5END_COVERAGE_MIN]))
#   # avoid deviding by zero
#   coverage[toLowCoverage] <- 1
#   # calc relative stop data
#   res <- stopsData / coverage
#   # this quenches unsepecific results
#   # offset to not cover the end of coverage
#   res[c(1,toLowCoverage)] <- NA
#   res
# }
# 
# # converts global to local positions and modifies data accoringly
# .convert_global_stops_to_local_stops <- function(gff,
#                                                  gr,
#                                                  data,
#                                                  posToBeRemoved){
#   # interest in read's 5' position
#   # reset to relative positions to gene start
#   if(.is_on_minus_strand(gr)){
#     stops <- BiocGenerics::end(data)
#   } else {
#     stops <- BiocGenerics::start(data)
#   }
#   # offset stops based on how many stops the read has passed from
#   # transcription start
#   stops <- .move_positions(stops, 
#                            posToBeRemoved, 
#                            .get_unique_strand(gr))
#   # reset to relative stops to gene start
#   if(.is_on_minus_strand(gr)){
#     stops <- BiocGenerics::end(gr) - stops + 1
#   } else {
#     stops <- stops - BiocGenerics::start(gr) + 1
#   }
#   # get stops per position
#   stops <- table(stops)
#   # get the transcript length
#   length <- width(gr) - length(unlist(posToBeRemoved))
#   # spread table with zero values to the length of transcript
#   stops <- stats::setNames(as.double(unlist(
#     lapply(1:length, 
#            function(i){
#              if(length(stops[names(stops) == i]) == 0) return(0)
#              stops[names(stops) == i]
#            }))),1:length)
#   return(stops)
# }
