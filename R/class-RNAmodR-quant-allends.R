#' @include class-RNAmodR-ident.R
#' @include class-RNAmodR-quant.R
NULL

#' @rdname RNAmodR-quant-class
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
#' ad <- new("RNAmodRquant_allends")
#' }
setClass("RNAmodRquant_allends",
         contains = "RNAmodRquant",
         prototype = list(dataType = "allends",
                          dataLabel = "mean(relative arrest rate)",
                          dataFormat = scales::percent))


#' @rdname quantifiyReadDataPerTranscript
#' @export
setMethod(
  f = "quantifiyReadDataPerTranscript",
  signature = signature(x = "RNAmodRquant_allends",
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
      FUN = .get_position_data_of_transcript_allends,
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
.get_position_data_of_transcript_allends <- function(data,
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
  respos <- .save_allends_counts_to_gpos(respos,
                                         data)
  respos <- .keep_non_intron_positions(respos,
                                       gff,
                                       gr$ID)
  return(respos)
}

.save_allends_counts_to_gpos <- function(gpos,
                                         data){
  tstart <- as.data.frame(table(start(data)))
  colnames(tstart) <- c("pos","value")
  tend <- as.data.frame(table(end(data)))
  colnames(tend) <- c("pos","value")
  t <- merge(tstart,
             tend,
             by = "pos",
             all = TRUE)
  t[is.na(t)] <- 0
  t$value <- rowSums(t[2:3])
  gpos[pos(gpos) %in% t$pos]$counts <- t$value
  return(gpos)
}

# # converts global to local positions and modifies data accoringly
# .convert_global_ends_to_local_ends <- function(gff,
#                                                gr,
#                                                data,
#                                                posToBeRemoved){
#   # interest in read's 5' position
#   # reset to relative positions to gene start
#   if(.is_on_minus_strand(gr)){
#     start <- BiocGenerics::end(data)
#     stops <- BiocGenerics::start(data)
#   } else {
#     start <- BiocGenerics::start(data)
#     stops <- BiocGenerics::end(data)
#   }
#   ends <- c(start,stops)
#   # offset stops based on how many stops the read has passed from
#   # transcription start
#   ends <- .move_positions(ends, 
#                           posToBeRemoved, 
#                           .get_unique_strand(gr))
#   # reset to relative stops to gene start
#   if(.is_on_minus_strand(gr)){
#     ends <- BiocGenerics::end(gr) - ends + 1
#   } else {
#     ends <- ends - BiocGenerics::start(gr) + 1
#   }
#   # get stops per position
#   ends <- table(ends)
#   # get the transcript length
#   length <- width(gr) - length(unlist(posToBeRemoved))
#   # spread table with zero values to the length of transcript
#   ends <- stats::setNames(as.double(unlist(
#     lapply(1:length, 
#            function(i){
#              if(length(ends[names(ends) == i]) == 0) return(0)
#              ends[names(ends) == i]
#            }))),1:length)
#   return(ends)
# }
