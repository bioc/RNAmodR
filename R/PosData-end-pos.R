#' @include RNAmodR.R
#' @include PosData-class.R
NULL

#' @name EndPosData
#' 
#' @title EndPosData
#' 
#' @description
#' title
NULL


#' @rdname EndPosData
#' @export
setClass(Class = "End5PosData",
         contains = "PosData",
         prototype = list(minQuality = 5L))

#' @rdname EndPosData
#' @export
setClass(Class = "End3PosData",
         contains = "PosData",
         prototype = list(minQuality = 5L))

#' @rdname EndPosData
#' @export
setClass(Class = "EndPosData",
         contains = "PosData",
         prototype = list(minQuality = 5L))



# End5PosData ------------------------------------------------------------------

#'@importFrom reshape2 acast
.get_position_data_of_transcript_ends <- function(bamFile,
                                                   ranges,
                                                   param,
                                                   type = "5prime",
                                                   args = list()){
  parentRanges <- .get_parent_annotations(ranges)
  # add some flags
  # bamWhat(param) <- c("flag","mapq")
  data <- GenomicAlignments::readGAlignments(bamFile, param = param)
  hits <- findOverlaps(data,parentRanges)
  # split results per transcript
  data <- split(subsetByOverlaps(data, parentRanges),
                subjectHits(hits))
  # factor for found and non found transcripts
  f <- as.integer(names(data))
  f_not_found <- as.integer(
    seq_along(parentRanges)[!(seq_along(parentRanges) %in% 
                                unique(subjectHits(hits)))])
  # get data for lapply
  starts <- start(data)
  ends <- end(data)
  widths <- width(parentRanges)
  strands <- as.character(strand(parentRanges)[f])
  # aggregate 5'-pos of reads based on strand information
  if(type == "5prime"){
    data <- IntegerList(lapply(seq_along(data),
                               function(i){
                                 if(strands[i] == "+"){
                                   starts[[i]]
                                 } else {
                                   ends[[i]]
                                 }
                               }))
  } else if(type == "3prime"){
    data <- IntegerList(lapply(seq_along(data),
                               function(i){
                                 if(strands[i] == "-"){
                                   starts[[i]]
                                 } else {
                                   ends[[i]]
                                 }
                               }))
  } else if(type == "all"){
    data <- IntegerList(lapply(seq_along(data),
                               function(i){
                                 c(starts[[i]],ends[[i]])
                               }))
  } else {
    stop("Something went wrong. Invalid type '", type, "'.")
  }
  # calculate tables and add empty positions
  data <- IntegerList(mapply(
    function(d,w){
      bg <- table(seq_len(w)) - 1
      d <- table(d)
      d <- reshape2::acast(data.frame(pos = as.integer(c(names(bg),names(d))),
                                     count = as.integer(c(bg,d))),
                 pos ~ .,
                 value.var = "count",
                 fun.aggregate = sum)
      as.integer(d)
    },
    data,
    widths[f],SIMPLIFY = FALSE))
  names(data) <- f
  # get data for empty transcripts
  data_not_found <- IntegerList(mapply(
    function(w){
      d <- table(seq_len(w)) - 1
      as.integer(d)
    },
    widths[f_not_found],SIMPLIFY = FALSE))
  names(data_not_found) <- f_not_found
  # merge and order based on factor numbers
  data <- c(data,data_not_found)
  data <- data[order(as.integer(names(data))),]
  # name results based on transcript ID
  names(data) <- parentRanges$ID
  data
}

#' @rdname EndPosData
#' @export
End5PosData <- function(bamfiles,
                        fasta,
                        gff,
                        ...){
  ans <- new("End5PosData",
             bamfiles,
             fasta,
             gff,
             ...)
  args <- .get_mod_data_args(...)
  ranges <- .load_annotation(ans@gff)
  sequences <- .load_transcript_sequences(ans@fasta,
                                          ranges)
  param <- .assemble_scanBamParam(ranges,
                                  ans@minQuality,
                                  ans@chromosomes)
  message("Loading 5'-end position data from BAM files...")
  data <- lapply(ans@bamfiles,
                 FUN = .get_position_data_of_transcript_ends,
                 ranges = ranges,
                 param = param,
                 type = "5prime",
                 args = args)
  names(data) <- paste0("end5.",
                        names(ans@bamfiles),
                        ".",
                        seq_along(ans@bamfiles))
  .postprocess_read_data(ans,
                         data,
                         ranges,
                         sequences)
}

#' @rdname EndPosData
#' @export
End3PosData <- function(bamfiles,
                        fasta,
                        gff,
                        ...){
  ans <- new("End3PosData",
             bamfiles,
             fasta,
             gff,
             ...)
  args <- .get_mod_data_args(...)
  ranges <- .load_annotation(ans@gff)
  sequences <- .load_transcript_sequences(ans@fasta,
                                          ranges)
  param <- .assemble_scanBamParam(ranges,
                                  ans@minQuality,
                                  ans@chromosomes)
  message("Loading 3'-end position data from BAM files...")
  data <- lapply(ans@bamfiles,
                 FUN = .get_position_data_of_transcript_ends,
                 ranges = ranges,
                 param = param,
                 type = "3prime",
                 args = args)
  names(data) <- paste0("end3.",
                        names(ans@bamfiles),
                        ".",
                        seq_along(ans@bamfiles))
  .postprocess_read_data(ans,
                         data,
                         ranges,
                         sequences)
}

#' @rdname EndPosData
#' @export
EndPosData <- function(bamfiles,
                       fasta,
                       gff,
                       ...){
  ans <- new("EndPosData",
             bamfiles,
             fasta,
             gff,
             ...)
  args <- .get_mod_data_args(...)
  ranges <- .load_annotation(ans@gff)
  sequences <- .load_transcript_sequences(ans@fasta,
                                          ranges)
  param <- .assemble_scanBamParam(ranges,
                                  ans@minQuality,
                                  ans@chromosomes)
  browser()
  message("Loading read end position data (5' and 3') from BAM files...")
  data <- lapply(ans@bamfiles,
                 FUN = .get_position_data_of_transcript_ends,
                 ranges = ranges,
                 param = param,
                 type = "all",
                 args = args)
  names(data) <- paste0("end.",
                        names(ans@bamfiles),
                        ".",
                        seq_along(ans@bamfiles))
  .postprocess_read_data(ans,
                         data,
                         ranges,
                         sequences)
}
