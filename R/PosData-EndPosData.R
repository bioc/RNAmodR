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
setClass(Class = "EndPosData",
         contains = "PosData",
         prototype = list(minQuality = 5L))

# EndPosData -------------------------------------------------------------------
.get_position_data_of_transcript_allends <- function(bamFile,
                                                     ranges,
                                                     param,
                                                     args = list()){
  browser()
  # skip if transcript does not have data
  if(length(bamData) == 0) return(NULL)
  # move position based on strand
  bamData <- bamData[.is_on_correct_strand(bamData,.get_unique_strand(range))]
  # discard reads out of boundaries
  bamData <- bamData[BiocGenerics::end(bamData) <= BiocGenerics::start(range),]
  bamData <- bamData[BiocGenerics::start(bamData) >= BiocGenerics::end(range),]
  # create result GRanges object
  return(respos)
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
  message("Loading read end position data (5' and 3') from BAM files...")
  data <- lapply(ans@bamfiles,
                 FUN = .get_position_data_of_transcript_allends,
                 ranges = ranges,
                 param = param,
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
