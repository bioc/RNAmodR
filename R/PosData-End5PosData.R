#' @include RNAmodR.R
#' @include PosData-class.R
NULL

#' @name End5PosData
#' 
#' @title PileupPosData
#' 
#' @description
#' title
NULL


#' @rdname End5PosData
#' @export
setClass(Class = "End5PosData",
         contains = "PosData",
         prototype = list(minQuality = 5L))

# End5PosData ------------------------------------------------------------------

.get_position_data_of_transcript_5ends <- function(bamFile,
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
  browser()
  return(respos)
}



#' @rdname End5PosData
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
                 FUN = .get_position_data_of_transcript_5ends,
                 ranges = ranges,
                 param = param,
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
