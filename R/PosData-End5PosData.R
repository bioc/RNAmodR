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
  .Object <- new("End5PosData",
                 bamfiles,
                 fasta,
                 gff,
                 ...)
  args <- .get_mod_data_args(...)
  ranges <- .load_annotation(.Object@gff)
  sequences <- .load_transcript_sequences(.Object@fasta,
                                          ranges)
  param <- .assemble_scanBamParam(ranges,
                                  .Object@minQuality,
                                  .Object@chromosomes)
  message("Loading data from BAM files...")
  data <- lapply(.Object@bamfiles,
                 FUN = .get_position_data_of_transcript_5ends,
                 ranges = ranges,
                 param = param,
                 args = args)
  names(data) <- paste0("end5.",
                        names(.Object@bamfiles),
                        ".",
                        seq_along(.Object@bamfiles))
  .Object <- .postprocess_read_data(.Object,
                                    data,
                                    ranges,
                                    sequences)
  return(.Object)
}
