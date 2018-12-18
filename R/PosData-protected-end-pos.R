#' @include RNAmodR.R
#' @include PosData-class.R
NULL

#' @name ProtectedEndPosData
#' 
#' @title ProtectedEndPosData
#' 
#' @description
#' title
NULL

#' @rdname ProtectedEndPosData
#' @export
setClass(Class = "ProtectedEndPosData",
         contains = "PosData",
         prototype = list(minQuality = 5L))


# ProtectedEndPosData ------------------------------------------------------------------

#' @rdname ProtectedEndPosData
#' @export
ProtectedEndPosData <- function(bamfiles,
                                fasta,
                                gff,
                                ...){
  ans <- new("ProtectedEndPosData",
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
  message("Loading protected end data from BAM files...")
  enddata <- lapply(ans@bamfiles,
                    FUN = .get_position_data_of_transcript_ends,
                    ranges = ranges,
                    param = param,
                    type = "all",
                    args = args)
  coverage <- lapply(ans@bamfiles,
                     FUN = .get_position_data_of_transcript_coverage,
                     ranges = ranges,
                     param = param,
                     args = args)
  data <- lapply(seq_along(ans@bamfiles),
                 function(i){
                   coverage[[i]] - enddata[[i]]
                 })
  names(data) <- paste0("protectedend.",
                        names(ans@bamfiles),
                        ".",
                        seq_along(ans@bamfiles))
  .postprocess_read_data(ans,
                         data,
                         ranges,
                         sequences)
}
