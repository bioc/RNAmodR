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
  parentRanges <- .get_parent_annotations(ranges)
  bamWhat(param) <- c("flag","mapq")
  data <- GenomicAlignments::readGAlignments(bamFile, param = param)
  hits <- findOverlaps(data,parentRanges)
  data <- split(subsetByOverlaps(data, parentRanges),
                subjectHits(hits))
  f <- as.integer(names(data))
  names <- parentRanges$ID[f]
  names(data) <- parentRanges$ID[f]
  starts <- start(data)
  ends <- end(data)
  widths <- width(parentRanges[f])
  strands <- as.character(strand(parentRanges)[f])
  data <- IntegerList(lapply(seq_along(data),
                             function(i){
                               if(strands[i] == "+"){
                                 starts[[i]]
                               } else {
                                 ends[[i]]
                               }
                             }))
  data <- IntegerList(mapply(
    function(d,w){
      bg <- table(seq_len(w)) - 1
      d <- table(d)
      d <- acast(data.frame(pos = as.integer(c(names(bg),names(d))),
                            count = as.integer(c(bg,d))),
                 pos ~ .,
                 value.var = "count",
                 fun.aggregate = sum)
      as.integer(d)
    },
    data,
    widths,SIMPLIFY = FALSE))
  names(data) <- names
  data
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
