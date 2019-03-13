#' @include RNAmodR.R
#' @import Gviz
#' @importFrom biovizBase getBioColor
#' @importFrom Biostrings subseq
NULL

################################################################################
# Since the DNA alphabet is hardcoded into the SequenceTrack class of Gviz
# it needs to be reimplmeneted a bit
################################################################################

setClass("ModifiedSequenceTrack",
         contains = "GdObject",
         representation = representation(seqType = "character"),
         prototype = prototype(seqType = "XString")
)

setMethod("initialize", 
          signature = signature("ModifiedSequenceTrack"),
          definition = function(.Object,
                                sequence,
                                chromosome = NA,
                                genome = NA,
                                ...) {
            # get sequence sanitized
            if(is.null(sequence)){
              sequence <- .Object@sequence
            }
            .Object@sequence <- sequence
            ## the diplay parameter defaults
            Gviz:::.makeParMapping()
            .Object <- Gviz:::.updatePars(.Object, "SequenceTrack")
            if(!missing(chromosome) && 
               !is.na(chromosome) && 
               !is.null(chromosome)){
              .Object@chromosome <- .chrName(names(sequence)[1])[1]
            }
            if(missing(genome) || is.na(genome) || is.null(genome)) {
              genome <- as.character(NA)
            }
            .Object@genome <- genome
            .Object <- callNextMethod(.Object, ...)
            return(.Object)
})

#' @rdname RNAmodR-internals
setMethod("show", "ModifiedSequenceTrack",
          function(object){
            cat(.sequenceTrackInfo(object))
          }
)


# constructor ------------------------------------------------------------------

.stringSet_to_ModifiedSequenceTrack <- function(modSeqTrackType,
                                                seqType,
                                                sequence = NA,
                                                chromosome = NA,
                                                genome = NA,
                                                fontcolor,
                                                ...){
  if(is.null(sequence)){
    return(new(modSeqTrackType,
               chromosome = chromosome,
               genome = genome,
               name = modSeqTrackType,
               ...))
  }
  if(!is(sequence, seqType)){
    stop("Argument sequence must be of class '",seqType,"'",
         call. = FALSE)
  }
  if(is.null(names(sequence))){
    stop("The sequences in the ",seqType," must be named",
         call. = FALSE)
  }
  if(any(duplicated(names(sequence)))){
    stop("The sequence names in the ",seqType," must be unique",
         call. = FALSE)
  }
  if(missing(chromosome) || is.na(chromosome) || is.null(chromosome)){
    chromosome <- names(sequence)[1]
  }
  obj <- new(modSeqTrackType,
             sequence = sequence,
             chromosome = chromosome,
             genome = genome,
             name = modSeqTrackType,
             ...)
  displayPars(obj)$fontcolor <- fontcolor
  return(obj)
}


################################################################################
## Gviz + RNAString ############################################################
################################################################################

#' @name RNASequenceTrack
#' @aliases RNASequenceTrack-class
#' 
#' @title RNASequenceTrack
#' 
#' @description
#' A \code{Gviz} compatible 
#' \code{\link[Gviz:SequenceTrack-class]{SequenceTrack}} for showing RNA 
#' sequences.
#' 
#' @export
setClass("RNASequenceTrack",
         contains = "ModifiedSequenceTrack",
         representation = representation(sequence = "RNAStringSet",
                                         chromosome = "character",
                                         genome = "character"),
         prototype = prototype(seqType = "RNAString",
                               sequence = RNAStringSet(),
                               name = "Sequence",
                               chromosome = "chr0",
                               genome = "all")
)

#' @rdname RNASequenceTrack
#' 
#' @param sequence A \code{character} vector or \code{RNAString} object of 
#' length one. The sequence to display. See 
#' \code{\link[Gviz:SequenceTrack-class]{SequenceTrack}}.
#' @param chromosome,genome,... See 
#' \code{\link[Gviz:SequenceTrack-class]{SequenceTrack}}.
#' 
#' @return a \code{RNASequenceTrack} object
#'   
#' @export
#'
#' @examples
#' seq <- RNAStringSet(c(chr1 = "AGCUAGCUAGCUAGCUAGCUAGCU"))
#' st <- RNASequenceTrack(seq)
#' Gviz::plotTracks(st, chromosome = "chr1",from = 1L, to = 20L)
RNASequenceTrack <- function(sequence, chromosome, genome, ...){
  .stringSet_to_ModifiedSequenceTrack("RNASequenceTrack",
                                      "RNAStringSet",
                                      sequence,
                                      chromosome,
                                      genome,
                                      .get_RNA_bio_color(),
                                      ...)
}

################################################################################
## Gviz + ModRNAString #########################################################
################################################################################

#' @name ModRNASequenceTrack
#' @aliases ModRNASequenceTrack-class
#' 
#' @title ModRNASequenceTrack
#' 
#' @description 
#' A \code{Gviz} compatible 
#' \code{\link[Gviz:SequenceTrack-class]{SequenceTrack}} for showing modified 
#' RNA sequences.
#' 
#' @export
setClass("ModRNASequenceTrack",
         contains = "ModifiedSequenceTrack",
         representation = representation(sequence = "ModRNAStringSet",
                                         chromosome = "character",
                                         genome = "character"),
         prototype = prototype(seqType = "ModRNAString",
                               sequence = ModRNAStringSet(),
                               name = "Sequence",
                               chromosome = "chr0",
                               genome = "all")
)

#' @rdname ModRNASequenceTrack
#' 
#' @param sequence A \code{character} vector or \code{ModRNAString} object of 
#' length one. The sequence to display. See 
#' \code{\link[Gviz:SequenceTrack-class]{SequenceTrack}}.
#' @param chromosome,genome,... See 
#' \code{\link[Gviz:SequenceTrack-class]{SequenceTrack}}.
#' 
#' @return a \code{ModRNASequenceTrack} object
#'   
#' @export
#' 
#' @examples
#' seq <- ModRNAStringSet(c(chr1 = paste0(alphabet(ModRNAString()),collapse = "")))
#' st <- ModRNASequenceTrack(seq)
#' Gviz::plotTracks(st, chromosome = "chr1",from = 1L, to = 20L)
ModRNASequenceTrack <- function(sequence, chromosome, genome, ...){
  .stringSet_to_ModifiedSequenceTrack("ModRNASequenceTrack",
                                      "ModRNAStringSet",
                                      sequence,
                                      chromosome,
                                      genome,
                                      .get_ModRNA_bio_color(),
                                      ...)
}
