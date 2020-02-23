#' @include RNAmodR.R
#' @import Gviz
NULL

.SequenceTrack <- Gviz:::.SequenceTrack

.get_ModDNA_bio_color <- function(){
  alphabetNames <- alphabet(ModDNAString())
  alphabet <- rep("#33FF00",length(alphabetNames))
  names(alphabet) <- alphabetNames
  dna_color <- c("A" = "#ABD9E9", "T" = "#2C7BB6", "G" = "#D7191C",
                 "C" = "#FDAE61", "N" = "#FFFFBF")
  alphabet[match(names(dna_color),names(alphabet))] <- dna_color
  alphabet
}

.get_ModRNA_bio_color <- function(){
  alphabetNames <- alphabet(ModRNAString())
  alphabet <- rep("#33FF00",length(alphabetNames))
  names(alphabet) <- alphabetNames
  rna_color <- c("A" = "#ABD9E9", "U" = "#2C7BB6", "G" = "#D7191C",
                 "C" = "#FDAE61", "N" = "#FFFFBF")
  alphabet[match(names(rna_color),names(alphabet))] <- rna_color
  alphabet
}


################################################################################
## Gviz + ModRNAString #########################################################
################################################################################


#' @name SequenceModRNAStringSetTrack-class
#' @aliases ModRNASequenceTrack SequenceModRNAStringSetTrack
#' 
#' @title ModRNASequenceTrack
#' 
#' @description 
#' A \code{Gviz} compatible 
#' \code{\link[Gviz:SequenceTrack-class]{SequenceTrack}} for showing modified 
#' RNA sequences.
#' 
#' @slot sequence A \code{ModRNAStringSet} object
#' 
#' @export
setClass("SequenceModRNAStringSetTrack",
         representation = representation(sequence = "ModRNAStringSet"),
         contains = "SequenceTrack",
         prototype = prototype(sequence = ModRNAStringSet(),
                               dp = DisplayPars(add53 = FALSE,
                                                background.title = "transparent",
                                                col = "darkgray",
                                                complement = FALSE,
                                                fontcolor = .get_ModRNA_bio_color(),
                                                fontface = 2,
                                                fontsize = 10,
                                                lwd = 2,
                                                min.width = 2,
                                                noLetters = FALSE,
                                                showTitle = FALSE,
                                                size = NULL)))

setMethod("initialize", "SequenceModRNAStringSetTrack",
          function(.Object, sequence, ...) {
              if(missing(sequence) || is.null(sequence))
                  sequence <- ModRNAStringSet()
              .Object@sequence <- sequence
              .Object <- callNextMethod(.Object, ...)
              return(.Object)
          }
)


#' @rdname SequenceModRNAStringSetTrack-class
#' 
#' @param sequence A \code{character} vector or \code{ModRNAString} object of 
#' length one. The sequence to display. See 
#' \code{\link[Gviz:SequenceTrack-class]{SequenceTrack}}.
#' @param chromosome,genome,name,... See 
#' \code{\link[Gviz:SequenceTrack-class]{SequenceTrack}}.
#' 
#' @return a \code{SequenceModRNAStringSetTrack} object
#'   
#' @export
#' 
#' @examples
#' seq <- ModRNAStringSet(c(chr1 = paste0(alphabet(ModRNAString()),
#'                                        collapse = "")))
#' st <- ModRNASequenceTrack(seq)
#' Gviz::plotTracks(st, chromosome = "chr1",from = 1L, to = 20L)
ModRNASequenceTrack <- function(sequence, chromosome, genome,
                                name = "SequenceTrack", ...){
  .SequenceTrack("SequenceModRNAStringSetTrack",
                 seqtype = "ModRNA",
                 sequence = sequence,
                 chromosome = chromosome,
                 genome = genome,
                 name = name,
                 ...)
}

#' @rdname SequenceModRNAStringSetTrack-class
#' @export
setMethod("seqnames", "SequenceModRNAStringSetTrack",
          function(x) as.character(names(x@sequence)))
#' @rdname SequenceModRNAStringSetTrack-class
#' @export
setMethod("seqlevels", "SequenceModRNAStringSetTrack",
          function(x) seqnames(x)[width(x@sequence)>0])
setMethod("show", "SequenceModRNAStringSetTrack",
          function(object) cat(.sequenceTrackInfo(object)))


################################################################################
## Gviz + ModDNAString #########################################################
################################################################################

#' @name SequenceModDNAStringSetTrack-class
#' @aliases ModDNASequenceTrack SequenceModDNAStringSetTrack
#' 
#' @title ModDNASequenceTrack
#' 
#' @description 
#' A \code{Gviz} compatible 
#' \code{\link[Gviz:SequenceTrack-class]{SequenceTrack}} for showing modified 
#' DNA sequences.
#' 
#' @slot sequence A \code{ModDNAStringSet} object
#' 
#' @export
setClass("SequenceModDNAStringSetTrack",
         representation = representation(sequence = "ModDNAStringSet"),
         contains = "SequenceTrack",
         prototype = prototype(sequence = ModDNAStringSet(),
                               dp = DisplayPars(add53 = FALSE,
                                                background.title = "transparent",
                                                col = "darkgray",
                                                complement = FALSE,
                                                fontcolor = .get_ModDNA_bio_color(),
                                                fontface = 2,
                                                fontsize = 10,
                                                lwd = 2,
                                                min.width = 2,
                                                noLetters = FALSE,
                                                showTitle = FALSE,
                                                size = NULL)))

setMethod("initialize", "SequenceModDNAStringSetTrack",
          function(.Object, sequence, ...) {
              if(missing(sequence) || is.null(sequence))
                  sequence <- ModDNAStringSet()
              .Object@sequence <- sequence
              .Object <- callNextMethod(.Object, ...)
              return(.Object)
          }
)


#' @rdname SequenceModDNAStringSetTrack-class
#' 
#' @param sequence A \code{character} vector or \code{ModDNAString} object of 
#' length one. The sequence to display. See 
#' \code{\link[Gviz:SequenceTrack-class]{SequenceTrack}}.
#' @param chromosome,genome,name,... See 
#' \code{\link[Gviz:SequenceTrack-class]{SequenceTrack}}.
#' 
#' @return a \code{SequenceModDNAStringSetTrack} object
#'   
#' @export
#' 
#' @examples
#' seq <- ModDNAStringSet(c(chr1 = paste0(alphabet(ModDNAString()),
#'                                        collapse = "")))
#' st <- ModDNASequenceTrack(seq)
#' Gviz::plotTracks(st, chromosome = "chr1",from = 1L, to = 20L)
ModDNASequenceTrack <- function(sequence, chromosome, genome,
                                name = "SequenceTrack", ...){
  .SequenceTrack("SequenceModDNAStringSetTrack",
                 seqtype = "ModDNA",
                 sequence = sequence,
                 chromosome = chromosome,
                 genome = genome,
                 name = name,
                 ...)
}

#' @rdname SequenceModDNAStringSetTrack-class
#' @export
setMethod("seqnames", "SequenceModDNAStringSetTrack",
          function(x) as.character(names(x@sequence)))
#' @rdname SequenceModDNAStringSetTrack-class
#' @export
setMethod("seqlevels", "SequenceModDNAStringSetTrack",
          function(x) seqnames(x)[width(x@sequence)>0])
setMethod("show", "SequenceModDNAStringSetTrack",
          function(object) cat(.sequenceTrackInfo(object)))
