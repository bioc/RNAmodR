#' @include RNAmodR.R
#' @include SequenceData-class.R
#' @include Modifier-utils.R
NULL

#' @name Modifier
#' @title Modifier
#' @description 
#' title
NULL

#' @rdname Modifier
#' @export
setClass("Modifier",
         contains = c("VIRTUAL"),
         slots = c(mod = "character", # this have to be populated by subclass
                   dataClass = "character", # this have to be populated by subclass
                   bamfiles = "BamFileList",
                   conditions = "factor",
                   fasta = "FaFile",
                   gff = "GFFFile",
                   data = "SequenceData",
                   modifications = "GRanges"))

setMethod(
  f = "initialize", 
  signature = signature(.Object = "Modifier"),
  definition = function(.Object,
                        bamfiles,
                        fasta,
                        gff) {
    className <- class(.Object)
    # check bam files
    bamfiles <- .norm_bamfiles(bamfiles,className)
    # check genome sequences
    fasta <- .norm_fasta(fasta,className)
    # check genome annotation
    gff <- .norm_gff(gff,className)
    # check modification ident
    mod <- .norm_mod(.Object@mod,className)
    # set clots
    .Object@bamfiles <- bamfiles
    .Object@conditions <- factor(names(bamfiles))
    .Object@fasta <- fasta
    .Object@gff <- gff
    .Object@mod <- mod
    return(.Object)
  }
)

#' @rdname Modifier
#' @export
setMethod(
  f = "show", 
  signature = signature(object = "Modifier"),
  definition = function(object) {
    cat("test")
  }
)

# accessors --------------------------------------------------------------------

# converts the genomic coordinates to transcript based coordinates
.get_modifications_per_transcript <- function(x){
  ranges <- .get_parent_annotations(ranges(x))
  modifications <- modifications(x)
  modRanges <- ranges[as.character(ranges$ID) %in% as.character(modifications$Parent),]
  modRanges <- modRanges[match(as.character(modRanges$ID),
                               as.character(modifications$Parent))]
  # modify modifcation positions from genome centric to transcript centric
  start(modifications[strand(modifications) == "+"]) <- 
    start(modifications[strand(modifications) == "+"]) - 
    start(modRanges[strand(modRanges) == "+"]) + 1L
  end(modifications[strand(modifications) == "+"]) <- 
    end(modifications[strand(modifications) == "+"]) - 
    start(modRanges[strand(modRanges) == "+"]) + 1L
  end(modifications[strand(modifications) == "-"]) <- 
    end(modRanges[strand(modRanges) == "-"]) - 
    end(modifications[strand(modifications) == "-"]) + 1L
  start(modifications[strand(modifications) == "-"]) <- 
    end(modRanges[strand(modRanges) == "-"]) - 
    start(modifications[strand(modifications) == "-"]) + 1L
  names(modifications) <- as.character(modifications$Parent)
  modifications
}

#' @name Modifier
#' @export
setMethod(f = "gff", 
          signature = signature(x = "Modifier"),
          definition = function(x){x@gff})
  
#' @name Modifier
#' @export
setMethod(f = "fasta", 
          signature = signature(x = "Modifier"),
          definition = function(x){x@fasta})
 
#' @name Modifier
#' @export
setMethod(f = "sequences", 
          signature = signature(x = "Modifier"),
          definition = 
            function(x,
                     modified = FALSE,
                     with.qualities = FALSE){
              if(!assertive::is_a_bool(modified)){
                stop("'modified' has to be a single logical value.")
              }
              if(!assertive::is_a_bool(with.qualities)){
                stop("'with.qualities' has to be a single logical value.")
              }
              if(modified == FALSE){
                ans <- x@data@sequences
              }
              if(modified == TRUE){
                mod <- .get_modifications_per_transcript(x)
                mod <- split(mod,names(mod))
                ans <- ModRNAStringSet(sequences(seqdata(x)))
                modSeqList <- ans[names(ans) %in% names(mod)]
                mod <- mod[match(names(mod),names(modSeqList))]
                ans[names(ans) %in% names(mod)] <- 
                  combineIntoModstrings(modSeqList,
                                        mod)
                if(with.qualities == TRUE){
                  browser()
                }
              }
              ans
            }
)
  
#' @name Modifier
#' @export
setMethod(f = "ranges", 
          signature = signature(x = "Modifier"),
          definition = function(x){ranges(seqdata(x))})

#' @name Modifier
#' @export
setMethod(f = "bamfiles", 
          signature = signature(x = "Modifier"),
          definition = function(x){x@bamfiles})

#' @name Modifier
#' @export
setMethod(f = "seqdata", 
          signature = signature(x = "Modifier"),
          definition = function(x){x@data})
  
#' @name Modifier
#' @export
setMethod(f = "modifications", 
          signature = signature(x = "Modifier"),
          definition = 
            function(x,
                     perTranscript = FALSE){
              if(!assertive::is_a_bool(perTranscript)){
                stop("'perTranscript' has to be a single logical value.")
              }
              if(perTranscript){
                return(.get_modifications_per_transcript(x))
              }
              x@modifications
            }
)
  
# dummy functions --------------------------------------------------------------
# these need to be implemented by each subclass

#' @name Modifier
#' @export
setMethod(f = "modify", 
          signature = signature(x = "Modifier"),
          definition = 
            function(x){
              stop("This functions needs to be implemented by '",class(x),"'.",
                   call. = FALSE)
            }
)