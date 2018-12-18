#' @include RNAmodR.R
#' @include PosData-class.R
NULL

#' @name ModExperiment
#' @title ModExperiment
#' @description 
#' title
NULL

#' @rdname ModExperiment
#' @export
setClass("ModExperiment",
         contains = c("VIRTUAL"),
         slots = c(mod = "character", # this have to be populated by subclass
                   dataClass = "character", # this have to be populated by subclass
                   bamfiles = "BamFileList",
                   conditions = "factor",
                   fasta = "FaFile",
                   gff = "GFFFile",
                   data = "PosData",
                   modifications = "GRanges"))

setMethod(
  f = "initialize", 
  signature = signature(.Object = "ModExperiment"),
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

#' @rdname ModExperiment
#' @export
setMethod(
  f = "show", 
  signature = signature(object = "ModExperiment"),
  definition = function(object) {
    cat("test")
  }
)

# accessors --------------------------------------------------------------------

#' @name ModExperiment
#' @export
setMethod(f = "gff", 
          signature = signature(x = "ModExperiment"),
          definition = function(x){x@gff})
  
#' @name ModExperiment
#' @export
setMethod(f = "fasta", 
          signature = signature(x = "ModExperiment"),
          definition = function(x){x@fasta})
 
#' @name ModExperiment
#' @export
setMethod(f = "sequences", 
          signature = signature(x = "ModExperiment"),
          definition = function(x){x@data@sequences})
  
#' @name ModExperiment
#' @export
setMethod(f = "ranges", 
          signature = signature(x = "ModExperiment"),
          definition = function(x){x@data@ranges})
  
#' @name ModExperiment
#' @export
setMethod(f = "bamfiles", 
          signature = signature(x = "ModExperiment"),
          definition = function(x){x@bamfiles})
  
#' @name ModExperiment
#' @export
setMethod(f = "modifications", 
          signature = signature(x = "ModExperiment"),
          definition = function(x){x@modifications})
  
