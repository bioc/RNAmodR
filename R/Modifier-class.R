#' @include RNAmodR.R
#' @include PosData-class.R
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
                   data = "PosData",
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
          definition = function(x){x@data@sequences})
  
#' @name Modifier
#' @export
setMethod(f = "ranges", 
          signature = signature(x = "Modifier"),
          definition = function(x){x@data@ranges})
  
#' @name Modifier
#' @export
setMethod(f = "bamfiles", 
          signature = signature(x = "Modifier"),
          definition = function(x){x@bamfiles})
  
#' @name Modifier
#' @export
setMethod(f = "modifications", 
          signature = signature(x = "Modifier"),
          definition = function(x){x@modifications})
  
