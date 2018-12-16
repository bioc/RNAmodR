#' @include RNAmodR.R
NULL

setGeneric( 
  name = "getMod",
  def = function(x,
                 ...) standardGeneric("getMod")
) 

# ModExperiment ----------------------------------------------------------------

setGeneric( 
  name = "match.genome",
  def = function(x,
                 ...) standardGeneric("match.genome")
)
setGeneric( 
  name = "match.annotation",
  def = function(x,
                 ...) standardGeneric("match.annotation")
) 

# ModExperiment accessors ------------------------------------------------------
setGeneric( 
  name = "gff",
  def = function(x) standardGeneric("gff")
)
setGeneric( 
  name = "fasta",
  def = function(x) standardGeneric("fasta")
)
setGeneric( 
  name = "getSeq",
  def = function(x) standardGeneric("getSeq")
)
setGeneric( 
  name = "bamfiles",
  def = function(x) standardGeneric("bamfiles")
)

# PosData ----------------------------------------------------------------------






# Functions --------------------------------------------------------------------
