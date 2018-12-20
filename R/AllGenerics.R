#' @include RNAmodR.R
NULL

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

# Modifier accessors ------------------------------------------------------
setGeneric( 
  name = "gff",
  def = function(x) standardGeneric("gff")
)
setGeneric( 
  name = "fasta",
  def = function(x) standardGeneric("fasta")
)
setGeneric( 
  name = "sequences",
  def = function(x,
                 ...) standardGeneric("sequences")
)
setGeneric( 
  name = "bamfiles",
  def = function(x) standardGeneric("bamfiles")
)
setGeneric( 
  name = "seqdata",
  def = function(x) standardGeneric("seqdata")
) 
setGeneric( 
  name = "modifications",
  def = function(x,
                 ...) standardGeneric("modifications")
) 

# Modifier functions -----------------------------------------------------------

setGeneric( 
  name = "modify",
  def = function(x,
                 ...) standardGeneric("modify")
) 


# PosData ----------------------------------------------------------------------






# Functions --------------------------------------------------------------------
