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
  name = "modifiertype",
  def = function(x) standardGeneric("modifiertype")
)
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
  name = "seqData",
  def = function(x) standardGeneric("seqData")
) 
setGeneric( 
  name = "aggregateData",
  def = function(x) standardGeneric("aggregateData")
) 
setGeneric( 
  name = "modifications",
  def = function(x,
                 ...) standardGeneric("modifications")
) 

# check functions --------------------------------------------------------------

setGeneric( 
  name = "hasAggregateData",
  def = function(x,
                 ...) standardGeneric("hasAggregateData")
) 

# Modifier functions -----------------------------------------------------------

setGeneric( 
  name = "modify",
  def = function(x,
                 ...) standardGeneric("modify")
)

setGeneric( 
  name = "visualizeData",
  def = function(x,
                 ...) standardGeneric("visualizeData")
) 

# SequenceData functions -------------------------------------------------------


# Modifier/SequenceData functions ----------------------------------------------

setGeneric( 
  name = "aggregate",
  def = function(x,
                 ...)
    standardGeneric("aggregate")
) 


# ModifierSet functions --------------------------------------------------------

setGeneric( 
  name = "ModifierSet",
  signature = "x",
  def = function(modifiertype,
                 x,
                 ...)
    standardGeneric("ModifierSet")
) 
