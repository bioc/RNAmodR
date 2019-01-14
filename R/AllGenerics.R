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

# ModifierSet functions --------------------------------------------------------
setGeneric( 
  name = "ModifierSet",
  signature = "x",
  def = function(modifiertype,
                 x,
                 ...)
    standardGeneric("ModifierSet")
)

# Modifier accessors ------------------------------------------------------
setGeneric( 
  name = "modifierType",
  def = function(x) standardGeneric("modifierType")
)
setGeneric( 
  name = "modType",
  def = function(x) standardGeneric("modType")
)
setGeneric( 
  name = "mainScore",
  def = function(x) standardGeneric("mainScore")
)
setGeneric( 
  name = "settings",
  def = function(x,name) standardGeneric("settings")
)
setGeneric( 
  name = "settings<-",
  def = function(x,name,value) standardGeneric("settings<-")
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
  signature = c("x"),
  def = function(x,
                 force = FALSE) standardGeneric("modify")
)

setGeneric(
  name = "visualizeDataByCoord",
  signature = c("x","coord"),
  def = function(x,
                 coord,
                 type,
                 window.size = 15L,
                 ...) standardGeneric("visualizeDataByCoord")
)
setGeneric(
  name = ".dataTracksByCoord",
  signature = c("x","data","seqdata","sequence"),
  def = function(x,
                 data,
                 seqdata,
                 sequence,
                 args) standardGeneric(".dataTracksByCoord")
)

# SequenceData functions -------------------------------------------------------


# Modifier/SequenceData functions ----------------------------------------------

setGeneric( 
  name = "aggregate",
  signature = c("x"),
  def = function(x,
                 ...)
    standardGeneric("aggregate")
) 

# Modifier/ModifierSet functions ----------------------------------------------
setGeneric( 
  name = "subsetByCoord",
  signature = c("x","coord"),
  def = function(x,
                 coord,
                 ...)
    standardGeneric("subsetByCoord")
)
setGeneric( 
  name = "plotROC",
  signature = c("x"),
  def = function(x,
                 coord,
                 redo = FALSE,
                 ...)
    standardGeneric("plotROC")
) 

# ModifierSet functions ----------------------------------------------

setGeneric( 
  name = "compareByCoord",
  signature = c("x","coord"),
  def = function(x,
                 coord,
                 ...)
    standardGeneric("compareByCoord")
) 

setGeneric( 
  name = "plotCompareByCoord",
  signature = c("x","coord"),
  def = function(x,
                 coord,
                 normalize,
                 ...)
    standardGeneric("plotCompareByCoord")
) 
