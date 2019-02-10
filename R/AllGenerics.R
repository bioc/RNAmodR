#' @include RNAmodR.R
NULL

# Class constructors -----------------------------------------------------------

setGeneric( 
  name = "SequenceData",
  signature = c("annotation","sequences"),
  def = function(dataType, bamfiles, annotation, sequences, seqinfo, ...)
    standardGeneric("SequenceData")
) 

setGeneric( 
  name = "Modifier",
  signature = c("x"),
  def = function(className, x, annotation, sequences, seqinfo, ...)
    standardGeneric("Modifier")
) 

setGeneric( 
  name = "ModifierSet",
  signature = c("x"),
  def = function(className, x, annotation, sequences, seqinfo, ...)
    standardGeneric("ModifierSet")
)

# Modifier accessors -----------------------------------------------------------

#' @rdname Modifier-functions
#' @title Modifier/ModifierSet functions
#' @export
setGeneric( 
  name = "modifierType",
  def = function(x) standardGeneric("modifierType")
)
#' @rdname Modifier-functions
#' @export
setGeneric( 
  name = "modType",
  def = function(x) standardGeneric("modType")
)
#' @rdname Modifier-functions
#' @export
setGeneric( 
  name = "mainScore",
  def = function(x) standardGeneric("mainScore")
)
#' @rdname Modifier-functions
#' @export
setGeneric( 
  name = "settings",
  def = function(x, name = NULL) standardGeneric("settings")
)
#' @rdname Modifier-functions
#' @export
setGeneric( 
  name = "settings<-",
  def = function(x, name, value) standardGeneric("settings<-")
)
#' @rdname Modifier-functions
#' @export
setGeneric( 
  name = "seqinfo",
  def = function(x) standardGeneric("seqinfo")
)
#' @rdname Modifier-functions
#' @export
setGeneric( 
  name = "sequences",
  def = function(x, ...) standardGeneric("sequences")
)
#' @rdname Modifier-functions
#' @export
setGeneric( 
  name = "bamfiles",
  def = function(x) standardGeneric("bamfiles")
)
#' @rdname Modifier-functions
#' @export
setGeneric( 
  name = "seqData",
  def = function(x) standardGeneric("seqData")
) 
#' @rdname modify
#' @name modify
#' @export
setGeneric( 
  name = "modifications",
  def = function(x, ...) standardGeneric("modifications")
) 

# Modifier functions -----------------------------------------------------------

#' @rdname RNAmodR-internals
#' @export
setGeneric( 
  name = ".constructModRanges",
  signature = c("range","data"),
  def = function(range, data, modType, scoreFun, source, type)
    standardGeneric(".constructModRanges")
)
#' @rdname modify
#' @export
setGeneric( 
  name = "modify",
  signature = c("x"),
  def = function(x, force = FALSE) standardGeneric("modify")
)
#' @rdname plotROC
#' @export
setGeneric( 
  name = "plotROC",
  signature = c("x"),
  def = function(x, coord, ...)
    standardGeneric("plotROC")
) 

# SequenceData functions -------------------------------------------------------

#' @name RNAmodR-internals
#' @export
setGeneric( 
  name = ".getData",
  signature = c("x","grl","sequences","param"),
  def = function(x, grl, sequences, param, args)
    standardGeneric(".getData")
) 

# Modifier/SequenceData functions ----------------------------------------------

#' @rdname aggregate
#' @export
setGeneric( 
  name = "aggregate",
  signature = c("x"),
  def = function(x, ...) standardGeneric("aggregate")
)
#' @rdname aggregate
#' @export
setGeneric( 
  name = "aggregateData",
  def = function(x) standardGeneric("aggregateData")
) 
#' @rdname aggregate
#' @export
setGeneric( 
  name = "hasAggregateData",
  signature = c("x"),
  def = function(x, ...) standardGeneric("hasAggregateData")
) 


# SequenceData/Modifier/ModifierSet functions ----------------------------------

#' @rdname subset
#' @export
setGeneric( 
  name = "subsetByCoord",
  signature = c("x", "coord"),
  def = function(x, coord, ...)
    standardGeneric("subsetByCoord")
)
#' @rdname visualizeData
#' @export
setGeneric(
  name = "visualizeDataByCoord",
  signature = c("x","coord"),
  def = function(x, coord, type, window.size = 15L, ...)
    standardGeneric("visualizeDataByCoord")
)
#' @rdname visualizeData
#' @export
setGeneric(
  name = "visualizeData",
  signature = c("x"),
  def = function(x, name, from = 1L, to = 30L, type, ...)
    standardGeneric("visualizeData")
)
#' @rdname visualizeData
#' @export
setGeneric(
  name = "getDataTrack",
  signature = c("x"),
  def = function(x, name, ...)
    standardGeneric("getDataTrack")
)

# ModifierSet functions --------------------------------------------------------

#' @rdname compare
#' @export
setGeneric( 
  name = "compare",
  signature = c("x"),
  def = function(x, name, from = 1L, to = 30L, ...)
    standardGeneric("compare")
) 
#' @rdname compare
#' @export
setGeneric( 
  name = "compareByCoord",
  signature = c("x","coord"),
  def = function(x, coord, ...)
    standardGeneric("compareByCoord")
) 
#' @rdname compare
#' @export
setGeneric( 
  name = "plotCompare",
  signature = c("x"),
  def = function(x, name, from = 1L, to = 30L, normalize, ...)
    standardGeneric("plotCompare")
) 
#' @rdname compare
#' @export
setGeneric( 
  name = "plotCompareByCoord",
  signature = c("x","coord"),
  def = function(x, coord, normalize, ...)
    standardGeneric("plotCompareByCoord")
) 
