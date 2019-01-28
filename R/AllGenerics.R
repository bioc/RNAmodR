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

# Modifier accessors -----------------------------------------------------------

#' @rdname Modifier
#' @export
setGeneric( 
  name = "modifierType",
  def = function(x) standardGeneric("modifierType")
)
#' @rdname Modifier
#' @export
setGeneric( 
  name = "modType",
  def = function(x) standardGeneric("modType")
)
#' @rdname Modifier
#' @export
setGeneric( 
  name = "mainScore",
  def = function(x) standardGeneric("mainScore")
)
#' @rdname Modifier
#' @export
setGeneric( 
  name = "settings",
  def = function(x,name = NULL) standardGeneric("settings")
)
#' @rdname Modifier
#' @export
setGeneric( 
  name = "settings<-",
  def = function(x,name,value) standardGeneric("settings<-")
)
#' @rdname Modifier
#' @export
setGeneric( 
  name = "seqinfo",
  def = function(x) standardGeneric("seqinfo")
)
#' @rdname Modifier
#' @export
setGeneric( 
  name = "sequences",
  def = function(x,
                 ...) standardGeneric("sequences")
)
#' @rdname Modifier
#' @export
setGeneric( 
  name = "bamfiles",
  def = function(x) standardGeneric("bamfiles")
)
#' @rdname Modifier
#' @export
setGeneric( 
  name = "seqData",
  def = function(x) standardGeneric("seqData")
) 
#' @rdname aggregate
#' @export
setGeneric( 
  name = "aggregateData",
  def = function(x) standardGeneric("aggregateData")
) 
#' @rdname modify
#' @export
setGeneric( 
  name = "modifications",
  def = function(x,
                 ...) standardGeneric("modifications")
) 

# check functions --------------------------------------------------------------

#' @rdname aggregate
#' @export
setGeneric( 
  name = "hasAggregateData",
  def = function(x,
                 ...) standardGeneric("hasAggregateData")
) 

# Modifier functions -----------------------------------------------------------

#' @rdname RNAmodR-internals
#' @export
setGeneric( 
  name = ".constructModRanges",
  signature = c("range","data"),
  def = function(range,
                 data,
                 modType,
                 scoreFun,
                 source,
                 type) standardGeneric(".constructModRanges")
)
#' @rdname modify
#' @export
setGeneric( 
  name = "modify",
  signature = c("x"),
  def = function(x,
                 force = FALSE) standardGeneric("modify")
)
#' @rdname plotROC
#' @export
setGeneric( 
  name = "plotROC",
  signature = c("x"),
  def = function(x,
                 coord,
                 redo = FALSE,
                 ...)
    standardGeneric("plotROC")
) 

# SequenceData functions -------------------------------------------------------


# Modifier/SequenceData functions ----------------------------------------------

#' @rdname aggregate
#' @export
setGeneric( 
  name = "aggregate",
  signature = c("x"),
  def = function(x,
                 ...)
    standardGeneric("aggregate")
)

# SequenceData/Modifier/ModifierSet functions ----------------------------------

#' @rdname subset
#' @export
setGeneric( 
  name = "subsetByCoord",
  signature = c("x","coord"),
  def = function(x,
                 coord,
                 ...)
    standardGeneric("subsetByCoord")
)
#' @rdname subset
#' @export
setGeneric( 
  name = "subsetByCoord",
  signature = c("x","coord"),
  def = function(x,
                 coord,
                 ...)
    standardGeneric("subsetByCoord")
)
#' @rdname visualizeData
#' @export
setGeneric(
  name = "visualizeDataByCoord",
  signature = c("x","coord"),
  def = function(x,
                 coord,
                 type,
                 window.size = 15L,
                 ...) standardGeneric("visualizeDataByCoord")
)
#' @rdname visualizeData
#' @export
setGeneric(
  name = "visualizeData",
  signature = c("x"),
  def = function(x,
                 name,
                 from = 1L,
                 to = 30L,
                 type,
                 ...) standardGeneric("visualizeData")
)
#' @rdname RNAmodR-internals
#' @export
setGeneric(
  name = ".dataTracks",
  signature = c("x","data","seqdata","sequence"),
  def = function(x,
                 data,
                 seqdata,
                 sequence,
                 args) standardGeneric(".dataTracks")
)

# ModifierSet functions --------------------------------------------------------

#' @rdname compare
#' @export
setGeneric( 
  name = "compare",
  signature = c("x"),
  def = function(x,
                 name,
                 from = 1L,
                 to = 30L,
                 ...)
    standardGeneric("compare")
) 
#' @rdname compare
#' @export
setGeneric( 
  name = "compareByCoord",
  signature = c("x","coord"),
  def = function(x,
                 coord,
                 ...)
    standardGeneric("compareByCoord")
) 
#' @rdname compare
#' @export
setGeneric( 
  name = "plotCompare",
  signature = c("x"),
  def = function(x,
                 name,
                 from = 1L,
                 to = 30L,
                 normalize,
                 ...)
    standardGeneric("plotCompare")
) 
#' @rdname compare
#' @export
setGeneric( 
  name = "plotCompareByCoord",
  signature = c("x","coord"),
  def = function(x,
                 coord,
                 normalize,
                 ...)
    standardGeneric("plotCompareByCoord")
) 
