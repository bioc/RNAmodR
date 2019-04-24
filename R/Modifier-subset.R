#' @include RNAmodR.R
#' @include Modifier-class.R
NULL

#' @name subsetByCoord
#' @aliases subsetByCoord labelByCoord
#' 
#' @title Subsetting data from a \code{Modifier} or \code{ModifierSet} object.
#' 
#' @description 
#' With \code{subsetByCoord} data from a \code{Modifier} or \code{ModifierSet}
#' object will be subset to position as defined in \code{coord}. If \code{coord}
#' contains a column \code{mod} and \code{x} of type \code{Modifier}, it will
#' be filtered to identifiers matching the
#' \code{\link[=Modifier-functions]{modType}} of \code{x}. To disable remove
#' the column \code{mod} from \code{coord} or set \code{type = NA}
#' 
#' \code{labelByCoord} functions the same. It will return a
#' \code{SplitDataFrameList}, which matches the dimensions of the aggregated
#' data plus the \code{labels} column.
#' 
#' @param x a \code{Modifier} or \code{ModifierSet} object.
#' @param coord coordinates of position to subset to. Either a \code{GRanges} or
#' a \code{GRangesList} object. For both types the Parent column is expected to
#' match the gene or transcript name.
#' @param ... optional parameters:
#' \itemize{
#' \item{\code{name}} {Limit results to one specific gene or transcript}
#' \item{\code{type}} {the modification type used for subsetting. By default this
#' is derived from the \code{modType(x)}, but it can be overwritten using 
#' \code{type}. It must be a valid shortName for a modification according to
#' \code{shortName(ModRNAString())} and of course present in metadata column 
#' \code{mod} of \code{coord}}
#' \item{flanking} {a single integer value how many flanking position should be
#' included in the subset (default = \code{flanking = 0L}).}
#' \item{merge} {\code{TRUE} or \code{FALSE}: Should the 
#' overlapping selections by merge? This is particular important, if flanking
#' value \code{!= 0L} are set. (default: \code{merge = TRUE}).}
#' \item{\code{perTranscript}} {\code{TRUE} or \code{FALSE}: Should the 
#' positions labeled per transcript and not per chromosome?
#' (default: \code{perTranscript = FALSE}).}
#' }
#' 
#' @return 
#' If 'x' is a
#' \itemize{
#' \item{\code{\link[=Modifier-class]{Modifier}}} {a \code{SplitDataFrameList}
#' with elments per transcript.}
#' \item{\code{\link[=ModifierSet-class]{ModifierSet}}} {a \code{SimpleList} of
#' \code{SplitDataFrameList} with elments per transcript.}
#' }
#' 
#' @examples
#' data(msi,package="RNAmodR")
#' mod <- modifications(msi)
#' coord <- unique(unlist(mod))
#' coord$score <- NULL
#' coord$sd <- NULL
#' subsetByCoord(msi,coord)
NULL

# subsetting Modifier ----------------------------------------------------------

.subset_Modifier_by_GRangesList <- function(x, coord, ...){
  args <- .norm_subset_args(list(...), x)
  if(args[["sequenceData"]]){
    return(subsetByCoord(sequenceData(x), coord, ...))
  }
  # converts everything to a GRangesList
  coord <- .norm_coord(coord, args[["type"]]) 
  data <- getAggregateData(x)
  .subset_SplitDataFrameList_by_GRangesList(data, coord, ...)
}

# subsetting ModifierSet -------------------------------------------------------

.subset_ModifierSet_by_GRangesList <- function(x, coord, ...){
  args <- .norm_subset_args(list(...),x)
  # converts everything to a GRangesList
  coord <- .norm_coord(coord, args[["type"]], args[["merge"]])
  if(args[["sequenceData"]]){
    lapply(x,
           function(z){
             subsetByCoord(sequenceData(z), coord, ...)
           })
  } else {
    lapply(x,
           function(z){
             data <- getAggregateData(z)
             names <- .get_element_names(data, coord, args[["name"]],
                                         args[["type"]])
             data <- data[match(names, names(data))]
             coord <- coord[match(names, names(coord))]
             .perform_subset(data, coord, args[["flanking"]], 
                             args[["perTranscript"]])
           })
  }
}

#' @rdname subsetByCoord
#' @export
setMethod("subsetByCoord",
          signature = c(x = "Modifier", coord = "GRanges"),
          function(x, coord, ...){
            .subset_Modifier_by_GRangesList(x, coord, ...)
          }
)
#' @rdname subsetByCoord
#' @export
setMethod("subsetByCoord",
          signature = c(x = "Modifier", coord = "GRangesList"),
          function(x, coord, ...){
            .subset_Modifier_by_GRangesList(x, coord, ...)
          }
)
#' @rdname subsetByCoord
#' @export
setMethod("subsetByCoord",
          signature = c(x = "ModifierSet", coord = "GRanges"),
          function(x, coord, ...){
            .subset_ModifierSet_by_GRangesList(x, coord, ...)
          }
)
#' @rdname subsetByCoord
#' @export
setMethod("subsetByCoord",
          signature = c(x = "ModifierSet", coord = "GRangesList"),
          function(x, coord, ...){
            .subset_ModifierSet_by_GRangesList(x, coord, ...)
          }
)

################################################################################
# This is used for ROC and shares functionality with subsetting

.label_Modifier_by_GRangesList <- function(x, coord, ...){
  args <- .norm_subset_args(list(...), x)
  # converts everything to a GRangesList
  coord <- .norm_coord(coord, args[["type"]])
  data <- getAggregateData(x)
  .label_SplitDataFrameList_by_GRangesList(data, coord, ...)
}


#' @rdname subsetByCoord
#' @export
setMethod("labelByCoord",
          signature = c(x = "Modifier", coord = "GRanges"),
          function(x, coord, ...){
            .label_Modifier_by_GRangesList(x, coord, ...)
          }
)
#' @rdname subsetByCoord
#' @export
setMethod("labelByCoord",
          signature = c(x = "Modifier", coord = "GRangesList"),
          function(x, coord, ...){
            .label_Modifier_by_GRangesList(x, coord, ...)
          }
)
#' @rdname subsetByCoord
#' @export
setMethod("labelByCoord",
          signature = c(x = "ModifierSet", coord = "GRanges"),
          function(x, coord, ...){
            .label_Modifier_by_GRangesList(x, coord, ...)
          }
)
#' @rdname subsetByCoord
#' @export
setMethod("labelByCoord",
          signature = c(x = "ModifierSet", coord = "GRangesList"),
          function(x, coord, ...){
            .label_Modifier_by_GRangesList(x, coord, ...)
          }
)
