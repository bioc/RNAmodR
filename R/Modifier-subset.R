#' @include RNAmodR.R
#' @include Modifier-class.R
#' @include settings.R
NULL

#' @name subsetByCoord
#' @aliases subsetByCoord labelByCoord
#' 
#' @title Subsetting data from a \code{SequenceData}, \code{SequenceDataSet},
#' \code{SequenceDataList}, \code{Modifier} or \code{ModifierSet} object.
#' 
#' @description 
#' With the \code{subsetByCoord} function data from a \code{SequenceData},
#' \code{SequenceDataSet}, \code{SequenceDataList}, \code{Modifier} or
#' \code{ModifierSet} object can be subset to positions as defined in
#' \code{coord}.
#' 
#' If \code{coord} contains a column \code{mod} and \code{x} is a
#' \code{Modifier} object, it will be filtered to identifiers matching the
#' \code{\link[=Modifier-functions]{modType}} of \code{x}. To disable this
#' behaviour remove the column \code{mod} from \code{coord} or set \code{type =
#' NA}
#' 
#' \code{labelByCoord} functions similarly. It will return a
#' \code{SplitDataFrameList}, which matches the dimensions of the aggregated
#' data plus the \code{labels} column, which contains logical values to indicate
#' selected positions.
#' 
#' @param x a \code{SequenceData}, \code{SequenceDataSet},
#' \code{SequenceDataList}, \code{Modifier} or \code{ModifierSet} object.
#' @param coord coordinates of position to subset to. Either a \code{GRanges} or
#' a \code{GRangesList} object. For both types the 'Parent' column is expected to
#' match the transcript name.
#' @param name Optional: Limit results to one specific transcript.
#' @param pos Optional: Limit results to a specific position.
#' @param ... Optional parameters:
#' \itemize{
#' \item{\code{type}:} {the modification type used for subsetting. By default 
#' this is derived from the \code{modType(x)}, but it can be overwritten using 
#' \code{type}. It must be a valid shortName for a modification according to
#' \code{shortName(ModRNAString())} or \code{shortName(ModDNAString())} 
#' (depending on the type of Modifier class) and of course be present in 
#' metadata column \code{mod} of \code{coord}. To disable subsetting based on 
#' type, set \code{type = NA}.}
#' \item{\code{flanking}:} {a single integer value to select how many flanking
#' position should be included in the subset (default: \code{flanking = 0L}).}
#' \item{\code{merge}:} {\code{TRUE} or \code{FALSE}: Should the 
#' overlapping selections be merged? This is particular important, if flanking
#' value \code{!= 0L} are set. (default: \code{merge = TRUE}).}
#' \item{\code{perTranscript}:} {\code{TRUE} or \code{FALSE}: Should the 
#' positions labeled per transcript and not per chromosome?
#' (default: \code{perTranscript = FALSE}).}
#' }
#' 
#' @return 
#' If 'x' is a
#' \itemize{
#' \item{\code{\link[=SequenceData-class]{SequenceData}} or 
#' \code{\link[=Modifier-class]{Modifier}}:} {a \code{SplitDataFrameList}
#' with elments per transcript.}
#' \item{\code{\link[=SequenceDataSet-class]{SequenceDataSet}},
#' \code{\link[=SequenceDataList-class]{SequenceDataList}} or
#' \code{\link[=ModifierSet-class]{ModifierSet}}:} {a \code{SimpleList} of
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
             element_names <- .get_element_names(data, coord, args[["name"]],
                                         args[["type"]])
             data <- data[match(element_names, names(data))]
             coord <- coord[match(element_names, names(coord))]
             .perform_subset(data, coord, args[["flanking"]], 
                             args[["perTranscript"]])
           })
  }
}

#' @rdname RNAmodR-internals
setMethod("subset",
          signature = c(x = "Modifier"),
          function(x, name, pos = 1L, ...){
            coord <- .construct_coord_from_name_from_to(x, name, pos)
            .subset_Modifier_by_GRangesList(x, coord, ...)
          }
)
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
setMethod("subset",
          signature = c(x = "ModifierSet"),
          function(x, name, pos = 1L, ...){
            coord <- .construct_coord_from_name_from_to(x, name, pos)
            .subset_ModifierSet_by_GRangesList(x, coord, ...)
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
