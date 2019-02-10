#' @include RNAmodR.R
#' @include Modifier-class.R
NULL

#' @name visualizeData
#' @aliases visualizeData visualizeDataByCoord
#'
#' @title visualizeData
#'
#' @description
#' title
#'
#' @param x a \code{Modifier} or \code{ModifierSet} object.
#' @param coord coordinates of a positions to subset to as a
#' \code{GRanges} object. The Parent column is expected to match the gene or
#' transcript name.
#' @param name Only for \code{visualizeData}: the transcript name
#' @param from Only for \code{visualizeData}: start position
#' @param to Only for \code{visualizeData}: end position
#' @param type the data type of data show as data tracks.
#' @param seqdata \code{TRUE} or \code{FALSE}: whould the sequence data be
#' shown? (default: \code{seqdata = FALSE})
#' @param window.size integer value for the number of positions on the left and
#' right site of the selected positions included in the plotting (default:
#' \code{window.size = 15L})
#' @param ... optional parameters:
#' \itemize{
#' \item{\code{modified.seq}}{\code{TRUE} or \code{FALSE}. Should the sequence
#' shown with modified nucleotide positions? (default:
#' \code{modified.seq = FALSE})}
#' \item{\code{additional.mod}}{other modifications, which should be shown
#' in the annotation and sequence track.}
#' \item{\code{annotation.track.pars}}{Parameters passed onto the A
#' nnotationTrack.}
#' \item{\code{sequence.track.pars}}{Parameters passed onto the SequenceTrack.}
#' \item{\code{data.track.pars}}{Parameters passed onto the DataTrack(s).}
#' }
NULL

.norm_seqdata_show <- function(seqdata){
  if(missing(seqdata) || !assertive::is_a_bool(seqdata)){
    seqdata <- FALSE
  }
  seqdata
}

.norm_type <- function(type){
  if(is.na(type) || !is.character(type)){
    stop("'type' must be a character vector.", call. = FALSE)
  }
  type
}

.norm_viz_args_Modifier <- function(input){
  modified.seq <- FALSE
  additional.mod <- GRanges()
  if(!is.null(input[["modified.seq"]])){
    modified.seq <- input[["modified.seq"]]
    if(!assertive::is_a_bool(modified.seq)){
      stop("'modified.seq' must be a single logical value.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["additional.mod"]])){
    additional.mod <- input[["additional.mod"]]
    if(!is(additional.mod,"GRanges") && !is(additional.mod,"GRangesList")){
      stop("'additional.mod' must be a GRanges or GRangesList object, which is",
           " compatible with combineIntoModstrings().",
           call. = FALSE)
    }
  }
  args <- c(.norm_viz_args_SequenceData(input),
            list(modified.seq = modified.seq,
                 additional.mod = additional.mod))
  args
}

# .norm_addition.mod_for_visualization <- function(additional.mod,
#                                                  coord){
#   if(is(additional.mod,"GRangesList")){
#     additional.mod <- unlist(additional.mod)
#   }
#   additional.mod <- additional.mod[additional.mod$Parent == coord$Parent]
#   if(length(additional.mod) == 0L){
#     return(additional.mod)
#   }
#   GRanges(seqnames = coord$Parent,
#           ranges = ranges(additional.mod),
#           strand = strand(additional.mod),
#           mcols(additional.mod))
# }

# .get_viz_sequence <- function(seq,coord,args,modifications){
#   if(args[["modified.seq"]]){
#     if(length(modifications) > 0L){
#       seq <- combineIntoModstrings(seq,modifications)
#     }
#   }
#   seq
# }


# ------------------------------------------------------------------------------

#' @rdname visualizeData
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "Modifier", coord = "GRanges"),
  definition = function(x, coord, type = NA, seqdata = FALSE, window.size = 15L,
                        ...) {
    # input check
    coord <- .norm_coord_for_visualization(ranges(x), coord)
    from_to <- .get_viz_from_to_coord(ranges(x), coord, window.size)
    visualizeData(x, name = coord$Parent, from = from_to$from,
                  to = from_to$to, type = type, seqdata = seqdata, ...)
  }
)

#' @rdname visualizeData
#' @export
setMethod(
  f = "visualizeData",
  signature = signature(x = "Modifier"),
  definition = function(x, name, from, to, type = NA, seqdata = FALSE, ...) {
    # get plotting arguments
    args <- .norm_viz_args_Modifier(list(...))
    chromosome <- .norm_viz_chromosome(ranges(x), name)
    from_to <- .get_viz_from_to(ranges(x), name, from, to)
    seqdata <- .norm_seqdata_show(seqdata)
    type <- .norm_type(type)
    # get tracks
    atm <- .get_viz_annotation_track(ranges(x), args[["annotation.track.pars"]],
                                     args[["alias"]])
    st <- .get_viz_sequence_track(sequences(x), ranges(x), chromosome,
                                  args[["sequence.track.pars"]])
    dt <- getDataTrack(x, name = name, type = type, ...)
    if(!is.list(dt)){
      dt <- list(dt)
    }
    if(seqdata){
      sdt <- getDataTrack(seqData(x), name = name,...)
      if(!is.list(sdt)){
        sdt <- list(sdt)
      }
      tracks <- c(dt,sdt,list(st,atm))
    } else {
      tracks <- c(dt,list(st,atm))
    }
    # plot tracks
    do.call(Gviz::plotTracks,
            c(list(tracks, from = from_to$from, to = from_to$to,
                   chromosome = chromosome),
              args[["plot.pars"]]))
  }
)

#' @rdname visualizeData
#' @export
setMethod(
  f = "getDataTrack",
  signature = signature(x = "Modifier"),
  definition = function(x, name = name, ...) {
    stop("'getDataTrack' needs to be implemented for class '", class(x)[[1]], "'")
  }
)
