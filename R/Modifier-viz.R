#' @include RNAmodR.R
#' @include Modifier-class.R
NULL

#' @name plotData
#' @aliases plotData plotDataByCoord getDataTrack
#' 
#' @title Visualizing data data from a \code{SequenceData}, 
#' \code{SequenceDataSet}, \code{SequenceDataList}, \code{Modifier} or 
#' \code{ModifierSet} object.
#' 
#' @description 
#' With the \code{plotData} and \code{plotDataByCoord} functions data
#' from a \code{SequenceData}, \code{SequenceDataSet}, \code{SequenceDataList},
#' \code{Modifier} or \code{ModifierSet} object can be visualized.
#' 
#' Internally the functionality of the \code{Gviz} package is used. For each
#' \code{SequenceData} and \code{Modifier} class the \code{getDataTrack} is
#' implemented returning a \code{\link[Gviz:DataTrack-class]{DataTrack}} object
#' from the \code{Gviz} package.
#' 
#' Positions to be visualized are selected by defining a genomic coordinate,
#' for which \code{x} has to contain data.
#' 
#' @param x a \code{SequenceData}, \code{SequenceDataSet},
#'   \code{SequenceDataList}, \code{Modifier} or \code{ModifierSet} object.
#' @param coord coordinates of a positions to subset to as a 
#' \code{GRanges} object. The 'Parent' column is expected to match the 
#' transcript name.
#' @param name Only for \code{plotData}: the transcript name
#' @param from Only for \code{plotData}: start position
#' @param to Only for \code{plotData}: end position
#' @param type the data type of data show as data tracks.
#' @param showSequenceData \code{TRUE} or \code{FALSE}: should the sequence data
#' be shown? (default: \code{seqdata = FALSE})
#' @param showSequence \code{TRUE} or \code{FALSE}: should a sequence track be 
#' shown? (default: \code{seqdata = TRUE})
#' @param showAnnotation \code{TRUE} or \code{FALSE}: should a annotation track 
#' be shown? (default: \code{seqdata = FALSE})
#' @param window.size integer value for the number of positions on the left and 
#' right site of the selected positions included in the plotting (default: 
#' \code{window.size = 15L})
#' @param perTranscript \code{TRUE} or \code{FALSE}: Should the positions shown
#' per transcript? (default: \code{perTranscript = FALSE})
#' @param ... optional parameters:
#' \itemize{
#' \item{\code{modified.seq}} {\code{TRUE} or \code{FALSE}. Should the sequence 
#' shown with modified nucleotide positions? (default: 
#' \code{modified.seq = FALSE})}
#' \item{\code{additional.mod}} {other modifications, which should be shown
#' in the annotation and sequence track. The must be a \code{GRanges} compatible
#' with \code{\link[Modstrings:separate]{combineIntoModstrings}}.}
#' \item{\code{annotation.track.pars}} {Parameters passed onto the 
#' \code{\link[Gviz:AnnotationTrack-class]{AnnotationTrack}}.}
#' \item{\code{sequence.track.pars}} {Parameters passed onto the 
#' \code{\link[Gviz:SequenceTrack-class]{SequenceTrack}}.}
#' }
#' 
#' @return a plot send to the active graphic device
#' 
#' @examples 
#' data(msi,package="RNAmodR")
#' plotData(msi[[1]], "2", from = 10L, to = 45L)
#' \dontrun{
#' plotData(msi, "2", from = 10L, to = 45L)
#' }
NULL

.norm_show_argument <- function(show_arg, default = FALSE){
  if(missing(show_arg) || !assertive::is_a_bool(show_arg)){
    show_arg <- default
  }
  show_arg
}

.norm_score_type <- function(type, colnames = NA, multiple = FALSE){
  if(missing(type) && !anyNA(colnames)){
    if(multiple){
      type <- colnames
    } else {
      type <- colnames[1]
    }
  } else if(missing(type)) {
    stop("'type' is missing.", call. = FALSE)
  }
  if(is.na(type) || !is.character(type)){
    stop("'type' must be a character vector.", call. = FALSE)
  }
  if(!anyNA(colnames)){
    if(any(!(type %in% colnames))){
      stop("'type' was not found in data.", call. = FALSE)
    }
  }
  type
}

.viz_Modifier_settings <- data.frame(
  variable = c("modified.seq",
               "additional.mod"),
  testFUN = c(".is_a_bool",
              ".is_not_GRanges_or_GRangesList"),
  errorValue = c(FALSE,
                 TRUE),
  errorMessage = c("'modified.seq' must be a single logical value.",
                   "'additional.mod' must be a GRanges or GRangesList object, which is compatible with combineIntoModstrings()."),
  stringsAsFactors = FALSE)
.norm_viz_args_Modifier <- function(input, x){
  modified.seq <- FALSE
  additional.mod <- GRanges()
  args <- .norm_settings(input, .viz_Modifier_settings, modified.seq,
                         additional.mod)
  args <- c(.norm_viz_args_SequenceData(input, x),
            args)
  args
}

.get_viz_sequence <- function(x,args){
  seq <- sequences(x)
  if(args[["modified.seq"]]){
    mod <- modifications(x)
    if(is(mod,"GRangesList")){
      mod <- unlist(mod)
    }
    mcols(mod) <- mcols(mod)[,c("mod","Parent")]
    mod <- unique(mod)
    add.mod <- args[["additional.mod"]]
    if(length(add.mod) > 0L){
      mcols(add.mod) <- mcols(add.mod)[,c("mod","Parent")]
      mod <- c(mod,add.mod)
      mod <- unique(mod)
    }
    mod <- .rebase_seqnames(mod, mod$Parent)
    if(length(modifications) > 0L){
      seq <- combineIntoModstrings(seq,mod)
    }
  }
  seq
}

# ------------------------------------------------------------------------------

#' @rdname plotData
setMethod(
  f = "plotDataByCoord",
  signature = signature(x = "Modifier", coord = "GRanges"),
  definition = function(x, coord, type = NA, window.size = 15L, ...) {
    # input check
    coord <- .norm_coord_for_visualization(ranges(x), coord)
    from_to <- .get_viz_from_to_coord(ranges(x), coord, window.size)
    plotData(x, name = coord$Parent, from = from_to$from,
                  to = from_to$to, type = type, ...)
  }
)

#' @rdname plotData
#' @export
setMethod(
  f = "plotData",
  signature = signature(x = "Modifier"),
  definition = function(x, name, from, to, type = NA, showSequenceData = FALSE, 
                        showSequence = TRUE, showAnnotation = FALSE, ...) {
    # get plotting arguments
    args <- .norm_viz_args_Modifier(list(...), x)
    chromosome <- .norm_viz_chromosome(ranges(x), name)
    from_to <- .get_viz_from_to(ranges(x), name, from, to)
    showSequenceData <- .norm_show_argument(showSequenceData, FALSE)
    showSequence <- .norm_show_argument(showSequence, TRUE)
    showAnnotation <- .norm_show_argument(showAnnotation, FALSE)
    type <- .norm_score_type(type)
    # get tracks
    atm <- NULL
    st <- NULL
    if(showAnnotation){
      atm <- .get_viz_annotation_track(x, args)
    }
    if(showSequence){
      seq <- .get_viz_sequence(x, args)[name]
      st <- .get_viz_sequence_track(seq, ranges(x)[name], chromosome, args)
    }
    dt <- getDataTrack(x, name = name, type = type, ...)
    if(!is.list(dt)){
      dt <- list(dt)
    }
    if(showSequenceData){
      sdt <- getDataTrack(sequenceData(x), name = name,...)
      if(!is.list(sdt)){
        sdt <- list(sdt)
      }
      tracks <- c(dt,sdt,list(st,atm))
    } else {
      tracks <- c(dt,list(st,atm))
    }
    # plot tracks
    tracks <- tracks[!vapply(tracks, is.null, logical(1))]
    do.call(Gviz::plotTracks,
            c(list(tracks, from = from_to$from, to = from_to$to,
                   chromosome = chromosome),
              args[["plot.pars"]]))
  }
)

#' @rdname plotData
#' @export
setMethod(
  f = "getDataTrack",
  signature = signature(x = "Modifier"),
  definition = function(x, name = name, ...) {
    stop("'getDataTrack' needs to be implemented for class '", class(x)[[1]], "'")
  }
)
