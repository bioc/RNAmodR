#' @include RNAmodR.R
#' @include Modifier-class.R
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

.norm_viz_args_Modifier <- function(input, x){
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
  args <- c(.norm_viz_args_SequenceData(input, x),
            list(modified.seq = modified.seq,
                 additional.mod = additional.mod))
  args
}

.get_viz_sequence <- function(seq,coord,args,modifications){
  if(args[["modified.seq"]]){
    if(length(modifications) > 0L){
      seq <- combineIntoModstrings(seq,modifications)
    }
  }
  seq
}

# ------------------------------------------------------------------------------

#' @rdname visualizeData
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "Modifier", coord = "GRanges"),
  definition = function(x, coord, type = NA, window.size = 15L, ...) {
    # input check
    coord <- .norm_coord_for_visualization(ranges(x), coord)
    from_to <- .get_viz_from_to_coord(ranges(x), coord, window.size)
    visualizeData(x, name = coord$Parent, from = from_to$from,
                  to = from_to$to, type = type, ...)
  }
)

#' @rdname visualizeData
#' @export
setMethod(
  f = "visualizeData",
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
      st <- .get_viz_sequence_track(x, chromosome, args)
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

#' @rdname visualizeData
#' @export
setMethod(
  f = "getDataTrack",
  signature = signature(x = "Modifier"),
  definition = function(x, name = name, ...) {
    stop("'getDataTrack' needs to be implemented for class '", class(x)[[1]], "'")
  }
)
