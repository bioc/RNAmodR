#' @include RNAmodR.R
#' @include SequenceDataList-class.R
NULL

#' @rdname visualizeData
#' @export
setMethod(
  f = "getDataTrack",
  signature = signature(x = "SequenceDataList"),
  definition = function(x, name = name, ...) {
    sdts <- lapply(x,
                   function(z){
                     getDataTrack(z, name = name, ...)
                   })
    sdts <- unname(sdts)
    do.call(c,sdts)
  }
)


#' @rdname visualizeData
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "SequenceDataList", coord = "GRanges"),
  definition = function(x, coord, type = NA, window.size = 15L,
                        perTranscript = FALSE, ...) {
    # input check
    coord <- .norm_coord_for_visualization(ranges(x), coord)
    from_to <- .get_viz_from_to_coord(ranges(x), coord, window.size)
    visualizeData(x, name = coord$Parent, from = from_to$from,
                  to = from_to$to, type = type, ...)
  }
)

#' @rdname visualizeData
#' @importFrom Gviz plotTracks
#' @export
setMethod(
  f = "visualizeData",
  signature = signature(x = "SequenceDataList"),
  definition = function(x, name, from, to, perTranscript = FALSE, 
                        ...) {
    # get plotting arguments
    args <- .norm_viz_args_SequenceData(list(...))
    chromosome <- .norm_viz_chromosome(ranges(x), name)
    from_to <- .get_viz_from_to(ranges(x), name, from, to)
    # get tracks
    atm <- .get_viz_annotation_track(ranges(x), args[["annotation.track.pars"]],
                                     args[["alias"]])
    st <- .get_viz_sequence_track(sequences(x), ranges(x), chromosome,
                                  args[["sequence.track.pars"]])
    dts <- getDataTrack(x, name = name, ...)
    tracks <- c(dts,list(st,atm))
    # plot tracks
    do.call(Gviz::plotTracks,
            c(list(tracks, from = from_to$from, to = from_to$to,
                   chromosome = chromosome),
              args[["plot.pars"]]))
  }
)
