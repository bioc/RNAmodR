#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

# This might change
.norm_viz_args_SequenceData <- .norm_viz_args_Modifier

#' @rdname visualizeData
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "SequenceData",
                        coord = "GRanges"),
  definition = function(x,
                        coord,
                        type = NA,
                        window.size = 15L,
                        ...) {
    requireNamespace("Gviz")
    # input check
    args <- .norm_viz_args_SequenceData(list(...))
    coord <- .norm_coord_for_visualization(coord)
    # get plotting data
    seq <- .get_viz_sequence(sequences(x)[[coord$Parent]],
                             coord,
                             args)
    range <- .norm_viz_range(ranges(x)[[coord$Parent]])
    chromosome <- .norm_viz_chromosome(range)
    seqdata <- .norm_seqdata_for_visualization(
      aggregate(x[coord$Parent])[[1]],
      chromosome)
    coordValues <- .get_viz_window(seqdata,coord,window.size)
    modAnnotation <- .norm_viz_mod_annotation(args[["additional.mod"]],
                                              coord)
    # get tracks
    st <- .get_viz_sequence_track(seq,
                                  chromosome,
                                  args[["sequence.track.pars"]])
    atm <- .get_viz_annotation_track(modAnnotation,
                                     chromosome,
                                     args[["annotation.track.pars"]])
    dt <- .dataTracks(x,
                      seqdata = seqdata,
                      sequence = seq,
                      args = args)
    if(!is.list(dt)){
      dt <- list(dt)
    }
    # plot tracks
    plotTracks(c(dt,
                 list(st,atm)),
               from = coordValues$start,
               to = coordValues$end,
               chromosome = chromosome)
  }
)

#' @rdname visualizeData
#' @export
setMethod(
  f = "visualizeData",
  signature = signature(x = "SequenceData"),
  definition = function(x,
                        name,
                        from,
                        to,
                        type = NA,
                        ...) {
    coord <- .create_coord_for_visualization(x,
                                             name,
                                             from,
                                             to)
    visualizeDataByCoord(x,
                         coord,
                         type = type,
                         window.size = 0L,
                         ...)
  }
)

#' @name RNAmodR-internals
#' @export
setMethod(
  f = ".dataTracks",
  signature = signature(x = "SequenceData",
                        seqdata = "ANY",
                        sequence = "ANY"),
  definition = function(x,
                        seqdata,
                        sequence,
                        args) {
    stop(".dataTracks needs to be implemented for class '",class(x)[[1]],
         "'")
  }
)
