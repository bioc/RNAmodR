#' @include RNAmodR.R
#' @include Modifier-viz.R
NULL

#' @importFrom grDevices col2rgb
.is_colour <- function(x) {
  vapply(x,
         function(z) {
           tryCatch(is.matrix(grDevices::col2rgb(z)),
                    error = function(e) FALSE)
         },
         logical(1))
}

.norm_viz_args_ModifierSet <- function(input){
  colours <- NA
  if(!is.null(input[["colours"]])){
    colours <- input[["colours"]]
    if(!is.character(colours) || any(!.is_colour(colours))){
      stop("'colours' must be valid colour representation, which can be ",
           "interpreted by col2rgb().",
           call. = FALSE)
    }
  }
  args <- c(.norm_viz_args_Modifier(input),
            list(colours = colours))
  args
}

#' @importFrom RColorBrewer brewer.pal
.norm_viz_colours <- function(x ,colours){
  n <- length(x)
  if(!is.na(colours)){
    if(n != length(colours)){
      stop("'colours' must have the same length as 'x'.", call. = FALSE)
    }
    return(colours)
  }
  # brewer.pal needs n >=3
  ans <- RColorBrewer::brewer.pal(max(3,n),"Set1")
  ans[seq.int(1,n)]
}

.add_viz_colours <- function(dts, colours){
  dts <- mapply(
    function(dt, colour){
      if(is.list(dt)){
        dt <- lapply(dt,
                     function(t){
                       Gviz::displayPars(t)$col <- colour
                       Gviz::displayPars(t)$fill <- colour
                       t
                     })
      } else {
        Gviz::displayPars(dt)$col <- colour
        Gviz::displayPars(dt)$fill <- colour
      }
      dt
    },
    dts,
    colours)
  dts
}

.add_viz_names <- function(dts, names){
  dts <- mapply(
    function(dt, name){
      if(is.list(dt)){
        dt <- lapply(dt,
                     function(t){
                       t@name <- paste0(name,"\n",t@name)
                       t
                     })
      } else {
        dt@name <- paste0(name,"\n",dt@name)
      }
      dt
    },
    dts,
    names)
  dts
}

#' @rdname visualizeData
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "ModifierSet",
                        coord = "GRanges"),
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
  signature = signature(x = "ModifierSet"),
  definition = function(x, name, from, to, type = NA, seqdata = FALSE, ...) {
    browser()
    # get plotting arguments
    args <- .norm_viz_args_ModifierSet(list(...))
    chromosome <- .norm_viz_chromosome(ranges(x), name)
    from_to <- .get_viz_from_to(ranges(x), name, from, to)
    colours <- .norm_viz_colours(x, args[["colours"]])
    # get tracks
    atm <- .get_viz_annotation_track(ranges(x), args[["annotation.track.pars"]],
                                     args[["alias"]])
    st <- .get_viz_sequence_track(sequences(x), ranges(x), chromosome,
                                  args[["sequence.track.pars"]])
    dts <- lapply(x, getDataTrack, ...)
    dts <- .add_viz_ylim(dts, chromosome, from_to, args[["ylim"]])
    dts <- .add_viz_colours(dts, colours)
    dts <- .add_viz_names(dts, names(x))
    if(seqdata){
      sdts <- lapply(x, function(z){getDataTrack(seqData(z), ...)})
      sdts <- .add_viz_ylim(sdts, chromosome, from_to, args[["ylim"]])
      sdts <- .add_viz_names(sdts, names(x))
      tracks <- c(do.call(c,dts), do.call(c,sdts), list(st,atm))
    } else {
      tracks <- c(do.call(c,dts), list(st,atm))
    }
    # plot tracks
    do.call(Gviz::plotTracks,
            c(list(tracks, from = from_to$from, to = from_to$to,
                   chromosome = chromosome),
              args[["plot.pars"]]))
  }
)
