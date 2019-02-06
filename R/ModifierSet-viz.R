#' @include RNAmodR.R
#' @include Modifier-viz.R
NULL

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

.get_viz_window_ModifierSet <- function(data, coord, window.size){
  window.size <- .norm_viz_windows.size(window.size)
  start <- start(coord) - window.size
  end <- end(coord) + window.size
  pos.min <- as.integer(min(unique(unlist(lapply(data,start)))))
  pos.max <- as.integer(max(unique(unlist(lapply(data,end)))))
  if(start < pos.min){
    start <- pos.min
  }
  if(end > pos.max){
    end <- pos.max
  }
  list(start = start,
       end = end)
}

.norm_viz_ylim <- function(data){
  max <- max(vapply(data,
                function(d){
                  max(mcols(d)[,1])
                },
                numeric(1)))
  c(0,max)
}

#' @importFrom RColorBrewer brewer.pal
.norm_viz_colours <- function(data,colours){
  n <- length(data)
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

.norm_data_track_types <- function(dtl,coordValues,args){
  type <- args[["data.track.pars"]][["type"]]
  if(is.null(type)){
    f <- vapply(lapply(lapply(dtl,displayPars),"[[","type"),
                "==",logical(1),"histogram")
    if((coordValues$end - coordValues$start) > 60){
      dtl[f] <- lapply(dtl[f],
                       function(dt){
                         displayPars(dt)$type <- "h"
                         displayPars(dt)$lwd <- 3L
                         dt
                       })
    }
  }
  dtl
}

#' @rdname visualizeData
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "ModifierSet",
                        coord = "GRanges"),
  definition = function(x, coord, type = NA, seqdata = FALSE, window.size = 15L,
                        ...) {
    requireNamespace("Gviz")
    # input check
    seqdata_show <- .norm_seqdata_show(seqdata)
    args <- .norm_viz_args_ModifierSet(list(...))
    coord <- .norm_coord_for_visualization(coord)
    # get plotting data
    modAnnotation <- .norm_viz_mod_annotation(args[["additional.mod"]], coord)
    seq <- .get_viz_sequence(sequences(x)[[coord$Parent]], coord, args,
                             modAnnotation)
    range <- .norm_viz_range(ranges(x)[[coord$Parent]])
    chromosome <- .norm_viz_chromosome(range)
    data <- lapply(x,
                   function(z){
                     .norm_data_for_visualization(
                       aggregateData(z)[[coord$Parent]],
                       type,
                       chromosome)
                   })
    seqdata <- lapply(x,
                      function(z){
                        .norm_seqdata_for_visualization(
                          aggregate(seqData(z)[coord$Parent])[[1]],
                          chromosome)
                      })
    coordValues <- .get_viz_window_ModifierSet(data, coord, window.size)
    # special stuff for data of a ModifierSet
    args[["ylim"]] <- .norm_viz_ylim(data)
    colours <- .norm_viz_colours(data,args[["colours"]])
    args[["colours"]] <- NULL
    # get tracks
    st <- .get_viz_sequence_track(seq,chromosome,args)
    atm <- .get_viz_annotation_track(modAnnotation,chromosome,args)
    dtl <- lapply(seq_along(x),
                 function(i){
                   a <- args
                   a[["colour"]] <- colours[i]
                   dt <- .dataTracks(x[[i]], data = data[[i]],
                                     seqdata = seqdata[[i]], sequence = seq,
                                     args = a)
                   if(seqdata_show){
                     return(dt[["seqdata"]])
                   }
                   dt[[1]]
                 })
    # add sample names
    dtl <- mapply(function(dt,name){
                    names(dt) <- paste0(name, "\n", names(dt))
                    dt
                  },
                  dtl,
                  names(x),
                  SIMPLIFY = FALSE)
    # based in the width change from histogram to h type of the DataTrack
    dtl <- .norm_data_track_types(dtl,coordValues,args)
    # plot tracks
    plotTracks(c(dtl, list(st,atm)), from = coordValues$start, 
               to = coordValues$end, chromosome = chromosome)
  }
)

#' @rdname visualizeData
#' @export
setMethod(
  f = "visualizeData",
  signature = signature(x = "ModifierSet"),
  definition = function(x, name, from, to, type = NA, seqdata = FALSE, ...) {
    coord <- .create_coord_for_visualization(x, name, from, to)
    visualizeDataByCoord(x, coord, type = type, seqdata = seqdata, 
                         window.size = 0L, ...)
  }
)
