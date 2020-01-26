#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

# normalization functions for visualizations -----------------------------------

.norm_coord_for_visualization <- function(ranges, coord){
  if(length(coord) > 1L){
    stop("'coord' must be contain a single range.",
         call. = FALSE)
  }
  if(is.null(coord$Parent)){
    stop("'coord' must contain a metadata column named 'Parent'.",
         call. = FALSE)
  }
  if(!(coord$Parent %in% names(ranges))){
    stop("Transcript identifier '",coord$Parent,"' not found in data of 'x'.",
         call. = FALSE)
  }
  coord
}

.norm_viz_windows.size <- function(window.size){
  if(!is.integer(window.size) || length(window.size) > 1L){
    stop("'window.size' must be a single integer value.",
         call. = FALSE)
  }
  window.size
}

.get_viz_from_to_coord <- function(ranges, coord, window.size){
  window.size <- .norm_viz_windows.size(window.size)
  start <- start(coord) - window.size + 1L
  end <- end(coord) + window.size
  pos.min <- as.integer(min(start(ranges[[as.character(coord$Parent)]])))
  pos.max <- as.integer(max(end(ranges[[as.character(coord$Parent)]])))
  if(start < pos.min){
    start <- pos.min
  }
  if(end > pos.max){
    end <- pos.max
  }
  if(start > end){
    end <- pos.max
  }
  list(from = start,
       to = end)
}

.get_viz_from_to <- function(ranges, name, from, to){
  start <- from
  end <- to
  pos.min <- as.integer(min(start(ranges[[name]])))
  pos.max <- as.integer(max(end(ranges[[name]])))
  if(start < pos.min){
    start <- pos.min
  }
  if(end > pos.max){
    end <- pos.max
  }
  if(start > end){
    end <- pos.max
  }
  list(from = start,
       to = end)
}

.norm_viz_colour <- function(colour, type = NA){
  if(!is.character(colour) || any(!.are_colours(colour))){
    stop("'colour' must be a character vector and contain valid colours, which",
         "can be interpreted by col2rgb().",
         call. = FALSE)
  }
  if(length(type) > 1L || !anyNA(type)){
    if(length(colour) != 1){
      if(is.null(names(colour)) || !all(names(colour) %in% type)  ){
           stop("'colour' must be a named character vector parallel to 'type' ",
                "or of length = 1.",
                call. = FALSE)
      }
    }
    if(length(colour) == 1){
      colour <- rep(colour,length(type))
      names(colour) <- type
    }
  }
  colour
}

.norm_viz_chromosome <- function(ranges, name){
  if(!(name %in% names(ranges))){
    stop("Transcript name '",name,"' not found in 'x'", call. = FALSE)
  }
  as.character(unique(unlist(seqnames(ranges[name]),use.names=FALSE)))
}

.norm_viz_args_SequenceData <- function(input, x){
  sequence.track.pars <- list(add53 = FALSE)
  annotation.track.pars <- list()
  plot.pars <- list()
  if(!is.null(input[["sequence.track.pars"]])){
    sequence.track.pars <- input[["sequence.track.pars"]]
    sequence.track.pars <- 
      sequence.track.pars[!(names(plot.pars) %in% c("sequence"))]
  }
  if(!is.null(input[["annotation.track.pars"]])){
    annotation.track.pars <- input[["annotation.track.pars"]]
    annotation.track.pars <- 
      annotation.track.pars[!(names(plot.pars) %in% c("range"))]
  }
  if(!is.null(input[["plot.pars"]])){
    plot.pars <- input[["plot.pars"]]
    plot.pars <- 
      plot.pars[!(names(plot.pars) %in% c("tracks","from","to","chromosome"))]
  }
  args <- c(.norm_alias(input, x),
            list(sequence.track.pars = sequence.track.pars,
                 annotation.track.pars = annotation.track.pars,
                 plot.pars = plot.pars))
  args
}

# track loading wraps ----------------------------------------------------------

#' @importFrom Gviz AnnotationTrack
.get_viz_annotation_track <- function(x, args){
  ranges <- ranges(x)
  alias <- args[["alias"]]
  args <- args[["annotation.track.pars"]]
  d <- mcols(ranges@unlistData)
  if(!is.null(d$type) && is.null(d$feature)){
    d$feature <- d$type
  } else {
    d$feature <- vector(mode = "character", length = length(ranges))
  }
  if(is.null(d$group)){
    group <- vector(mode = "character", length = length(ranges))
    if(!is.null(alias)){
      m <- match(as.character(alias$tx_id),names(ranges))
      group[m[!is.na(m)]] <- as.character(alias$name)[!is.na(m)]
    }
    f <- group == ""
    group[f] <- names(ranges)[f]
    d$group <- rep(group,lengths(ranges))
  }
  if(is.null(d$id)){
    d$id <- d$exon_name
  }
  mcols(ranges@unlistData) <- d
  at <- Gviz::AnnotationTrack(range = unlist(ranges),
                              collapse = TRUE,
                              collapse = TRUE,
                              background.title = "#FFFFFF",
                              fontcolor.legend = "#000000",
                              featureAnnotation = "group",
                              fontsize.group = 10)
  if(!is.null(args$annotation.track.pars)){
    Gviz::displayPars(at) <- args$annotation.track.pars
  }
  at
}

.stitch_chromosome <- function(seq, ranges, chromosome){
  ranges <- ranges[seqnames(ranges) == chromosome]
  ranges <- ranges[!vapply(ranges,function(r){length(r) == 0L},logical(1))]
  if(length(ranges) == 0L){
    stop("No ranges with seqnames = '",chromosome,"' found.")
  }
  ranges_names <- names(ranges)
  seq <- seq[names(seq) %in% ranges_names]
  if(length(seq) == 0L){
    stop("No sequences for seqnames = '",chromosome,"' found.")
  }
  if(any(sum(width(ranges))!= width(seq))){
    stop("width() or sequences and ranges does not match.")
  }
  unlisted_ranges <- unlist(ranges)
  hits <- findOverlaps(unlisted_ranges)
  # if overlapping ranges exist, merge em
  if(length(hits) > length(unlisted_ranges)){
    stop("")
  }
  # get gaps in ranges
  GenomeInfoDb::seqlevels(unlisted_ranges) <- 
    GenomeInfoDb::seqlevelsInUse(unlisted_ranges)
  gaps <- gaps(unlisted_ranges)
  if(length(gaps) == 0L){
    names(seq) <- chromosome
    return(seq)
  } 
  gaps <- gaps[end(gaps) <= max(end(unlisted_ranges))]
  if(length(gaps) == 0L){
    names(seq) <- chromosome
    return(seq)
  }
  # get N sequence for gaps
  FUN <- match.fun(class(seq))
  N <- FUN(rep(paste(rep("N",width(gaps)),collapse = ""),length(gaps)))
  seq <- relist(unlist(seq),
                IRanges::PartitioningByWidth(width(unlisted_ranges)))
  # assemble result
  ans <- pc(N,seq)
  ans <- FUN(unlist(ans))
  names(ans) <- chromosome
  ans
}

.get_viz_sequence_track <- function(seq, ranges, chromosome, args){
  args <- args[["sequence.track.pars"]]
  FUN <- function(trackClass, seqClass, seq, args){
    set <- do.call(seqClass,list(seq))
    track <- do.call(trackClass,
                     c(list(sequence = set),
                       args[["sequence.track.pars"]]))
    track
  }
  if(!is(seq,"DNAStringSet") &&  !is(seq,"RNAStringSet") && 
     !is(seq,"ModRNAStringSet") && !is(seq,"ModDNAStringSet")){
    stop("Invalid sequence type '",class(seq),"'. sequences(x) must be a ",
         "RNA/ModRNA/DNA/ModDNA*StringSet.",
         call. = FALSE)
  }
  # reconstruct the chromosomal sequences for plotting
  seq <- .stitch_chromosome(seq, ranges, chromosome)
  if(is(seq,"RNAStringSet")){
    st <- FUN("RNASequenceTrack","RNAStringSet", seq, args)
  } else if(is(seq,"ModRNAStringSet")){
    st <- FUN("ModRNASequenceTrack","ModRNAStringSet", seq, args)
  } else if(is(seq,"DNAStringSet")){
    st <- FUN("SequenceTrack","DNAStringSet", seq, args)
  } else if(is(seq,"ModDNAStringSet")){
    st <- FUN("ModDNASequenceTrack","ModDNAStringSet", seq, args)
  } else {
    stop("")
  }
  st
}


# helper functions for visualization -------------------------------------------

.get_data_for_visualization <- function(x, name){
  ranges <- ranges(x)
  if(!(name %in% names(ranges))){
    stop("Transcript name '",name,"' not found in 'x'", call. = FALSE)
  }
  ranges <- ranges[name]
  if(is(x,"Modifier")){
    data <- getAggregateData(x)[name]
  } else if(is(x,"SequenceData")){
    x <- x[name]
    data <- aggregate(x)
  } else {
    stop("")
  }
  strand_u <- .get_strand_u_GRangesList(ranges)
  seqs <- .seqs_rl(ranges)
  seqs[strand_u == "-"] <- rev(seqs[strand_u == "-"])
  ranges_seqnames <- .seqnames_rl(ranges)
  strands <- .strands_rl(ranges)
  ans <- GenomicRanges::GRanges(seqnames = unlist(ranges_seqnames),
                                ranges = IRanges::IRanges(start = unlist(seqs),
                                                          width = 1),
                                strand = unlist(strands),
                                unlist(data, use.names = FALSE))
  ans <- relist(ans, IRanges::PartitioningByEnd(data))
  metadata(ans) <- metadata(ranges)
  ans
}

################################################################################

#' @rdname plotData
#' @export
setMethod(
  f = "plotDataByCoord",
  signature = signature(x = "SequenceData", coord = "GRanges"),
  definition = function(x, coord, type = NA, window.size = 15L, ...) {
    # input check
    coord <- .norm_coord_for_visualization(ranges(x), coord)
    from_to <- .get_viz_from_to_coord(ranges(x), coord, window.size)
    plotData(x, name = coord$Parent, from = from_to$from,
                  to = from_to$to, type = type, ...)
  }
)

#' @rdname plotData
#' @importFrom Gviz plotTracks
#' @export
setMethod(
  f = "plotData",
  signature = signature(x = "SequenceData"),
  definition = function(x, name, from, to, perTranscript = FALSE, 
                        showSequence = TRUE, showAnnotation = FALSE, ...) {
    # get plotting arguments
    args <- .norm_viz_args_SequenceData(list(...), x)
    chromosome <- .norm_viz_chromosome(ranges(x), name)
    from_to <- .get_viz_from_to(ranges(x), name, from, to)
    showSequence <- .norm_show_argument(showSequence, TRUE)
    showAnnotation <- .norm_show_argument(showAnnotation, FALSE)
    # get tracks
    atm <- NULL
    st <- NULL
    if(showAnnotation){
      atm <- .get_viz_annotation_track(x,args)
    }
    if(showSequence){
      st <- .get_viz_sequence_track(sequences(x)[name], ranges(x)[name], 
                                    chromosome, args)
    }
    dt <- getDataTrack(x, name = name, ...)
    if(!is.list(dt)){
      dt <- list(dt)
    }
    tracks <- c(dt,list(st,atm))
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
  signature = signature(x = "SequenceData"),
  definition = function(x, name = name, ...) {
    stop("This functions needs to be implemented by '",class(x),"'.",
         call. = FALSE)
  }
)
