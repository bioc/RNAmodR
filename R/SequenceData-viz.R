#' @include RNAmodR.R
#' @include SequenceData-class.R
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

.norm_viz_name <- function(name){
  if(missing(name)){
    name <- NULL
  }
  name
}

.get_viz_from_to_coord <- function(ranges, coord, window.size){
  window.size <- .norm_viz_windows.size(window.size)
  start <- start(coord) - window.size
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


#' @importFrom grDevices col2rgb
.is_colour <- function(x) {
  vapply(x,
         function(z) {
           tryCatch(is.matrix(grDevices::col2rgb(z)),
                    error = function(e) FALSE)
         },
         logical(1))
}

.norm_viz_colour <- function(colour, type = NA){
  if(!is.character(colour) || !.is_colour(colour)){
    stop("'colour' must be a character vector and contain valid colours, which",
         "can be interpreted by col2rgb().",
         call. = FALSE)
  }
  if(!is.na(type)){
    if(length(colour) != 1 && !all(names(colour) %in% type)){
      stop("'colour' must be a named character vector parallel to 'type' ",
           "or of length = 1.",
           call. = FALSE)
    }
    if(length(colour) == 1){
      colour <- rep(colour,length(type))
      names(colour) <- type
    }
  }
  colour
}

.norm_viz_args_SequenceData <- function(input){
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
  args <- c(.norm_alias(input),
            list(sequence.track.pars = sequence.track.pars,
                 annotation.track.pars = annotation.track.pars,
                 plot.pars = plot.pars))
  args
}

.norm_viz_chromosome <- function(ranges, name){
  if(!(name %in% names(ranges))){
    stop("Transcript name '",name,"' not found in 'x'", call. = FALSE)
  }
  as.character(seqnames(ranges[name]))
}

# track loading wraps ----------------------------------------------------------

#' @importFrom Gviz AnnotationTrack
.get_viz_annotation_track <- function(ranges, args, alias){
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
  names <- names(ranges)
  seq <- seq[names(seq) %in% names]
  ranges <- unlist(ranges)
  hits <- findOverlaps(ranges)
  if(length(hits) > length(ranges)){
    h <- split(subjectHits(hits),
                      queryHits(hits))
    f <- lengths(h) > 1L
    ranges <- c(list(ranges[!f]),
                lapply(h[f],
                       function(i){
                         r <- ranges[i]
                         GenomicRanges::GRanges(seqnames = unique(seqnames(r)),
                                                ranges = IRanges::IRanges(min(start(r)),
                                                                          max(end(r))),
                                                strand = "+")
                       }))
    ranges <- ranges[!duplicated(ranges)]
    ranges <- unlist(GRangesList(ranges))
  }
  max_end <- max(end(ranges))
  starts <- c(1L,end(ranges) + 1L)
  starts <- starts[starts <= max_end]
  ends <- start(ranges) - 1L
  start_N <- ends > 0L
  starts <- starts[start_N]
  ends <- ends[start_N]
  if(length(ends) != length(starts)){
    stop("Something went wrong.")
  }
  if(length(starts) == 0L){
    names(seq) <- chromosome
    return(seq)
  }
  start_N <- start_N[1]
  starts <- starts[order(starts)]
  ends <- ends[order(ends)]
  f <- starts > ends
  starts <- starts[!f]
  ends <- ends[!f]
  N <- paste0(rep("N",max_end),collapse = "")
  N <- unlist(do.call(class(seq),list(unlist(N))))
  Ns <- as(Views(N,starts,ends),class(seq))
  common_length <- min(length(seq),length(Ns))
  common_seq <- seq_len(common_length)
  if(start_N){
    comb <- pc(Ns[common_seq],seq[common_seq])
  } else {
    missing_length_seq <- seq_along(seq)
    missing_length_seq <- missing_length_seq[missing_length_seq > common_length]
    comb <- pc(seq[common_seq],Ns[common_seq])
    comb <- Biostrings::xscat(comb,seq[missing_length_seq])
  }
  comb <- as(unlist(comb),class(seq))
  names(comb) <- chromosome
  comb
}

.get_viz_sequence_track <- function(seq, ranges, chromosome, args){
  FUN <- function(trackClass, seqClass, seq, args){
    set <- do.call(seqClass,list(seq))
    track <- do.call(trackClass,
                     c(list(sequence = set),
                       args[["sequence.track.pars"]]))
    track
  }
  if(!is(seq,"DNAStringSet") && 
     !is(seq,"RNAStringSet") && 
     !is(seq,"ModRNAStringSet")){
    stop("Invalid sequence type '",class(seq),"'. sequences(x) must be a ",
         "RNA/ModRNA/DNAStringSet.",
         call. = FALSE)
  }
  if(is(seq,"DNAStringSet")){
    seq <- as(seq,"RNAStringSet")
  }
  # reconstruct the chromosomal sequences for plotting
  seq <- .stitch_chromosome(seq, ranges, chromosome)
  if(is(seq,"RNAStringSet")){
    st <- FUN("RNASequenceTrack","RNAStringSet", seq, args)
  } else if(is(seq,"ModRNAStringSet")){
    st <- FUN("ModRNASequenceTrack","ModRNAStringSet", seq, args)
  } else {
    stop("Something went wrong.")
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
    data <- aggregateData(x)[name]
  } else if(is(x,"SequenceData")){
    x <- x[name]
    data <- aggregate(x)
  } else {
    stop("Something went wrong.")
  }
  strand_u <- .get_strand_u_GRangesList(ranges)
  seqs <- .seqs_rl(ranges)
  seqs[strand_u == "-"] <- rev(seqs[strand_u == "-"])
  seqnames <- .seqnames_rl(ranges)
  strands <- .strands_rl(ranges)
  ans <- GRanges(seqnames = unlist(seqnames),
                 ranges = IRanges(start = unlist(seqs),
                                  width = 1),
                 strand = unlist(strands),
                 data@unlistData)
  ans <- GRangesList(ans)
  ans@partitioning <- data@partitioning
  ans@metadata <- ranges@metadata
  ans
}

################################################################################

#' @rdname visualizeData
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "SequenceData", coord = "GRanges"),
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
  signature = signature(x = "SequenceData"),
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
    dt <- getDataTrack(x, name = name, ...)
    if(!is.list(dt)){
      dt <- list(dt)
    }
    tracks <- c(dt,list(st,atm))
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
  signature = signature(x = "SequenceData"),
  definition = function(x, name = name, ...) {
    stop("This functions needs to be implemented by '",class(x),"'.",
         call. = FALSE)
  }
)
