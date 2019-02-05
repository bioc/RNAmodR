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

.norm_viz_args_Modifier <- function(input){
  modified.seq <- FALSE
  additional.mod <- GRanges()
  colour <- NA
  sequence.track.pars <- list(add53 = FALSE)
  annotation.track.pars <- list()
  data.track.pars <- NULL
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
  if(!is.null(input[["colour"]])){
    colour <- input[["colour"]]
    if(!is.character(colour)){
      stop("'colour' must be a character vector parallel to the types selected",
           "or of length = 1, if now type selection is available.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["sequence.track.pars"]])){
    sequence.track.pars <- input[["sequence.track.pars"]]
  }
  if(!is.null(input[["annotation.track.pars"]])){
    annotation.track.pars <- input[["annotation.track.pars"]]
  }
  if(!is.null(input[["data.track.pars"]])){
    data.track.pars <- input[["data.track.pars"]]
  }
  args <- list(modified.seq = modified.seq,
               additional.mod = additional.mod,
               colour = colour,
               sequence.track.pars = sequence.track.pars,
               annotation.track.pars = annotation.track.pars,
               data.track.pars = data.track.pars)
  args
}

.norm_viz_windows.size <- function(window.size){
  if(!is.integer(window.size) || length(window.size) > 1L){
    stop("'window.size' must be a single integer value.",
         call. = FALSE)
  }
  window.size
}

.norm_viz_range <- function(range){
  GRanges(seqnames = .norm_viz_chromosome(range),
          ranges = ranges(range),
          strand = strand(range),
          mcols(range))
}

.norm_viz_chromosome <- function(range){
  "chr1"
}

.norm_viz_mod_annotation <- function(additional.mod,
                                     coord){
  FUN <- function(l){
    if(length(l) == 0L) return(NULL)
    l <- .norm_viz_range(l)
    strand(l) <- "*"
    l <- l[!is.na(l$mod) && !is.null(l$mod) && !(l$mod == "")]
    l$ID <- l$mod
    l$Activity <- NULL
    l
  }
  additional.mod <- .norm_addition.mod_for_visualization(additional.mod,
                                                         coord)
  additional.mod <- FUN(additional.mod)
  coord <- FUN(coord)
  list <- list(additional.mod,coord)
  list <- list[!vapply(list,is.null,logical(1))]
  ans <- unlist(GRangesList(list))
  ans[!duplicated(ans)]
}

.norm_coord_for_visualization <- function(coord){
  if(length(coord) > 1L){
    stop("'coord' must be contain a single range.",
         call. = FALSE)
  }
  if(is.null(coord$Parent)){
    stop("'coord' must contain a metadata column named 'Parent'.",
         call. = FALSE)
  }
  coord
}

.norm_addition.mod_for_visualization <- function(additional.mod,
                                                 coord){
  if(is(additional.mod,"GRangesList")){
    additional.mod <- unlist(additional.mod)
  }
  additional.mod <- additional.mod[additional.mod$Parent == coord$Parent]
  if(length(additional.mod) == 0L){
    return(additional.mod)
  }
  GRanges(seqnames = coord$Parent,
          ranges = ranges(additional.mod),
          strand = strand(additional.mod),
          mcols(additional.mod))
}

.norm_data_for_visualization <- function(data,type,chromosome){
  data <- data[,colnames(data) %in% type,drop=FALSE]
  if(ncol(data) == 0L){
    stop("No data found for the selected type(s): '",
         paste0(type,collapse = "','"),"'",
         call. = FALSE)
  }
  n <- nrow(data)
  ans <- GRanges(seqnames = rep(chromosome,n),
                 ranges = IRanges(start = seq_len(n),
                                  width = 1),
                 strand = "*",
                 data)
  ans
}

.norm_seqdata_for_visualization <- function(data,chromosome){
  n <- nrow(data)
  ans <- GRanges(seqnames = rep(chromosome,n),
                 ranges = IRanges(start = seq_len(n),
                                  width = 1),
                 strand = "*",
                 data)
  ans
}

.get_viz_window_Modifier <- function(data,coord,window.size){
  window.size <- .norm_viz_windows.size(window.size)
  start <- start(coord) - window.size
  end <- end(coord) + window.size
  pos.min <- as.integer(min(start(data)))
  pos.max <- as.integer(max(end(data)))
  if(start < pos.min){
    start <- pos.min
  }
  if(end > pos.max){
    end <- pos.max
  }
  list(start = start,
       end = end)
}

.get_viz_sequence <- function(seq,coord,args,modifications){
  if(args[["modified.seq"]]){
    if(length(modifications) > 0L){
      seq <- combineIntoModstrings(seq,modifications)
    }
  }
  seq
}

.get_viz_sequence_track <- function(seq,chromosome,args){
  FUN <- function(trackClass,seqClass,seq,chromosome,args){
    set <- do.call(seqClass,list(seq))
    names(set) <- chromosome
    track <- do.call(trackClass,list(set))
    if(!is.null(args[["sequence.track.pars"]])){
      Gviz::displayPars(track) <- args[["sequence.track.pars"]]
    }
    track
  }
  if(!is(seq,"DNAString") && !is(seq,"RNAString") && !is(seq,"ModRNAString")){
    stop("Invalid sequence type '",class(seq),"'. sequences(x) must be a ",
         "RNA/ModRNA/DNAString.",
         call. = FALSE)
  }
  if(is(seq,"DNAString")){
    seq <- as(seq,"RNAString")
  }
  if(is(seq,"RNAString")){
    st <- FUN("RNASequenceTrack","RNAStringSet",seq,chromosome,args)
  } else if(is(seq,"ModRNAString")){
    st <- FUN("ModRNASequenceTrack","ModRNAStringSet",seq,chromosome,args)
  } else {
    stop("Something went wrong.")
  }
  st
}

.get_viz_annotation_track <- function(modAnnotation,
                                      chromosome,
                                      args){
  at <- Gviz::AnnotationTrack(modAnnotation,
                              chromosome = chromosome,
                              collapse = TRUE,
                              name = "",
                              id = modAnnotation$ID,
                              showFeatureId = TRUE,
                              collapse = TRUE,
                              background.title = "#FFFFFF",
                              fontcolor.legend = "#000000")
  if(!is.null(args[["annotation.track.pars"]])){
    Gviz::displayPars(at) <- args[["annotation.track.pars"]]
  }
  at
}

#' @rdname visualizeData
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "Modifier", coord = "GRanges"),
  definition = function(x, coord, type = NA, window.size = 15L, ...) {
    requireNamespace("Gviz")
    # input check
    args <- .norm_viz_args_Modifier(list(...))
    coord <- .norm_coord_for_visualization(coord)
    # get plotting data
    modAnnotation <- .norm_viz_mod_annotation(args[["additional.mod"]],
                                              coord)
    seq <- .get_viz_sequence(sequences(x)[[coord$Parent]],
                             coord,
                             args,
                             modAnnotation)
    range <- .norm_viz_range(ranges(x)[[coord$Parent]])
    chromosome <- .norm_viz_chromosome(range)
    data <- .norm_data_for_visualization(aggregateData(x)[[coord$Parent]],
                                         type,
                                         chromosome)
    seqdata <- .norm_seqdata_for_visualization(
      aggregate(seqData(x)[coord$Parent])[[1]],
      chromosome)
    coordValues <- .get_viz_window_Modifier(data, coord, window.size)
    # get tracks
    st <- .get_viz_sequence_track(seq,chromosome,args)
    atm <- .get_viz_annotation_track(modAnnotation,chromosome,args)
    dt <- .dataTracks(x,
                      data = data,
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

.create_coord_for_visualization <- function(x,
                                            name,
                                            from,
                                            to){
  if(is(x,"Modifier")){
    if(!(name %in% names(seqData(x)))){
      stop("Element '",name,"' not present in data.", call. = FALSE)
    }
  } else if(is(x,"ModifierSet")){
    names <- lapply(lapply(x,seqData),names)
    if(any(!vapply(names,function(n){name %in% n},logical(1)))){
      stop("Element '",name,"' not present in data.", call. = FALSE)
    }
  } else {
    stop("Something went wrong.")
  }
  GRanges(seqnames = "chr1",
          ranges = IRanges::IRanges(start = as.integer(from),
                                    end = as.integer(to)),
          strand = "*",
          DataFrame(Parent = name))
}

#' @rdname visualizeData
#' @export
setMethod(
  f = "visualizeData",
  signature = signature(x = "Modifier"),
  definition = function(x, name, from, to, type = NA, ...) {
    coord <- .create_coord_for_visualization(x, name, from, to)
    visualizeDataByCoord(x, coord, type = type,  window.size = 0L, ...)
  }
)

#' @rdname RNAmodR-internals
setMethod(
  f = ".dataTracks",
  signature = signature(x = "Modifier",
                        data = "ANY",
                        seqdata = "ANY",
                        sequence = "ANY"),
  definition = function(x, data, seqdata, sequence, args) {
    stop(".dataTracks needs to be implemented for class '", class(x)[[1]], "'")
  }
)


# viz utils --------------------------------------------------------------------

#' @importFrom grDevices col2rgb
.is_colour <- function(x) {
  vapply(x,
         function(z) {
           tryCatch(is.matrix(grDevices::col2rgb(z)), 
                    error = function(e) FALSE)
         },
         logical(1))
}
