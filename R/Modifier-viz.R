#' @include RNAmodR.R
#' @include Modifier-class.R
NULL


#' @name visualizeDataByCoord
#' 
#' @title visualizeDataByCoord
#' 
#' @description 
#' title
#' 
#' @param x a \code{Modifier} or \code{ModifierSet} object.
#' @param coord coordinates of a single position to subset to as a 
#' \code{GRanges} object. The Parent column is expected to match the gene or 
#' transcript name.
#' @param type the data type of data show as data tracks.
#' @param window.size integer value for the number of positions on the left and 
#' right site of the selected positions included in the plotting (default: 
#' \code{window.size = 15L})
#' @param ... optional parameters:
#' \itemize{
#' \item{\code{modified.seq}}{\code{TRUE} or \code{FALSE}. Should the sequence 
#' shown with modified positions modified? (default: 
#' \code{modified.seq = FALSE})}
#' }
NULL

.norm_viz_args <- function(input){
  modified.seq <- FALSE
  additional.mod <- GRanges()
  colour <- NA
  if(!is.null(input[["modified.seq"]])){
    modified.seq <- input[["modified.seq"]]
    if(!assertive::is_a_bool(modified.seq)){
      stop("'modified.seq' must be a single logical value.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["additional.mod"]])){
    additional.mod <- input[["additional.mod"]]
    if(!is(additional.mod,"GRanges")){
      stop("'additional.mod' must be a GRanges object, which is compatible ",
           "with combineIntoModstrings().",
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
  args <- list(modified.seq = modified.seq,
               additional.mod = additional.mod,
               colour = colour)
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

.norm_viz_mod_annotation <- function(list){
  list <- lapply(list,
                 function(l){
                   if(length(l) == 0L) return(NULL)
                   l <- .norm_viz_range(l)
                   strand(l) <- "*"
                   l$ID <- l$mod
                   l$Activity <- NULL
                   l
                 })
  list <- list[!vapply(list,is.null,logical(1))]
  ans <- unlist(GRangesList(list))
  ans[!duplicated(ans)]
}

.norm_GRanges_for_visualization <- function(coord){
  if(length(coord) > 1L){
    stop("'coord' must be contain a single range.",
         call. = FALSE)
  }
  if(is.null(coord$Parent)){
    stop("'coord' must be a metadata column named 'Parent'.",
         call. = FALSE)
  }
  coord
}

.norm_data_for_visualization <- function(data,type,chromosome){
  data <- data[,colnames(data) %in% type]
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

.get_viz_window <- function(data,coord,window.size){
  window.size <- .norm_viz_windows.size(window.size)
  start <- start(coord) - window.size
  end <- start(coord) + window.size
  pos <- start(data)
  if(start < min(pos)){
    start <- as.integer(min(pos))
  }
  if(end > max(pos)){
    end <- as.integer(max(pos))
  }
  list(start = start,
       end = end)
}

.get_viz_sequence <- function(seq,coord,args){
  if(args[["modified.seq"]]){
    modifications <- args[["additional.mod"]]
    modifications <- modifications[modifications$Parent == coord$Parent]
    if(length(modifications) > 0L){
      seq <- combineIntoModstrings(seq,modifications)
    }
  }
  seq
}

.get_viz_sequence_track <- function(seq,chromosome){
  if(is(seq,"RNAString")){
    set <- RNAStringSet(c(chromosome = seq))
    names(set) <- chromosome
    st <- RNASequenceTrack(set)
  } else if(is(seq,"ModRNAString")){
    set <- ModRNAStringSet(c(chromosome = seq))
    names(set) <- chromosome
    st <- ModRNASequenceTrack(set)
  } else {
    stop("Something went wrong.")
  }
  st
}

#' @rdname visualizeDataByCoord
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "Modifier",
                        coord = "GRanges"),
  definition = function(x,
                        coord,
                        type = NA,
                        window.size = 15L,
                        ...) {
    requireNamespace("Gviz")
    browser()
    # input check
    args <- .norm_viz_args(list(...))
    coord <- .norm_GRanges_for_visualization(coord)
    # get plotting data
    seq <- .get_viz_sequence(sequences(x)[[coord$Parent]],
                             coord,
                             args)
    range <- .norm_viz_range(ranges(x)[[coord$Parent]])
    chromosome <- .norm_viz_chromosome(range)
    data <- .norm_data_for_visualization(aggregateData(x)[[coord$Parent]],
                                         type,
                                         chromosome)
    coordValues <- .get_viz_window(data,coord,window.size)
    modAnnotation <- .norm_viz_mod_annotation(list(args[["additional.mod"]],
                                                   coord))
    # get tracks
    st <- .get_viz_sequence_track(seq,chromosome)
    atm <- AnnotationTrack(modAnnotation,
                           chromosome = chromosome,
                           collapse = TRUE,
                           name = "",
                           id = modAnnotation$ID,
                           showFeatureId = TRUE,
                           collapse = TRUE)
    dt <- .dataTracksByCoord(x,data,args)
    # plot tracks
    plotTracks(c(dt,
                 list(st,atm)),
               from = coordValues$start,
               to = coordValues$end,
               chromosome = chromosome)
  }
)

setMethod(
  f = ".dataTracksByCoord",
  signature = signature(x = "Modifier",
                        data = "ANY"),
  definition = function(x,
                        data,
                        args) {
    stop(".dataTracksByCoord needs to be implemented for class '",class(x)[[1]],
         "'")
  }
)
