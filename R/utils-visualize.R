#' @include RNAmodR.R
NULL

.norm_viz_coord <- function(data,
                            start,
                            end){
  if(missing(start)){
    start <- 1L
  }
  if(missing(end)){
    end <- start + 99L
  }
  pos <- seq_len(nrow(data))
  if(start < min(pos)){
    start <- as.integer(min(pos))
  }
  if(end > max(pos)){
    end <- as.integer(max(pos))
  }
  list(start = start,
       end = end)
}

.norm_viz_positions <- function(data,
                                coord){
  pos <- coord$start + seq(coord$start, coord$end) - 1L
  data[pos,]
}

.norm_viz_sequences<- function(seq,
                               coord){
  XVector::subseq(seq,coord$start,coord$end)
}

.norm_viz_chromosome <- function(range){
  chrom <- as.character(seqnames(range))
  "chrNA"
}

.norm_viz_genome <- function(range){
  as.character(genome(range))
}

.get_sequence_track <- function(range,
                                sequence){
  browser()
  SequenceTrack()
}