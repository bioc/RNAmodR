


# sequence handling ------------------------------------------------------------

# get transcript sequence of multiple annotations
.get_seq_for_transcripts <- function(gff,fafile){
  seq <- getSeq(fafile,gff)
  seq[.is_on_minus_strand(gff)] <- rev(seq[.is_on_minus_strand(gff)])
  return(seq)
}

# get transcript sequence of one annotations
.get_seq_for_unique_transcript <- function(gr,fafile){
  if(length(gr) > 1){
    stop("Not a GRanges object of length = 1.",
         call. = FALSE)
  }
  return(unlist(.get_seq_for_transcripts(gr,fafile)))
}