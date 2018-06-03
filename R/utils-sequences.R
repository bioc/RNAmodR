


# sequence handling ------------------------------------------------------------

# retrieve concats DNAStrings based on the strand
.get_seq_for_unique_transcript <- function(gff,fafile){
  if(length(unique(as.character(BiocGenerics::strand(gff)))) != 1) {
    stop("Ambigeous type of GRanges given. Expects only one strand.",
         call. = FALSE)
  }
  seq <- getSeq(fafile,gff)
  if(as.character(BiocGenerics::strand(gff)) == "-"){
    return(unlist(rev(seq)))
  }
  return(unlist(seq))
}