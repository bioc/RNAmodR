#' @include RNAmodR.R
NULL

# converts a single GRanges object into a GPos object
.create_per_position_grange <- function(gr,
                                        seq){
  gpos <- GPos(gr)
  attr(gpos, "ID") <- mcols(gr)$ID
  attr(gpos, "Name") <- mcols(gr)$Name
  attr(gpos, "Alias") <- mcols(gr)$Alias
  attr(gpos, "gene") <- mcols(gr)$gene
  mcols(gpos)$nucleotide <- split(seq, 1:length(seq))
  mcols(gpos)$counts <- 0
  return(gpos)
}

# remove data, which lies in intro sequences
.keep_non_intron_positions <- function(gpos,
                                       gff,
                                       id){
  # get a list of introns and the position which are void
  posToBeRemoved <- unlist(.get_intron_positions(gff,
                                                 id))
  if(length(posToBeRemoved) > 0){
    gpos <- gpos[!(pos(gpos) %in% posToBeRemoved)]
  }
  return(gpos)
}
