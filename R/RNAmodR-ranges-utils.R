#' @include class-RNAmodR-analysis-type.R
NULL

# common function handling positions -------------------------------------------

# converts the position on a transcript from start to end, to the correct
# genomic position
.convert_local_to_global_locations <- function(gff,
                                               loc){
  
  
  
  # browser()
  # Intron handling needs to be added
  strand <- unique(as.character(strand(gff)))
  if(strand == "-"){
    locations <- (end(gff) - loc) + 1
  } else {
    locations <- (start(gff) + loc) - 1
  }
  
  
  return(locations)
}

# converts global to local positions and modifies data accoringly
.convert_global_to_local_position <- function(gr,
                                              data){
  # interest in read's 5' position
  # reset to relative positions to gene start
  if(.is_on_minus_strand(gr)){
    positions <- BiocGenerics::end(data)
    positions <- BiocGenerics::end(gr) - positions + 1
  } else {
    positions <- BiocGenerics::start(data)
    positions <- positions - BiocGenerics::start(gr) + 1
  }
  return(positions)
}

# converts global to local positions and modifies data accoringly
.convert_global_to_local_position <- function(gr,
                                              data){
  # interest in read's 5' position
  # reset to relative positions to gene start
  if(.is_on_minus_strand(gr)){
    positions <- BiocGenerics::end(data)
    positions <- BiocGenerics::end(gr) - positions + 1
  } else {
    positions <- BiocGenerics::start(data)
    positions <- positions - BiocGenerics::start(gr) + 1
  }
  return(positions)
}

# offset positions based on how many positions to be removed the read has passed 
# from transcription start
.move_positions <- function(positions,
                            posToBeRemoved,
                            strand){
  x <- unlist(posToBeRemoved)
  positions <- positions[!(positions %in% x)]
  unlist(lapply(positions, function(position){
    if(.is_minus_strand(strand)){
      position <- position + length(x[x>position])
    } else {
      position <- position - length(x[x<position])
    }
    position
  }))
}

# offset positions based on how many positions to be removed the read has passed 
# from transcription start
.move_reads <- function(ga,
                        posToBeRemoved){
  browser()
  x <- unlist(posToBeRemoved)
  ga <- ga[!(BiocGenerics::start(ga) %in% x) &
             !(BiocGenerics::end(ga) %in% x)]
  unlist(lapply(ga, function(read){
    if(.is_minus_strand(read)){
      BiocGenerics::start(read) <- BiocGenerics::start(read) + length(x[x>BiocGenerics::start(read)])
      BiocGenerics::end(read) <- BiocGenerics::end(read) + length(x[x>BiocGenerics::end(read)])
    } else {
      BiocGenerics::start(read) <- BiocGenerics::start(read) - length(x[x>BiocGenerics::start(read)])
      BiocGenerics::end(read) <- BiocGenerics::end(read) - length(x[x>BiocGenerics::end(read)])
    }
    read
  }))
}

# offset positions based on how many positions the read has passed from
# transcription start. used the name for identification
.move_positions_named <- function(positions,
                                  posToBeRemoved,
                                  strand){
  x <- unlist(posToBeRemoved)
  names(positions) <- as.numeric(names(positions))
  positions <- positions[!(names(positions) %in% x)]
  names(positions) <- .move_positions(as.numeric(names(positions)),
                                      posToBeRemoved,
                                      strand)
  positions
}

# aggregate the positions occupied by introns
.get_intron_positions <- function(gff, 
                                  ID){
  # get childs
  childs <- .subset_gff_for_unique_transcript(gff, 
                                              ID, 
                                              wo.childs = FALSE)
  lapply(childs[grepl("intron",childs$type),], 
         function(intron){
           Biostrings::start(intron):Biostrings::end(intron)
         })
}

# get transcript sequence without removed
.get_transcript_sequence <- function(gff,ID,seq){
  # get gr
  gr <- .subset_gff_for_unique_transcript(gff, ID)
  # get a list of introns and the position which are void
  posToBeRemoved <- .get_intron_positions(gff,
                                          gr$ID)
  # convert DNAString into character vector
  seq <- .get_single_position_letters(seq)
  # set positions as names
  if(.is_on_minus_strand(gr)){
    names(seq) <- BiocGenerics::end(gr):BiocGenerics::start(gr)
  } else{
    names(seq) <- BiocGenerics::start(gr):BiocGenerics::end(gr)
  }
  # remove intron sequence
  seq <- seq[!(names(seq) %in% unlist(posToBeRemoved))]
  # set local position
  names(seq) <- 1:length(seq)
  return(seq)
}