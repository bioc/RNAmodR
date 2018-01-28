#' @include class-RNAmod-analysis-type.R
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
.convert_global_to_local_position <- function(gff,
                                              gr,
                                              data){
  # discard reads out of boundaries
  data <- data[BiocGenerics::end(data) <= BiocGenerics::end(gr),]
  data <- data[BiocGenerics::start(data) >= BiocGenerics::start(gr),]
  # get a list of introns and the position which are void
  posToBeRemoved <- .get_intron_positions(gff,
                                          gr$ID)
  # get the genomic distance
  length <- width(gr) - length(unlist(posToBeRemoved))
  #
  #
  # GA specific handling
  # move position based on strand
  strand <- .get_unique_strand(gr)
  data <- data[.is_on_correct_strand(data,strand)]
  # interest in read's 5' position
  if(.is_on_minus_strand(gr)){
    positions <- BiocGenerics::end(data)
  } else {
    positions <- BiocGenerics::start(data)
  }
  # offset positions based on how many positions the read has passed from
  # transcription start
  positions <- unlist(lapply(positions, 
                             FUN = .move_postion,
                             posToBeRemoved, 
                             strand))
  # reset to relative positions to gene start
  if(.is_on_minus_strand(gr)){
    positions <- BiocGenerics::end(gr) - positions + 1
  } else {
    positions <- positions - BiocGenerics::start(gr) + 1
  }
  # discard reads out of boundaries
  # positions <- positions[positions>0]
  # positions <- positions[positions<=length]
  return(positions)
}

# offset positions based on how many positions the read has passed from
# transcription start
.move_postion <- function(positions,
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