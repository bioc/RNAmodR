#' @include class-RNAmodR-quant.R
NULL

# BiocGeneric helper functions -------------------------------------------------

# strand related functions
.get_strand <- function(x){
  as.character(BiocGenerics::strand(x))
}
.get_unique_strand <- function(x){
  unique(.get_strand(x))
}
.is_minus_strand <- function(x) {
  all(as.logical(.get_strand(x) == "-"))
}
.is_on_minus_strand <- function(x) {
  all(.is_on_correct_strand(x,"-"))
}
# is on the minus strand?
.is_on_correct_strand <- function(x, strand) {
  as.logical(.get_strand(x) == strand) 
}
.is_on_correct_strand2 <- function(x, gr) {
  as.logical(.get_strand(x) == .get_unique_strand(gr)) 
}

# GRanges helper functions -----------------------------------------------------

# get aggregated identifiers
.get_gr_ids <- function(gr,
                        na.rm = TRUE){
  IDs <- unique(c(gr$ID,gr$Name))
  if(na.rm){
    IDs <- IDs[!is.na(IDs)]
  }
  IDs
}

# get unique ids
# returns all identifiers which can be used as a identifier in the Parent
# column aggrgated by precident: ID > Name > gene
.get_unique_identifiers <- function(gr){
  ids <- as.character(gr$ID)
  ids[is.na(ids)] <- as.character(gr[is.na(ids)]$Name)
  ids[is.na(ids)] <- as.character(gr[is.na(ids)]$gene)
  ids
}

# subset to types present in RNAMODR_MOD_SEQ_FEATURES
.subset_rnamod_transcript_features <- function(gr){
  gr[S4Vectors::mcols(gr)$type %in% RNAMODR_MOD_TRANSCRIPT_FEATURES,]
}
# subset to types present in RNAMODR_MOD_SEQ_FEATURES
.subset_rnamod_containing_features <- function(gr){
  gr[S4Vectors::mcols(gr)$type %in% RNAMODR_MOD_CONTAINING_FEATURES,]
}

# subsets a gRanges object for entries concerning a single ID
# - Name/ID for parent entry
# - all childs based on Parent name equals ID
.subset_gff_for_unique_transcript <- function(gff, 
                                              ID,
                                              wo.childs = TRUE){
  gr <- gff[(is.na(S4Vectors::mcols(gff)$ID) & 
               S4Vectors::mcols(gff)$Name == ID) |
              (!is.na(S4Vectors::mcols(gff)$ID) & 
                 S4Vectors::mcols(gff)$ID == ID),]
  if(!wo.childs){
    gr <- append(gr, .get_childs(gff,
                                 as.character(GenomeInfoDb::seqnames(gr)),
                                 ID))
  }
  .order_GRanges(gr)
}
# recursive search for all childs of an ID
.get_childs <- function(gff,
                        chrom,
                        ID){
  gr <- gff[as.character(GenomeInfoDb::seqnames(gff)) == chrom & 
              !is.na(as.character(gff$Parent)) & 
              as.character(gff$Parent) == ID,]
  if(length(gr) == 0) return(GRanges())
  # Use ID and Name for child search
  IDs <- .get_gr_ids(gr)
  IDs <- IDs[IDs != ID]
  # Walk up the parent chain
  grl <- append(list(gr),lapply(IDs, function(x){.get_childs(gff,chrom,x)}))
  grl <- grl[!vapply(grl,is.null,logical(1))]
  return(unlist(GRangesList(grl)))
}

# subset GRanges for highest ranking parent
.get_parent_annotations <- function(gr,
                                    forceSingle = FALSE,
                                    doRecursiveSearch = FALSE,
                                    IDs){
  if(!doRecursiveSearch){
    res <- gr[is.na(as.character(gr$Parent)),]
    if(length(res) > 1 && forceSingle) 
      stop("No distinct single parent annotation detectable.")
    return(res)
  }
  if(missing(IDs)) IDs <- .get_gr_ids(gr)
  l <- lapply(IDs, .get_parent, gr = gr)
  l <- unique(unlist(l[!vapply(l,is.null,logical(1))]))
  .order_GRanges(gr[(!is.na(gr$ID) & gr$ID %in% l) |
                      (!is.na(gr$Name) & gr$Name %in% l),])
}
# quasi recursive search for all possible parents
.get_parent <- function(ID,
                        gr){
  parent <- as.character(gr[(!is.na(gr$ID) & gr$ID == ID) |
                              (!is.na(gr$Name) & gr$Name == ID),]$Parent)
  parent <- parent[!is.na(parent)]
  if(length(parent) == 0) return(ID)
  return(unlist(lapply(parent, .get_parent, gr)))
}

# returns an order based on the named of the chromosome and start and end 
# coordinates
.order_GRanges <- function(gr){
  gr[order(GenomeInfoDb::seqnames(gr),
           BiocGenerics::start(gr),
           BiocGenerics::end(gr))]
}

# common function handling positions -------------------------------------------

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

# offset positions based on how many positions the read has passed from
# transcription start. used the name for identification
.move_positions_named <- function(positions,
                                  posToBeRemoved,
                                  strand){
  x <- unlist(posToBeRemoved)
  positions <- positions[!(as.numeric(names(positions)) %in% x)]
  names(positions) <- .move_positions(as.numeric(names(positions)),
                                      posToBeRemoved,
                                      strand)
  return(positions)
}

# sequences --------------------------------------------------------------------

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

.get_single_position_letters <- function(x) {
  x <- as.character(x)
  substring(x, 1:nchar(x), 1:nchar(x))  
}

# Global to Local --------------------------------------------------------------

# converts global to local positions and modifies data accoringly
.convert_global_coverage_to_local_coverage <- function(gff,
                                                       gr,
                                                       data,
                                                       posToBeRemoved){
  # get coverage
  coverage <- as.numeric(GenomicAlignments::coverage(data, 
                                                     shift = -BiocGenerics::start(gr)+1,
                                                     width = BiocGenerics::width(gr),
                                                     method = "hash")[GenomeInfoDb::seqnames(gr)][[1]])
  names(coverage) <- BiocGenerics::start(gr):BiocGenerics::end(gr)
  # take care of intron positions
  coverage <- .move_positions_named(coverage, 
                                    posToBeRemoved, 
                                    .get_unique_strand(gr))
  # reset to relative positions to gene start
  if(.is_on_minus_strand(gr)){
    names(coverage) <- BiocGenerics::end(gr) - as.numeric(names(coverage)) + 1
  } else {
    names(coverage) <- as.numeric(names(coverage)) - BiocGenerics::start(gr) + 1
  }
  coverage <- coverage[order(as.numeric(names(coverage)))]
  return(coverage)
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

# Local to Global --------------------------------------------------------------

# converts the position on a transcript from start to end, to the correct
# genomic position
.convert_local_to_global_locations <- function(gff,
                                               loc){
  # Intron handling needs to be added
  strand <- unique(as.character(strand(gff)))
  if(strand == "-"){
    locations <- (end(gff) - loc) + 1
  } else {
    locations <- (start(gff) + loc) - 1
  }
  return(locations)
}
