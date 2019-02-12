#' @include RNAmodR.R
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

# GRanges/GRangesList helper functions -----------------------------------------

# per element entry
.get_strands_GRangesList <- function(grl){
  IRanges::CharacterList(strand(grl))
}
.get_column_GRangesList <- function(grl,column){
  columns <- IRanges::CharacterList(mcols(grl@unlistData)[,column])
  columns@partitioning <- grl@partitioning
  columns
}

# per element of GRangesList unique 
.get_strand_u_GRangesList <- function(grl){
  strand_u <- unique(IRanges::CharacterList(strand(grl)))
  ans <- unlist(strand_u)
  if(length(strand_u) != length(ans)){
    stop("Non unqiue strands per GRangesList element.")
  }
  ans
}

# per positions of each element

.seqnames_rl <- function(rl){
  seqnames <- as.character(seqnames(rl@unlistData))
  width <- as.integer(width(rl@unlistData))
  seqnames <- mapply(rep,seqnames,width,SIMPLIFY = FALSE)
  seqnames <- IRanges::CharacterList(lapply(mapply(seq.int,
                                                   start(rl@partitioning),
                                                   end(rl@partitioning),
                                                   SIMPLIFY = FALSE),
                                            function(i){
                                              unname(unlist(seqnames[i]))
                                            }))
  seqnames
}

.strands_rl <- function(rl){
  strands <- as.character(strand(rl@unlistData))
  width <- as.integer(width(rl@unlistData))
  strands <- mapply(rep,strands,width,SIMPLIFY = FALSE)
  strands <- IRanges::CharacterList(lapply(mapply(seq.int,
                                                  start(rl@partitioning),
                                                  end(rl@partitioning),
                                                  SIMPLIFY = FALSE),
                                           function(i){
                                             unname(unlist(strands[i]))
                                           }))
  strands
}


# Vectorize version of seq specific for start/ends from a RangesList
.seqs_rl <- function(rl){
  .seqs_rl_by(rl)
}

# Vectorize version of seq specific for start/ends from a RangesList with a by 
# option
.seqs_rl_by <- function(rl, by = 1L){
  starts <- unlist(start(rl))
  ends <- unlist(end(rl))
  .seqs_l_by(starts, ends, by)
}

#' @importFrom IRanges PartitioningByWidth PartitioningByEnd
#' @importClassesFrom IRanges PartitioningByWidth PartitioningByEnd
# Vectorize version of seq using to input lists
.seqs_l_by <- function(from, to, by = 1L){
  if(is.null(names(from)) || is.null(names(to))){
    stop("Inputs must be named.")
  }
  if(length(from) != length(to)){
    stop("Inputs must have the same length.")
  }
  if(by == 0L){
    stop("by cannot be zero.")
  }
  if(any(names(from) != names(to))){
    stop("Unmatched names.")
  }
  if(by < 0L){ # switch from to around if negative
    tmp <- to
    to <- from
    from <- tmp
    rm(tmp)
  }
  ans <- mapply(
    function(f,t){
      seq.int(f,t,by)
    },
    from,
    to,
    SIMPLIFY = FALSE)
  ans <- IRanges::IntegerList(ans)
  width_x <- IRanges::IntegerList(split(width(ans@partitioning),
                                        names(ans@partitioning)))
  m <- match(unique(names(from)),names(width_x))
  width_x <- width_x[m]
  width_ans <- sum(width_x)
  part <- IRanges::PartitioningByWidth(width_ans, names = names(width_ans))
  ans@partitioning <- as(part,"PartitioningByEnd")
  ans
}
