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

# seqnames related functions

.get_seqnames <- function(x){
  as.character(GenomicRanges::seqnames(x))
}
.get_unique_seqnames <- function(x){
  unique(.get_seqnames(x))
}

.rebase_seqnames <- function(gr, seqnames){
  GenomicRanges::GRanges(seqnames = seqnames,
                         ranges = ranges(gr),
                         strand = strand(gr),
                         mcols(gr))
}

# seqlengths related functions
.rebase_GRanges <- function(gr){
  usn <- .get_unique_seqnames(gr)
  seqnames <- Rle(factor(GenomicRanges::seqnames(gr), levels = usn))
  seqlengths <- GenomeInfoDb::seqlengths(gr)[usn]
  seqinfo <- GenomeInfoDb::Seqinfo(usn, seqlengths)
  GenomicRanges::GRanges(seqnames = seqnames,
                         ranges = IRanges::ranges(gr),
                         strand = BiocGenerics::strand(gr),
                         seqinfo = seqinfo,
                         mcols(gr))
}

# GRanges/GRangesList helper functions -----------------------------------------

# per element entry
.get_strands_GRangesList <- function(grl){
  IRanges::CharacterList(strand(grl))
}
.get_column_GRangesList <- function(grl,column){
  relist(mcols(grl@unlistData)[,column],grl@partitioning)
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
.seqs_rl_strand <- function(rl, force_continous = FALSE, 
                            minus_decreasing = FALSE){
  strand_u <- .get_strand_u_GRangesList(rl)
  strand_minus <- strand_u == "-"
  ansP <- .seqs_rl_by(rl[!strand_minus])
  ansM <- .seqs_rl_by(rl[strand_minus], by = -1L)
  if(force_continous){
    ansM <- ansM[IRanges::IntegerList(lapply(ansM, order,
                                             decreasing = minus_decreasing))]
  }
  ans <- c(ansP, ansM)
  ans[match(names(rl),names(ans))]
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

# DataFrame like helper functions ----------------------------------------------

# splits x along x$which_labe. However, x$which_labe is restructured to reflect 
# length GRanges elements in a GRangesList. This is helpful to split data along
# transcripts instead of exons
#' @importFrom IRanges splitAsList
.splitPileupAsList_transcript <- function(x, grl, drop = FALSE){
  ans <- IRanges::splitAsList(x, x$which_label, drop)
  names(ans) <- vapply(strsplit(names(ans),"\\."),"[[",character(1),1)
  ugrl <- unlist(grl)
  f_order <- paste0(seqnames(ugrl),":",start(ugrl),"-",end(ugrl))
  f_order_match <- match(f_order,names(ans))
  if(anyNA(f_order_match)){
    f_order_match <- f_order_match[!is.na(f_order_match)]
  }
  ans <- ans[f_order_match]
  f_target <- unlist(mapply(rep, names(grl), lengths(grl)))
  f_target <- f_target[!is.na(f_order_match)]
  f_target <- factor(unname(f_target), levels = unique(f_target))
  f_target <- vapply(split(width(PartitioningByWidth(ans)), f_target),
                     sum, integer(1))
  f_target <- cumsum(f_target)
  f_target <- PartitioningByEnd(f_target)
  relist(unlist(ans),f_target)
}

# SequenceData helper functions ------------------------------------------------

# subset to conditions
.subset_to_condition <- function(conditions, condition){
  if(condition != "both"){
    f <- conditions == condition
    if(all(f == FALSE)){
      stop("No data for condition '",condition,"' found.")
    }
  } else {
    f <- rep(TRUE,length(conditions))
  }
  f
}

# Modstrings related helper functions ------------------------------------------

.is_valid_modType <- function(modType){
  modType %in% Modstrings::shortName(Modstrings::ModRNAString())
}
