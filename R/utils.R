#' @include RNAmodR.R
NULL

# seqnames related functions
.rebase_seqnames <- function(gr, seqnames){
  GenomicRanges::GRanges(seqnames = seqnames,
                         ranges = ranges(gr),
                         strand = strand(gr),
                         mcols(gr))
}

.get_which_label <- function(irl){
  which <- Map(
    function(ir,n){
      ans <- paste0(n,":",start(ir),"-",end(ir))
      names(ans) <- names(ir)
      ans
    },
    irl, names(irl))
  which <- IRanges::CharacterList(which)
  if(sum(vapply(which,anyDuplicated,integer(1))) != 0L){
    unlisted_which <- unlist(which, use.names = FALSE)
    names <- names(unlisted_which) 
    unlisted_which <- paste0(unlisted_which,".",
                             as.character(seq_along(unlisted_which)))
    names(unlisted_which) <- names
    which <- relist(unlisted_which,which)
  }
  which
}

# GRanges/GRangesList helper functions -----------------------------------------

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
  partitioning <- IRanges::PartitioningByEnd(rl)
  unlisted_rl <- unlist(rl)
  ul_seqnames <- as.character(seqnames(unlisted_rl))
  width <- as.integer(width(unlisted_rl))
  ul_seqnames <- Map(rep,ul_seqnames,width)
  ul_seqnames <- IRanges::CharacterList(lapply(Map(seq.int,
                                                   start(partitioning),
                                                   end(partitioning)),
                                               function(i){
                                                 unname(unlist(ul_seqnames[i]))
                                               }))
  ul_seqnames
}

.strands_rl <- function(rl){
  partitioning <- IRanges::PartitioningByEnd(rl)
  unlisted_rl <- unlist(rl)
  strands <- as.character(strand(unlisted_rl))
  width <- as.integer(width(unlisted_rl))
  strands <- Map(rep,strands,width)
  strands <- IRanges::CharacterList(lapply(Map(seq.int,
                                               start(partitioning),
                                               end(partitioning)),
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
  ans <- Map(
    function(f,t){
      seq.int(f,t,by)
    },
    from,
    to)
  ans <- IRanges::IntegerList(ans)
  partitioning <- IRanges::PartitioningByEnd(ans)
  width_x <- IRanges::IntegerList(split(width(partitioning),
                                        names(partitioning)))
  m <- match(unique(names(from)),names(width_x))
  width_x <- width_x[m]
  width_ans <- sum(width_x)
  ans <- relist(unname(unlist(ans)),
                IRanges::PartitioningByWidth(width_ans,
                                             names = names(width_ans)))
  unique(ans)
}

# DataFrame like helper functions ----------------------------------------------

# splits x along x$which_label. However, x$which_label is restructured to reflect 
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
  f_target <- vapply(split(width(IRanges::PartitioningByWidth(ans)), f_target),
                     sum, integer(1))
  f_target <- cumsum(f_target)
  f_target <- IRanges::PartitioningByEnd(f_target)
  relist(unlist(ans,use.names = FALSE),f_target)
}

# SequenceData helper functions ------------------------------------------------

# subset to conditions
.subset_to_condition <- function(conditions, condition){
  if(condition != "both"){
    f <- conditions == condition
    if(all(!f)){
      stop("No data for condition '",condition,"' found.")
    }
  } else {
    f <- rep(TRUE,length(conditions))
  }
  f
}

# partitioning object ----------------------------------------------------------

.seqs_partitioning <- function(partitioning){
  from <- rep.int(1,length(partitioning))
  to <- width(partitioning)
  names(from) <- names(partitioning)
  names(to) <- names(partitioning)
  .seqs_l_by(from,to)
}
