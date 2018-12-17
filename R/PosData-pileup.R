#' @include RNAmodR.R
#' @include PosData-class.R
NULL

#' @name PileupPosData
#' 
#' @title PileupPosData
#' 
#' @description
#' title
NULL

#' @rdname PileupPosData
#' @export
setClass(Class = "PileupPosData",
         contains = "PosData",
         prototype = list(minQuality = 5L))

# PileupPosData ----------------------------------------------------------------

.pileup_colnames <- c("pos","-","G","A","T","C")
.pileup_measure_colnames <- c("-","G","A","T","C")

.fill_up_pileup_data <- function(d,r){
  pos <- seq_len(width(r))
  if(is.null(d)){
    colnames <- .pileup_colnames
    d <- data.frame(pos = pos,
                    "-" = 0,
                    "G" = 0,
                    "A" = 0,
                    "T" = 0,
                    "C" = 0)
  } else {
    missingPos <- pos[!(pos %in% d$pos)]
    missingCols <- .pileup_colnames[!(.pileup_colnames %in% colnames(d))]
    d <- reshape2::melt(d,
                        id.vars = "pos",
                        measure.vars = .pileup_measure_colnames[
                          .pileup_measure_colnames %in% colnames(d)])
    if(length(missingPos) > 0){
      d <- rbind(d,
                 data.frame(pos = missingPos,
                            variable = "-",
                            value = 0))
    }
    mpos <- max(d$pos)
    if(length(missingCols) > 0){
      d <- do.call(
        rbind,
        c(list(d),
          lapply(
            missingCols,
            function(c){
              data.frame(pos = mpos,
                         variable = c,
                         value = ifelse(is.null(d[d$pos == mpos,"value"]),
                                        0,
                                        d[d$pos == mpos,"value"]))
            })))
    }
    d <- reshape2::dcast(d, pos ~ variable, fun.aggregate = sum, fill = 0)
    colnames <- colnames(d)
  }
  df <- DataFrame(d)
  colnames(df) <- colnames
  df <- df[,.pileup_colnames]
  df
}

.get_position_data_of_transcript_pileup <- function(bamFile,
                                                    ranges,
                                                    sequences,
                                                    param,
                                                    args = list()){
  # get user set argumenst
  pileupArgs <- args[c("max_depth", 
                       "min_base_quality", 
                       "min_mapq", 
                       "min_nucleotide_depth", 
                       "min_minor_allele_depth", 
                       "distinguish_strands", 
                       "distinguish_nucleotides", 
                       "ignore_query_Ns", 
                       "include_deletions", 
                       "include_insertions", 
                       "left_bins", 
                       "query_bins", 
                       "cycle_bins")]
  parentRanges <- .get_parent_annotations(ranges)
  pileupArgs <- pileupArgs[!vapply(pileupArgs,is.null,logical(1))]
  # get data per chromosome
  pileupParam <- do.call("PileupParam",pileupArgs)
  pileup <- pileup(bamFile, scanBamParam = param, pileupParam = pileupParam)
  pileup <- DataFrame(pileup)
  # reformat data
  pileup$seqnames <- as.character(pileup$seqnames)
  pileup$nucleotide <- as.character(pileup$nucleotide)
  # split into data per transcript which is defined by the which_label column
  # format: chromosome:start-end
  pileup <- split(pileup,
                  pileup$which_label)
  # sanitize pilup data
  # - keep only data for correct strand
  # - fillup empty positions with zero
  strands <- as.character(strand(parentRanges))
  pileup <- SplitDataFrameList(
    mapply(
      function(d,r,strand){
        ans <- NULL
        d <- d[d$strand == strand,]
        if(nrow(d) > 0) {
          ans <- reshape2::dcast(
            as.data.frame(d),
            pos ~ nucleotide,
            sum,
            value.var = "count")
        }
        ans <- .fill_up_pileup_data(ans,r)
        # remove pos column since we don't need this anymore. seq_along == pos
        ans$pos <- NULL
        ans
      },
      pileup,
      split(parentRanges,seq_along(parentRanges)),
      strands,
      SIMPLIFY = FALSE))
  pileup
}

#' @name PileupPosData
#' @importFrom reshape2 dcast melt
#' 
#' @export
PileupPosData <- function(bamfiles,
                          fasta,
                          gff,
                          ...){
  ans <- new("PileupPosData",
             bamfiles,
             fasta,
             gff,
             ...)
  args <- .get_mod_data_args(...)
  ranges <- .load_annotation(ans@gff)
  sequences <- .load_transcript_sequences(ans@fasta,
                                          ranges)
  param <- .assemble_scanBamParam(ranges,
                                  ans@minQuality,
                                  ans@chromosomes)
  message("Loading Pileup data from BAM files...")
  data <- lapply(ans@bamfiles,
                 FUN = .get_position_data_of_transcript_pileup,
                 ranges = ranges,
                 sequences = sequences,
                 param = param,
                 args = args)
  names(data) <- paste0("pileup.",
                        names(ans@bamfiles),
                        ".",
                        seq_along(ans@bamfiles))
  .postprocess_read_data(ans,
                         data,
                         ranges,
                         sequences)
}
