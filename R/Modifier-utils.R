#' @include RNAmodR.R
NULL

.norm_bamfiles <- function(x,
                           className,
                           .xname = assertive::get_name_in_parent(x)){
  # check bam files
  if(!is(x,"BamFileList")){
    x <- try(BamFileList(x))
    if (is(x, "try-error")){
      stop("To create a ",className," object, '",.xname,"' must be a ",
           "BamFileList object or be coercible to one.",
           call. = FALSE)
    }
  }
  if(!all(file.exists(path(x)))){
    stop("Some Bam files do not exists at the given locations.",
         call. = FALSE)
  }
  names <- tolower(unique(names(x)))
  if(!all(names %in% SAMPLE_TYPES)){
    stop("Names of BamFileList must either be 'Treated' or 'Control'",
         call. = FALSE)
  }
  names(x) <- tolower(names(x))
  x
}

.norm_fasta <- function(x,
                        className,
                        .xname = assertive::get_name_in_parent(x)){
  if(!is(x,"FaFile")){
    x <- try(FaFile(x))
    if (is(x, "try-error")){
      stop("To create a ",className," object, '",.xname,"' must be a FaFile ",
           "object or be coercible to one.",
           call. = FALSE)
    }
  }
  if(!all(file.exists(path(x)))){
    stop("The fasta file does not exist at the given location.",
         call. = FALSE)
  }
  indexFa(x)
  x
}

.norm_gff <- function(x,
                      className,
                      .xname = assertive::get_name_in_parent(x)){
  if(!is(x,"GFF3File")){
    x <- try(GFF3File(x))
    if (is(x, "try-error")){
      stop("To create a ",className," object, '",.xname,"' must be a GFF3File ",
           "object or be coercible to one.",
           call. = FALSE)
    }
  }
  if(!all(file.exists(path(x)))){
    stop("The gff3 file does not exist at the given location.",
         call. = FALSE)
  }
  x
}

.norm_mod <- function(mod,
                      className){
  f <- which(mod == shortName(ModRNAString()))
  if(length(f) != 1){
    stop("Modification '",mod,"' as defined for ",className," does not exist ",
         "in the Modstrings dictionary for modified RNA sequences.",
         call. = FALSE)
  }
  mod
}

.norm_data_type <- function(ans,pd){
  if(!is(pd,ans@dataClass)){
    stop("Data class '",ans@dataClass,"' is required by '",class(.Object),"'.",
         "\n'",class(pd),"' was provided. Aborting...",
         call. = FALSE)
  }
  pd
}

.norm_modifications <- function(ans,modifications){
  # ToDo: implement sanity check for given modifications
  modifications
}


.get_ModExperiment_objects <- function(...){
  args <- list(...)
  ModExperiments <- args[lapply(args,is,"ModExperiment")]
  ModExperiments
}

# data aggregation functions ---------------------------------------------------

#' @importFrom matrixStats rowSds
.aggregate_pile_up <- function(data){
  if(is(data,"CompressedSplitDataFrameList")){
    df <- data@unlistData
  } else {
    df <- data
  }
  replicates <- unique(data@replicate)
  for(i in seq_along(replicates)){
    df[,data@replicate == i] <- 
      as.data.frame(df[,data@replicate == i]) / 
      rowSums(as.data.frame(df[,data@replicate == i]))
  }
  ncol <- ncol(df[,data@replicate == 1L])
  seqAdd <- seq.int(from = 0, to = ncol(df) - 1, by = ncol)
  colNames <- as.character(
    data.frame(strsplit(colnames(df)[seq_len(ncol)],"\\."),
               stringsAsFactors = FALSE)[4,])
  means <- NumericList(lapply(seq_len(ncol),
                              function(i){
                                rowMeans(as.data.frame(df[,i + seqAdd]),
                                         na.rm = TRUE)
                              }))
  names(means) <- paste0("means.",colNames)
  sds <- NumericList(lapply(seq_len(ncol),
                            function(i){
                              matrixStats::rowSds(as.matrix(df[,i + seqAdd]),
                                                  na.rm = TRUE)
                            }))
  names(sds) <- paste0("sds.",colNames)
  ans <- cbind(do.call(DataFrame, means),
               do.call(DataFrame, sds))
  if(is(data,"CompressedSplitDataFrameList")){
    ans <- SplitDataFrameList(ans)
    ans@partitioning <- data@partitioning
  }
  ans
}

#' @importFrom matrixStats rowSds
.aggregate_pile_up_to_coverage <- function(data){
  if(is(data,"CompressedSplitDataFrameList")){
    df <- data@unlistData
  } else {
    df <- data
  }
  replicates <- unique(data@replicate)
  ans  <- IntegerList(lapply(seq_along(replicates),
                             function(i){
                               rowSums(as.data.frame(df[,data@replicate == i]))
                             }))
  names(ans) <- paste0("replicate.",replicates)
  ans <- do.call("DataFrame",ans)
  if(is(data,"CompressedSplitDataFrameList")){
    ans <- SplitDataFrameList(ans)
    ans@partitioning <- data@partitioning
  }
  ans
}

.construct_mod_ranges <- function(range,
                                  data,
                                  modType){
  positions <- as.integer(rownames(data))
  if(as.character(strand(range)) == "-"){
    positions <- end(range) - positions + 1L
  } else {
    positions <- start(range) + positions - 1L
  }
  mranges <- GRanges(seqnames = rep(as.character(seqnames(range)),
                                    nrow(data)),
                     ranges = IRanges::IRanges(start = positions,
                                               width = 1L),
                     strand = strand(range),
                     modType = rep(modType,nrow(data)),
                     score = .get_inosine_score(data))
  mranges
}

.get_inosine_score <- function(data){
  data$means.G / data$means.A
}