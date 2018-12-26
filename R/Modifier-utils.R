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
  f <- mod %in% shortName(ModRNAString())
  if(length(which(f)) != length(mod)){
    stop("Modification '",mod[!f],"' as defined for ",className," does not ",
         "exist in the Modstrings dictionary for modified RNA sequences.",
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

.norm_modifications <- function(ans,args){
  # ToDo: implement sanity check for given modifications
  ans
}


.get_ModExperiment_objects <- function(...){
  args <- list(...)
  ModExperiments <- args[lapply(args,is,"ModExperiment")]
  ModExperiments
}


# construct GRanges object from found modifications ----------------------------

.construct_mod_ranges <- function(range,
                                  data,
                                  modType,
                                  scoreFun,
                                  source,
                                  type){
  positions <- as.integer(rownames(data))
  if(as.character(strand(range)) == "-"){
    positions <- end(range) - positions + 1L
  } else {
    positions <- start(range) + positions - 1L
  }
  mranges <- do.call("GRanges",
                     c(list(seqnames = rep(as.character(seqnames(range)),
                                           nrow(data)),
                            ranges = IRanges::IRanges(start = positions,
                                                      width = 1L),
                            strand = strand(range),
                            seqinfo = seqinfo(range),
                            mod = rep(modType,nrow(data))),
                       source = source,
                       type = type,
                       do.call(scoreFun,list(data)),
                       list(Parent = range$ID)))
  mranges
}

# error propagation ------------------------------------------------------------

################################################################################
# modified from Lee Pang, 2015,
# (https://oddhypothesis.blogspot.com/2015/01/easy-error-propagation-in-r.html)
# (https://www.r-bloggers.com/easy-error-propagation-in-r/)
################################################################################

#' @importFrom dplyr mutate_
#' @importFrom stats D
.mutate_with_error <- function(.data, 
                               f){
  exprs = list(
    # expression to compute new variable values
    deparse(f[[3]]),
    # expression to compute new variable errors
    sprintf('sqrt(%s)',
            paste(sapply(all.vars(f[[3]]),
                         function(v) {
                           dfdp = deparse(stats::D(f[[3]], 
                                                   v))
                           sprintf('(d%s*(%s))^2', 
                                   v,
                                   dfdp)
                         }),
                  collapse = '+'))
  )
  names(exprs) = c(
    deparse(f[[2]]),
    sprintf('d%s', deparse(f[[2]]))
  )
  if( nrow(.data) > 0 ){
    # the standard evaluation alternative of mutate()
    .data <- dplyr::mutate_(.data, .dots = exprs)
  }
  return(.data)
}