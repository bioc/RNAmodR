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

.get_ModExperiment_objects <- function(...){
  args <- list(...)
  ModExperiments <- args[lapply(args,is,"ModExperiment")]
  ModExperiments
}
