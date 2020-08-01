#' @include RNAmodR.R
NULL

# Annotation -------------------------------------------------------------------

#' @importFrom rtracklayer GFF3File
.norm_gff <- function(x, className, .xname = .get_name_in_parent(x)){
  if(!is(x,"GFF3File")){
    x <- try(rtracklayer::GFF3File(x), silent = TRUE)
    if (is(x, "try-error")){
      stop("To create a ",className," object, '",.xname,"' must be a GFF3File ",
           "object or be coercible to one.",
           call. = FALSE)
    }
  }
  if(!all(file.exists(BiocGenerics::path(x)))){
    stop("The gff3 file does not exist at the given location.",
         call. = FALSE)
  }
  x
}

#' @importFrom GenomicFeatures makeTxDbFromGFF
# Returns a TxDb or a GRangesList object
.norm_annotation <- function(annotation, className, 
                  .annotationname = .get_name_in_parent(annotation)){
  if(!is(annotation,"GRangesList")){
    if(!is(annotation,"GFFFile") && !is(annotation,"TxDb")){
      annotation <- .norm_gff(annotation, className, .annotationname)
    } else if(is(annotation,"GFFFile")) {
      if(!.all_are_existing_files(c(BiocGenerics::path(annotation)))){
        stop("annotation files don't exist or cannot be accessed.",
             call. = FALSE)
      }
    } else if(is(annotation,"TxDb")) {
      if(!all(validObject(annotation))){
        stop(".")
      }
    } else {
      stop("Something went wrong. Unrecognized annotation input during ",
           "creation of class '",className,"'.",
           call. = FALSE)
    }
    if(!is(annotation,"TxDb")){
      annotation <- GenomicFeatures::makeTxDbFromGFF(annotation)
    }
  } else {
    annotation <- .norm_annotation_GRangesList(annotation)
  }
  annotation
}

.norm_annotation_GRangesList <- function(annotation){
  if(is.null(names(annotation))){
    stop("Elements of 'annotation' GRangesList must be named.")
  }
  if(any(duplicated(names(annotation)))){
    stop("Names of elements in 'annotation' GRangesList must be unique.")
  }
  if("*" %in% unique(unlist(IRanges::CharacterList(strand(annotation))))){
    stop("Invalid strand information. Strand must either be '+' or '-'.")
  }
  annotation
}

# Sequences --------------------------------------------------------------------

#' @importFrom Rsamtools FaFile
# Either return a FaFile or BSgenome object
.norm_sequences <- function(seq, className){
  if(!is(seq,"FaFile") && !is(seq,"BSgenome")){
    tmp <- try(Rsamtools::FaFile(seq))
    if(is(tmp,"try-error")){
      stop("Input is not a FaFile and could not be coerced to one during ",
           "creation of class '",className,"'.",
           call. = FALSE)
    }
    seq <- tmp
    if(!.all_are_existing_files(c(BiocGenerics::path(seq)))){
      stop("sequence files don't exist or cannot be accessed.",
           call. = FALSE)
    }
  } else if(is(seq,"FaFile")) {
    if(!.all_are_existing_files(c(BiocGenerics::path(seq)))){
      stop("sequence files don't exist or cannot be accessed.",
           call. = FALSE)
    }
  } else if(is(seq,"BSgenome")) {
    if(!all(validObject(seq))){
      stop(".")
    }
  } else {
    stop("Something went wrong. Unrecognized sequence input during creation of",
         " class '",className,"'.",
         call. = FALSE)
  }
  seq
}

# BamFiles ---------------------------------------------------------------------

SAMPLE_TYPES <- c("treated","control")

#' @importFrom Rsamtools BamFileList
.norm_bamfiles <- function(x, className, 
                           .xname = .get_name_in_parent(x)){
  if(is.list(x)){
    if(!is.character(x[[1]]) && !is(x[[1]],"BamFile")){
      ans <- lapply(x, .norm_bamfiles, className, .xname)
      names(ans) <- names(x)
      return(ans)
    }
  }
  if(!is(x,"BamFileList")){
    tmp <- try(Rsamtools::BamFileList(x), silent = TRUE)
    if (is(tmp, "try-error")){
      stop("To create a ",className," object, '",.xname,"' must be a ",
           "BamFileList object or be coercible to one.",
           call. = FALSE)
    }
    x <- tmp
  }
  if(length(x) == 0L){
    stop("BamFileList is empty.",
         call. = FALSE)
  }
  if(!all(file.exists(BiocGenerics::path(x)))){
    stop("Bam files do not exists at the given locations.",
         call. = FALSE)
  }
  if(is.null(names(x))){
    stop("Names of BamFileList must either be 'treated' or 'control' (case ",
         "insensitive). No names found.")
  }
  x_names <- tolower(unique(names(x)))
  if(!all(x_names %in% SAMPLE_TYPES)){
    stop("Names of BamFileList must either be 'treated' or 'control' (case ",
         "insensitive).",
         call. = FALSE)
  }
  names(x) <- tolower(names(x))
  x
}

#' @importFrom Rsamtools BamFileList scanBamHeader
# retrieve a Seqinfo object from the bam headers
.bam_header_to_seqinfo <- function(bfl){
  if(is(bfl,"BamFile")){
    bfl <- Rsamtools::BamFileList(bfl)
  } else if(!is(bfl,"BamFileList")){
    message <- "BamFileList or list of BamFileList required."
    if(!is.list(bfl)){
      stop(message)
    }
    test <- !vapply(bfl,is,logical(1),"BamFileList")
    if(any(test)){
      stop(message)
    }
  }
  if(is.list(bfl)){
    headers <- lapply(bfl,
                      function(l){
                        lapply(l,Rsamtools::scanBamHeader)
                      })
    headers <- unlist(headers)
  } else {
    headers <- lapply(bfl,Rsamtools::scanBamHeader)
  }
  targets <- lapply(headers,"[[","targets")
  targets <- unique(do.call(c,lapply(unname(targets),names)))
  seqinfo <- GenomeInfoDb::Seqinfo(targets)
  seqinfo
}


# Misc -------------------------------------------------------------------------

# try to coerce the input to a Seqinfo object
.norm_seqinfo <- function(seqinfo){
  if(!is(seqinfo,"Seqinfo")){
    tmp <- try(GenomeInfoDb::Seqinfo(seqinfo), silent = TRUE)
    if(is(tmp,"try-error")){
      stop("Input is not a Seqinfo object and could not be coerced to ",
           "one.",
           call. = FALSE)
    }
    seqinfo <- tmp
  }
  seqinfo
}

# Retrieve the intersection of seqnames in annotation, sequence and seqinfo
# data
#' @importFrom BSgenome seqnames
#' @importFrom Rsamtools scanFa
.norm_seqnames <- function(bamfiles, annotation, sequences, seqinfo, className){
  if(missing(seqinfo)){
    seqinfo <- NULL
  }
  # norm seqinfo
  if(is.null(seqinfo) || 
     (!is(seqinfo,"Seqinfo") && (is.na(seqinfo)))){
    seqinfo <- .bam_header_to_seqinfo(bamfiles)
  }
  if(!is(seqinfo,"Seqinfo") && 
     (is(seqinfo,"BamFile") | is(seqinfo,"BamFileList"))){
    seqinfo <- .bam_header_to_seqinfo(seqinfo)
  }
  if(!is(seqinfo,"Seqinfo")){
    seqinfo <- .norm_seqinfo(seqinfo)
  }
  # norm annotation
  if(!is(annotation,"GRangesList") & !is(annotation,"TxDb")){
    stop("")
  }
  if(!is(sequences,"FaFile") & !is(sequences,"BSgenome")){
    stop("")
  }
  # norm sequences input
  if(is(sequences,"FaFile")){
    seq_seqnames <- names(Rsamtools::scanFa(sequences))
  } else {
    seq_seqnames <- BSgenome::seqnames(sequences)
  }
  seq_seqnames <- 
    seq_seqnames[seq_seqnames %in% GenomeInfoDb::seqlevels(annotation)]
  seq_seqnames <- 
    seq_seqnames[seq_seqnames %in% GenomeInfoDb::seqnames(seqinfo)]
  if( length(seqnames) == 0L ) {
    stop("No intersection between chromosome names in fasta, ",
         "annotation and seqinfo data.", 
         call. = FALSE)
  }
  return(seqinfo)
}

# Modifiertype -----------------------------------------------------------------

# check if a class of type x exists
.norm_modifiertype <- function(x){
  if(x == ""){
    stop("Empty string.")
  }
  class <- try(getClass(x), silent = TRUE)
  if (is(class, "try-error")){
    stop("Class '",x,"' is not implemented.",
         call. = FALSE)
  }
  if(isVirtualClass(class)){
    stop("Class '",x,"' is virtual.")
  }
  if(!("Modifier" %in% extends(class))){
    stop("Class '",x,"' does not extend the 'Modifier' class.")
  }
  nameId <- gregexpr("Mod",class@className)[[1]]
  if(any(nameId < 0L) || length(nameId) == 0L){
    stop("Invalid class name of Modifier class: the string 'Mod' must be ",
         "present once at the front of the class name.",
         call. = FALSE)
  }
  if(length(nameId) > 1L || nameId[1L] != 1L || 
     attr(nameId,"match.length")[1L] != 3L){
    stop("Invalid class name of Modifier class: the string 'Mod' can only be ",
         "present once at the front of the class name.",
         call. = FALSE)
  }
  x
}

.norm_mod <- function(x){
  modType <- modType(x)
  f <- .is_valid_modType(modType, seqtype(x))
  if(length(which(f)) != length(modType)){
    stop("Modification '",modType[!f],"' as defined for ",class(x)," does not ",
         "exist in the Modstrings dictionary for modified RNA/DNA sequences.",
         call. = FALSE)
  }
  modType
}

# check data validity ----------------------------------------------------------

.norm_sequence_data <- function(data){
  if(any(lengths(rownames(data)) == 0L)){
    stop("Sequence data does not contain rownames.")
  }
  data
}

.norm_aggregate_data <- function(data){
  if(any(lengths(rownames(data)) == 0L)){
    stop("Aggregated data does not contain rownames.")
  }
  data
}

# conditions and replicates ----------------------------------------------------

.norm_conditions <- function(x){
  conditions <- conditions(x)
  if(is(x,"SequenceDataList")){
    conditions <- conditions[[1L]]
  }
  conditions
}

.norm_replicates <- function(x){
  replicates <- replicates(x)
  if(is(x,"SequenceDataList")){
    replicates <- replicates[[1L]]
  }
  replicates
}
