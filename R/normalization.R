#' @include RNAmodR.R
NULL

# Annotation -------------------------------------------------------------------

#' @importFrom rtracklayer GFF3File 
#' @importFrom BiocGenerics path
.norm_gff <- function(x, className, .xname = assertive::get_name_in_parent(x)){
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
#' @importFrom BiocGenerics path
# Return a TxDb object
.norm_annotation <- function(annotation, className, 
                  .annotationname = assertive::get_name_in_parent(annotation)){
  if(!is(annotation,"GFFFile") && !is(annotation,"TxDb")){
    annotation <- .norm_gff(annotation, className, .annotationname)
  } else if(is(annotation,"GFFFile")) {
    assertive::assert_all_are_existing_files(c(BiocGenerics::path(annotation)))
  } else if(is(annotation,"TxDb")) {
    assertive::assert_all_are_true(validObject(annotation))
  } else {
    stop("Something went wrong. Unrecognized annotation input during ",
         "creation of class '",className,"'.",
         call. = FALSE)
  }
  if(!is(annotation,"TxDb")){
    annotation <- GenomicFeatures::makeTxDbFromGFF(annotation)
  }
  annotation
}

# Sequences --------------------------------------------------------------------

#' @importFrom Rsamtools FaFile
#' @importFrom BiocGenerics path
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
    assertive::assert_all_are_existing_files(c(BiocGenerics::path(seq)))
  } else if(is(seq,"FaFile")) {
    assertive::assert_all_are_existing_files(c(BiocGenerics::path(seq)))
  } else if(is(seq,"BSgenome")) {
    assertive::assert_all_are_true(validObject(seq))
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
#' @importFrom BiocGenerics path
.norm_bamfiles <- function(x, className, 
                           .xname = assertive::get_name_in_parent(x)){
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
    stop("Some Bam files do not exists at the given locations.",
         call. = FALSE)
  }
  names <- tolower(unique(names(x)))
  if(!all(names %in% SAMPLE_TYPES)){
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
  }
  if(!is(bfl,"BamFileList")){
    stop("BamFileList required.")
  }
  headers <- lapply(bfl,Rsamtools::scanBamHeader)
  targets <- lapply(headers,"[[","targets")
  targets <- unique(do.call(c,lapply(unname(targets),names)))
  seqinfo <- GenomeInfoDb::Seqinfo(targets)
  seqinfo
}


# Misc -------------------------------------------------------------------------

# try to coerce the input to a Seqinfo object
.norm_seqinfo <- function(seqinfo){
  if(!is(seqinfo,"Seqinfo")){
    tmp <- try(GenomeInfoDb::Seqinfo(seqinfo))
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
  # norm seqinfo
  if(missing(seqinfo) || 
     (!is(seqinfo,"Seqinfo") && (is.na(seqinfo) || is.null(seqinfo)))){
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
  if(!is(annotation,"TxDb")){
    annotation  <- .norm_annotation(annotation, className)
  }
  # norm sequences input
  if(is(sequences,"FaFile")){
    seqnames <- names(Rsamtools::scanFa(sequences))
  } else if(is(sequences,"BSgenome")) {
    seqnames <- BSgenome::seqnames(sequences)
  }
  seqnames <- seqnames[seqnames %in% GenomeInfoDb::seqlevels(annotation)]
  seqnames <- seqnames[seqnames %in% GenomeInfoDb::seqnames(seqinfo)]
  if( length(seqnames) == 0L ) {
    stop("No intersection between chromosome names in fasta, ",
         "annotation and seqinfo data.", 
         call. = FALSE)
  }
  return(seqinfo)
}

# Modifiertype -----------------------------------------------------------------

#'@importFrom stringr str_locate
# check if a class of type x exists
.norm_modifiertype <- function(x){
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
  nameId <- stringr::str_locate(class@className,"Mod")
  if(nrow(nameId) == 0L){
    stop("The string 'Mod' must be present once at the front of the class ",
         "name.",
         call. = FALSE)
  }
  if(nrow(nameId) > 1L || nameId[,"start"] != 1L || nameId[,"end"] != 3L){
    stop("The string 'Mod' can only be present once at the front of the class ",
         "name.",
         call. = FALSE)
  }
  x
}

.norm_mod <- function(mod, className){
  f <- mod %in% Modstrings::shortName(Modstrings::ModRNAString())
  if(length(which(f)) != length(mod)){
    stop("Modification '",mod[!f],"' as defined for ",className," does not ",
         "exist in the Modstrings dictionary for modified RNA sequences.",
         call. = FALSE)
  }
  mod
}

.norm_modifications <- function(ans, args){
  # ToDo: implement sanity check for given modifications
  ans
}
