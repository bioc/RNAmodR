#' @include AllGenerics.R
#' @include AllGenerics-export.R
NULL

#' @name RNAmodR-class
#' 
#' @title RNAmodR Class
#' 
#' @description 
#' \code{RNAmodR}
NULL

#' @title RNAmodR-class
#' 
#' @rdname RNAmodR-class 
#'
#' @description 
#' title
#'
#' @return a RNAmodR object
#' @export
#' 
#' @import Rsamtools
#' @importClassesFrom GenomicRanges GRanges
#'
#' @examples
#' \donttest{
#' unzip(system.file("extdata", package = "RNAmodR", file = "RNAmodR.zip" ))
#' mod <- RNAmodR("test",
#'                "test_layout.csv",
#'                "test_gff.gff3",
#'                "test_masked.fasta")
#' }
setClass("RNAmodR",
         slots = c(
           experimentName = "character",
           inputFile = "character",
           inputGFF = "character",
           inputFasta = "character",
           .dataSamples = "DataFrame",
           .dataFasta = "FaFile",
           .dataGFF = "GRanges",
           .wd = "character",
           .inputFolder = "character",
           .outputFolder = "character",
           .mapQuality = "numeric"
         ),
         prototype=list(
           .dataSamples = S4Vectors::DataFrame(),
           .dataGFF = GenomicRanges::GRanges(),
           .wd = "",
           .inputFolder = "data/",
           .outputFolder = "results/",
           .mapQuality = 20
         )
)


#' @rdname RNAmodR-class
#'
#' @importFrom utils read.csv
#'
#' @export
RNAmodR <- function(experimentName, 
                    experimentLayout, 
                    experimentGFF, 
                    experimentFasta){
  
  # Input check
  if(missing(experimentName)) stop("Experimental name not provided")
  if(missing(experimentLayout)) stop("Experimental layout not provided")
  if(missing(experimentGFF)) stop("experimentGFF information not provided")
  if(missing(experimentFasta)) 
    stop("experimentFasta information not provided")
  
  if(!is.character(experimentName)) 
    stop("experimentName can only be assigned a character vector of length 1")
  if(!is.character(experimentLayout)) 
    stop("experimentLayout can only be assigned a character vector of length 1")
  if(!is.character(experimentGFF)) 
    stop("experimentGFF can only be assigned a character vector of length 1")
  if(!is.character(experimentFasta)) 
    stop("experimentFasta can only be assigned a character vector of length 1")
  
  class <- new("RNAmodR", 
               experimentName, 
               experimentLayout, 
               experimentGFF, 
               experimentFasta)
  return(class)
}

setMethod(
  f = "initialize", 
  signature = signature(.Object = "RNAmodR"),
  definition = function(.Object, 
                        experimentName, 
                        experimentLayout, 
                        experimentGFF, 
                        experimentFasta) {
    requireNamespace("rtracklayer", quietly = TRUE)
    # Save data to slots
    .Object@experimentName <- experimentName
    .Object@.wd <- paste0("./", .Object@experimentName, "/")
    
    .Object@inputFile <- paste0(.Object@.wd, 
                                .Object@.inputFolder, 
                                experimentLayout)
    .Object@inputGFF <- paste0(.Object@.wd, 
                               .Object@.inputFolder, 
                               experimentGFF)
    .Object@inputFasta <- paste0(.Object@.wd, 
                                 .Object@.inputFolder, 
                                 experimentFasta)
    # check for existing and readable files
    assertive::assert_all_are_existing_files(c(.Object@inputFile,
                                               .Object@inputGFF,
                                               .Object@inputFasta))
    # Load fasta and gff file
    .Object@inputFasta <- .matchFastaToGFF(.Object@inputFasta, 
                                           .Object@inputGFF)
    .Object@.dataFasta <- Rsamtools::FaFile(.Object@inputFasta)
    Rsamtools::indexFa(.Object@inputFasta)
    .Object@.dataGFF <- rtracklayer::import.gff3(.Object@inputGFF,
                                                 colnames = RNAMODR_GFF_COLNAMES)
    if( any( names(Rsamtools::scanFa(.Object@.dataFasta)) != unique(rtracklayer::chrom(.Object@.dataGFF)) ) ) {
      stop("Names of sequences in the fasta file do not match the chromosome names in the gff file")
    }
    #load sample data
    samples <- .check_sample_data(utils::read.csv(.Object@inputFile, 
                                                  header = TRUE, 
                                                  stringsAsFactors = FALSE, 
                                                  sep = ";",
                                                  fileEncoding = "UTF-8-BOM"),
                                  getInputFolder(.Object))
    .Object@.dataSamples <- samples
    return(.Object)
  }
)

# Checks experimental data for validity
.check_sample_data <- function(samples, inputFolder){
  samples <- S4Vectors::DataFrame(samples)
  # Names of the first seven cols will be overwritten and are expected to fit
  # to the data type
  cols <- colnames(samples)
  if( length(cols) < length(RNAMODR_DEFAULT_COLNAMES) |
      !all(cols[1:length(RNAMODR_DEFAULT_COLNAMES)] == RNAMODR_DEFAULT_COLNAMES) ) {
    stop("Please check the format of the csv file. It should be ; delimited ",
         "and contain at least the following columns: ",
         paste(RNAMODR_DEFAULT_COLNAMES, collapse = ","),
         call. = FALSE)
  }
  cols[1:length(RNAMODR_DEFAULT_COLNAMES)] <- RNAMODR_DEFAULT_COLNAMES
  colnames(samples) <- cols
  # Check if ExperimentNo and ExperimentName match
  tmpExperimentNo <- as.numeric(samples$ExperimentNo)
  tmpSampleName <- paste0(as.character(samples$SampleName),
                          "_",
                          as.character(samples$Replicate))
  if( !assertive::are_same_length(tmpExperimentNo, tmpSampleName)){
    stop("ExperimentNo and SampleName/Replicate are ambigeous. Check layout ",
         "file for errors in these columns.",
         call. = FALSE)
  }
  # Check inputs
  assertive::assert_all_are_positive(samples$ExperimentNo)
  assertive::assert_all_are_whole_numbers(samples$ExperimentNo)
  assertive::assert_all_are_non_missing_nor_empty_character(samples$SampleName)
  assertive::assert_all_are_non_missing_nor_empty_character(samples$ShortName)
  assertive::assert_all_are_positive(samples$Replicate)
  assertive::assert_all_are_whole_numbers(samples$Replicate)
  assertive::assert_all_are_existing_files(paste0(inputFolder,samples$BamFile))
  assertive::assert_all_are_whole_numbers(samples$MapQuality)
  # check modifications
  samples$Modifications <- lapply(samples$Modifications,
                                  function(str){str_split(str, ",")})
  modifications <- unique(unlist(samples$Modifications))
  modClasses <- .load_mod_classes(modifications)
  # check for valid condition types
  types <- unique(samples$Conditions)
  if(!all(types %in% RNAMODR_SAMPLE_TYPES)){
    stop("Not all condition types are valid. Valid types are:\n",
         paste(RNAMODR_SAMPLE_TYPES, collapse = ","),
         call. = FALSE)
  }
  return(samples)
}

# matches the names in fasta file to the chromosome identifier in the GFF file
# to allow easy subsetting with GenomicRanges
.matchFastaToGFF <- function(inputFasta, 
                             inputGFF){
  assertive::assert_all_are_existing_files(c(inputFasta,inputGFF))
  
  fsa <- Biostrings::readDNAStringSet(inputFasta)
  gff <- rtracklayer::import.gff3(inputGFF)
  
  # number of chromosomes and sequences have to match
  if( !assertive::are_same_length(names(fsa), unique(rtracklayer::chrom(gff)))){
    stop("Fasta and GFF file don't have matching number of sequences/",
         "sequence identifiers. Please make sure that the number of fasta ",
         "sequences matches the number of unique chromosome identifiers in ",
         "the GFF file. The names of the fasta sequences will be overridden ",
         "with the chromosomal identifier from the GFF file.",
         call. = FALSE)
  }
  
  # if identifies already match skip renaming
  if( all(names(fsa) == unique(rtracklayer::chrom(gff)))){
    return(inputFasta)
  }
  
  names(fsa) <- unique(rtracklayer::chrom(gff))
  
  fileName <- gsub(".fsa","_fixed.fsa",inputFasta)
  fileName <- gsub(".fasta","_fixed.fasta",fileName)
  Biostrings::writeXStringSet(fsa, fileName)
  
  message(paste0("Fasta file ", fileName, " with modified sequence names",
                 " written."))
  return(fileName)
}


# default class methods --------------------------------------------------------

#' @rdname RNAmodR-class
#' 
#' @param object a RNAmod object.
setMethod(
  f = "show",
  signature = signature(object = "RNAmodR"),
  definition = function(object) {
    print("test")
  }
)


# getters and setters ----------------------------------------------------------

#' @rdname RNAmodR-Accessors
#'
#' @export
setMethod(
  f = "getInputFolder",
  signature = signature(.Object = "RNAmodR"),
  definition = function(.Object){
    return(paste0(.Object@.wd, .Object@.inputFolder))
  }
)
#' @rdname RNAmodR-Accessors
#'
#' @export
setMethod(
  f = "getOutputFolder",
  signature = signature(.Object = "RNAmodR"),
  definition = function(.Object){
    return(paste0(.Object@.wd, .Object@.outputFolder))
  }
)
#' @rdname RNAmodR-Accessors
#'
#' @export
setMethod(
  f = "getExperimentName",
  signature = signature(.Object = "RNAmodR"),
  definition = function(.Object){
    return(.Object@experimentName)
  }
)
#' @rdname RNAmodR-Accessors
#'
#' @export
setMethod(
  f = "getGFFFile",
  signature = signature(.Object = "RNAmodR"),
  definition = function(.Object){
    return(.Object@.dataGFF)
  }
)
#' @rdname RNAmodR-Accessors
#'
#' @export
setMethod(
  f = "getFastaFile",
  signature = signature(.Object = "RNAmodR"),
  definition = function(.Object){
    return(.Object@.dataFasta)
  }
)