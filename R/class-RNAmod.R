#' @include AllGenerics.R
#' @include AllGenerics-export.R
NULL

#' @name RNAmod-class
#' 
#' @title RNAmod Class
#' 
#' @description 
#' \code{RNAmod}
NULL

#' @title RNAmod-class
#' 
#' @rdname RNAmod-class 
#'
#'
#' @return a RNAmod object
#' @export
#' 
#' @importClassesFrom GenomicRanges GRanges
#'
#' @examples
setClass("RNAmod",
         slots = c(
           experimentName = "character",
           inputFile = "character",
           inputGFF = "character",
           inputFasta = "character",
           .dataSamples = "data.frame",
           .dataFasta = "FaFile",
           .dataGFF = "GRanges",
           .wd = "character",
           .inputFolder = "character",
           .outputFolder = "character",
           .mapQuality = "numeric"
         ),
         prototype=list(
           .dataSamples = data.frame(),
           .dataGFF = GenomicRanges::GRanges(),
           .wd = "",
           .inputFolder = "data/",
           .outputFolder = "results/",
           .mapQuality = 20
         )
)


#' @rdname RNAmod-class
#'
#' @import Rsamtools
#'
#' @export
RNAmod <- function(experimentName, 
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
  
  class <- new("RNAmod", 
               experimentName, 
               experimentLayout, 
               experimentGFF, 
               experimentFasta)
  return(class)
}

setMethod(
  f = "initialize", 
  signature = signature(.Object = "RNAmod"),
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
    .Object@.dataGFF <- rtracklayer::import.gff3(.Object@inputGFF)
    if( any( names(Rsamtools::scanFa(.Object@.dataFasta)) != unique(rtracklayer::chrom(.Object@.dataGFF)) ) ) {
      stop("Names of sequences in the fasta file do not match the chromosome names in the gff file")
    }
    
    #load sample data
    samples <- .check_sample_data( read.csv(.Object@inputFile, 
                                            header = TRUE, 
                                            stringsAsFactors = FALSE, 
                                            sep = ";"),
                                   getInputFolder(.Object))
    
    .Object@.dataSamples <- samples
    
    return(.Object)
  }
)

.check_for_column <- function(sampleData, column){
  # Does the column exist?
  if( !(column %in% colnames(sampleData) ) ) return(FALSE)
  # are they all NA?
  if( all( is.na(sampleData[,column]) ) ) return(FALSE)
  return(TRUE)
}

# Checks experimental data for validity
.check_sample_data <- function(samples, inputFolder){
  standardCols <- c("ExperimentNo","SampleName",
                    "ShortName","Replicate","BamFile")
  # Names of the first seven cols will be overwritten and are expected to fit
  # to the data type
  cols <- colnames(samples)
  if( length(cols) < length(standardCols) ) {
    stop("Please check the format of the csv file. It should be ; delimited ",
         "and contain at least the following columns: ",
         paste(standardCols, collapse = ","))
  }
  cols[1:length(standardCols)] <- standardCols
  colnames(samples) <- cols
  
  # Check if ExperimentNo and ExperimentName match
  tmpExperimentNo <- as.numeric(samples$ExperimentNo)
  tmpSampleName <- paste0(as.character(samples$SampleName),"_",as.character(samples$Replicate))
  
  if( !assertive::are_same_length(tmpExperimentNo, tmpSampleName)){
    stop("ExperimentNo and SampleName/Replicate are ambigeous. Check layout ",
         "file for errors in these columns.")
  }
  
  # Check inputs
  assertive::assert_all_are_positive(samples$ExperimentNo)
  assertive::assert_all_are_whole_numbers(samples$ExperimentNo)
  assertive::assert_all_are_non_missing_nor_empty_character(samples$SampleName)
  assertive::assert_all_are_non_missing_nor_empty_character(samples$ShortName)
  assertive::assert_all_are_positive(samples$Replicate)
  assertive::assert_all_are_whole_numbers(samples$Replicate)
  assertive::assert_all_are_existing_files(paste0(inputFolder,samples$BamFile))
  
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

#' @rdname RNAmod-class
#' 
#' @param object a RNAmod object.
setMethod(
  f = "show",
  signature = signature(object = "RNAmod"),
  definition = function(object) {
    print("test")
  }
)



# getters and setters ----------------------------------------------------------

#' @rdname RNAmod-Accessors
#'
#' @export
#'
#' @examples
#' \donttest{
#' data("RNAmod_test", package = "RNAmod")
#'
#' getInputFolder(experiment)
#' getOutputFolder(experiment)
#' getExperimentName(experiment)
#' head(getGFF(experiment))
#' getFasta(experiment)
#' }
setMethod(
  f = "getInputFolder",
  signature = signature(.Object = "RNAmod"),
  definition = function(.Object){
    return(paste0(.Object@.wd, .Object@.inputFolder))
  }
)
#' @rdname RNAmod-Accessors
#'
#' @export
setMethod(
  f = "getOutputFolder",
  signature = signature(.Object = "RNAmod"),
  definition = function(.Object){
    return(paste0(.Object@.wd, .Object@.outputFolder))
  }
)
#' @rdname RNAmod-Accessors
#'
#' @export
setMethod(
  f = "getExperimentName",
  signature = signature(.Object = "RNAmod"),
  definition = function(.Object){
    return(.Object@experimentName)
  }
)
#' @rdname RNAmod-Accessors
#'
#' @export
setMethod(
  f = "getGFF",
  signature = signature(.Object = "RNAmod"),
  definition = function(.Object){
    return(.Object@.dataGFF)
  }
)
#' @rdname RNAmod-Accessors
#'
#' @export
setMethod(
  f = "getFasta",
  signature = signature(.Object = "RNAmod"),
  definition = function(.Object){
    return(.Object@.dataFasta)
  }
)