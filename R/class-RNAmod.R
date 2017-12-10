#' @include AllGenerics.R
#' @include AllGenerics-export.R
NULL

#' @name RNAmod-class
#' 
#' @title RNAmod Class
#' 
#' @description 
#' 
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
setClass("RNAmod",
         slots = c(
           experimentName = "character",
           inputFile = "character",
           inputGFF = "character",
           inputFasta = "character",
           .dataExperiments = "data.frame",
           .dataFasta = "FaFile",
           .dataGFF = "GRanges",
           .wd = "character",
           .inputFolder = "character",
           .outputFolder = "character"
         ),
         prototype=list(
           .dataExperiments = data.frame(),
           .dataGFF = GenomicRanges::GRanges(),
           .wd = "",
           .inputFolder = "data/",
           .outputFolder = "results/"
         )
)


#' @rdname RNAmod-class
#'
#'
#' @export
#'
#' @import Rsamtools
#'
#' @examples
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
                        experimentFasta, 
                        experimentGFFuORF = "",
                        outputExists) {
    requireNamespace("rtracklayer", quietly = TRUE)
    
    # Save data to slots
    .Object@experimentName <- experimentName
    .Object@.wd <- paste0("./", .Object@experimentName, "/")
    
    .Object@inputFile <- paste0(.Object@.wd, 
                                .Object@.inputFolder, 
                                experimentLayout)
    .Object@inputGFF <- paste0(.Object@.wd, 
                               .Object@.sourceFolder, 
                               experimentGFF)
    .Object@inputFasta <- paste0(.Object@.wd, 
                                 .Object@.sourceFolder, 
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
    .Object@.dataGFF <- import.gff3(.Object@inputGFF)
    if( any( names(Rsamtools::scanFa(.Object@.dataFasta)) != unique(chrom(.Object@.dataGFF)) ) ) {
      stop("Names of sequences in the fasta file do not match the chromosome names in the gff file")
    }
    
    #load experiments
    experiments <- .check_experiment_data( read.csv(.Object@inputFile, 
                                                    header = TRUE, 
                                                    stringsAsFactors = FALSE, 
                                                    sep = ";"),
                                           getInputFolder(.Object))
    
    .Object@.dataExperiments <- experiments
    
    return(.Object)
  }
)

.check_for_column <- function(experimentData, column){
  # Does the column exist?
  if( !(column %in% colnames(experimentData) ) ) return(FALSE)
  # are they all NA?
  if( all( is.na(experimentData[,column]) ) ) return(FALSE)
  return(TRUE)
}

# Checks experimental data for validity
.check_experiment_data <- function(experiments, inputFolder){
  standardCols <- c("ExperimentNo","ExperimentName","Sample",
                    "ShortName","Replicate","Condition","BamFile")
  # Names of the first seven cols will be overwritten and are expected to fit
  # to the data type
  cols <- colnames(experiments)
  if( length(cols) < 7 ) stop("Please check the format of the csv file. ",
                              "It should be ; delimited and contain at ",
                              "least the following columns: ",
                              paste(standardCols, collapse = ","))
  cols[1:7] <- standardCols
  colnames(experiments) <- cols
  
  # Check if ExperimentNo and ExperimentName match
  tmpExperimentNo <- as.numeric(unique(experiments$ExperimentNo))
  tmpExperimentName <- as.character(unique(experiments$ExperimentName))
  
  if( !assertive::are_same_length(tmpExperimentNo, tmpExperimentName)){
    stop("ExperimentNo and ExperimentName are ambigeous. Check layout ",
         "file for errors in these columns.")
  }
  
  # check for conditions
  checkTypes <- c("Control","Treated")
  if( !all(experiments$Condition %in% checkTypes)) 
    stop("Conditions are not set correctly. ",
         "Following types are supported: ",
         paste(checkTypes, collapse = ", "))
  
  # Check inputs
  assertive::assert_all_are_positive(experiments$ExperimentNo)
  assertive::assert_all_are_whole_numbers(experiments$ExperimentNo)
  assertive::assert_all_are_non_missing_nor_empty_character(experiments$ExperimentName)
  assertive::assert_all_are_non_missing_nor_empty_character(experiments$Sample)
  assertive::assert_all_are_non_missing_nor_empty_character(experiments$ShortName)
  assertive::assert_all_are_positive(experiments$Replicate)
  assertive::assert_all_are_whole_numbers(experiments$Replicate)
  assertive::assert_all_are_existing_files(experiments$BamFile)
  
  return(experiments)
}

# matches the names in fasta file to the chromosome identifier in the GFF file
# to allow easy subsetting with GenomicRanges
.matchFastaToGFF <- function(inputFasta, 
                             inputGFF){
  assertive::assert_all_are_existing_files(c(inputFasta,inputGFF))
  
  fsa <- Biostrings::readDNAStringSet(inputFasta)
  gff <- rtracklayer::import.gff(inputGFF)
  
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
  if(names(fsa) == unique(rtracklayer::chrom(gff))){
    return(inputFasta)
  }
  
  names(fsa) <- unique(rtracklayer::chrom(gff))
  
  fileName <- gsub(".fsa","_fixed.fsa",inputFasta)
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
#' getParamFolder(experiment)
#' getSourceFolder(experiment)
#' getExperimentName(experiment)
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