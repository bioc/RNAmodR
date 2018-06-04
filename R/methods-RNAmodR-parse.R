#' @include class-RNAmodR.R
#' @include class-RNAmodR-args.R
NULL

#' @rdname parseForModifications
#' 
#' @title parseForModifications
#' 
#' @description 
#' \code{parseForModifications} uses the input information from GRanges, Fafile
#' and bam input files to look for known patterns of RNA post-transcriptional
#' modifications. 
#'
#' @param x a RNAmod object. 
#' @param number a number defining an experiment. 
#' @param args arguments for modification parsing as an RNAmodRargs object
#' @param name a character string for identification
#' @param gff a gff3 file name
#' @param fasta a fasta file name
#' @param files file names for the bam input files
#' @param conditions treated or control
#' 
#' @return TRUE if the analysis does finish without error
#' @export
#' 
#' @examples
#' \donttest{
#' unzip(system.file("extdata", package = "RNAmodR", file = "RNAmodR.zip" ))
#' mod <- RNAmodR("test",
#'                "test_layout.csv",
#'                "test_gff.gff3",
#'                "test_masked.fasta")
#' parseForModifications(mod,1)               
#' }
setMethod(
  f = "parseForModifications", 
  signature = signature(x = "RNAmodR",
                        number = "numeric",
                        param = "character"),
  definition = function(x,
                        number,
                        param){
    # Check input
    assertive::is_a_number(number)
    # get experiment name
    experiment <- getExperimentData(x,number)
    name <- unique(experiment$SampleName)
    # get gff annotation
    gff <- x@.dataGFF
    fasta <- x@.dataFasta
    # combine path and file name and create argument class
    files <- paste0(getInputFolder(x),experiment$BamFile)
    conditions <- experiment$Conditions
    args <- RNAmodRargs(param,
                        files,
                        conditions)
    # pass assembled arguments on
    res <- parseForModificationsWithArgs(x = args,
                                         name = name,
                                         gff = gff,
                                         fasta = fasta)
    browser()
    # save aggregated scores to file
    saveScores(x, number, res$scores)
    message("Saved scores as gff3 and csv file.")
    # Save found modifications as gff file
    saveModifications(x, number, res$modifications)
    message("Saved detected modifications as gff3 and csv file.")
    return(invisible(TRUE))
  }
)

#' @rdname parseForModifications
#' 
#' @export
setMethod(
  f = "parseForModificationsWithArgs", 
  signature = signature(x = "RNAmodRargs",
                        name = "character",
                        gff = "GRanges",
                        fasta = "FaFile"),
  definition = function(x,
                        name,
                        gff,
                        fasta){
    browser()
    # Input checks
    assertive::assert_is_a_non_missing_nor_empty_string(name)
    # a little startup message
    message("Searching for modifications in sample '",
            name,
            "'...")
    # subset to relevant annotations 
    gff_subset <- .get_parent_annotations(
      .subset_rnamod_containing_features(gff)
    )
    # assemble param for scanBam
    scanBamParam <- .assemble_scanBamParam(gff_subset, 
                                           .get_min_map_quality(),
                                           .get_acceptable_chrom_ident(x@files))
    # retrieve the quantifier and identifier needed
    identifier <- loadIdentifier(x)
    quantifier <- loadQuantifier(x)
    # get transcript quantification 
    data <- 
      sapply(quantifier,
             function(quant){
               quantifyReadData(quant,
                                x,
                                gff,
                                fasta,
                                scanBamParam)
    }, simplify = FALSE, USE.NAMES = TRUE)
    scores <- 
      sapply(identifier, 
             function(ident){
               scoreModifications(ident,
                                  subsetData(ident,
                                             data),
                                  x)
    }, simplify = FALSE, USE.NAMES = TRUE)
    modifications <- 
      sapply(identifier,
             function(ident){
               identifyModifications(ident,
                                     subsetResults(ident,
                                                   scores),
                                     x)
    }, simplify = FALSE, USE.NAMES = TRUE)
    return(list(scores = scores,
                modifications = modifications))
  }
)


#' 
#' #' @rdname calculateScoreTable
#' #' @export
#' setMethod(
#'   f = "calculateScoreTable", 
#'   signature = signature(x = "RNAmodR",
#'                         number = "numeric",
#'                         param = "character"),
#'   definition = function(x,
#'                         number,
#'                         param){
#'   }
#' )