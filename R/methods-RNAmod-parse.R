#' @include RNAmod.R
NULL


#' @name analyzeModifications
#' 
#' @title analyzeModifications
#'
#' @param .Object RNAmod. 
#' @param number numeric. 
#' @param modifications character. 
#'
#' @return
#' @export
#' 
#' @importFrom Rsamtools scanBam
#'
#' @examples
setMethod(
  f = "analyzeModifications", 
  signature = signature(.Object = "RNAmod",
                        number = "numeric",
                        modifications = "character"),
  definition = function(.Object,
                        number,
                        modifications){
    # Check input
    assertive::is_a_number(number)
    assertive::assert_all_are_non_empty_character(modifications)
    
    # get experiment data and extract bam file path
    experiment <- getExperimentData(.Object,number)
    
    # load classes for modification analysis
    modClasses <- vector(mode = "list", length = length(modifications))
    for(i in seq_along(modifications)){
      className <- paste0("mod_",modifications[[i]])
      tryCatch(
        class <- new(className),
        error = function(e) stop("Class for detecting ",modifications[[1]],
                                 " does not exist (",className,").")
      )
      if( !existsMethod("parseMod",class(class)) )
        stop("Function parseMod() not defined for ",class(class))
      
      modClasses[[i]] <- class
    }
    names(modClasses) <- modifications
    
    # extract gene boundaries
    gff <- .Object@.dataGFF
    gff <- gff[S4Vectors::mcols(gff)$type %in% RNAMOD_MOD_CONTAINING_FEATURES,]
    # combine path and file name
    files <- paste0(getInputFolder(.Object),experiment$BamFile)
    
    # assemble param for scanBam
    param <- .assemble_scanBamParam(gff, 
                                    .Object@.mapQuality,
                                    .get_acceptable_chrom_ident(files))
    
    data <- lapply(files,
                   FUN = .detect_mod,
                   gff,
                   .Object@.dataFasta,
                   param,
                   modClasses)
    
    
  }
)


.detect_mod <- function(bamFile,
                        gff,
                        fasta,
                        param,
                        mods){
  # Construct Dataframe from scanBam data
  bamData <- Rsamtools::scanBam(bamFile, param=param)
  bamData <- .convert_bam_to_DataFrame(bamData, param=param)
  
  # process result
  bamData <- S4Vectors::split(bamData,
                              bamData$ID)
  
  # for testing
  browser()
  bamData <- bamData[grepl("RDN18",names(bamData))]
  
  res <- BiocParallel::bplapply(bamData,
                                FUN = .detect_mod_in_transcript,
                                gff,
                                mods)
  
  .save_mod_result(res)
}

.detect_mod_in_transcript <- function(data,gff,mods){
  ID <- unique(data$ID)
  resData <- vector(mode="list",length = length(mods))
  for(i in seq_along(mods)){
    resData[[i]] <- parseMod(mods[[i]],data)
  }
  names(resData) <- names(mods)
  
  resGRange <-.construct_GRanges_from_mod_result(res,ID,gff)
  res <- list(data = resData,
              GRange = resGRange)
  return(res)
}

.construct_GRanges_from_mod_result <- function(res,ID,gff){
  
}

.save_mod_result <- function(res){
  
}