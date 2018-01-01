#' @include RNAmod.R
NULL


#' @name analyzeModifications
#' 
#' @title analyzeModifications
#'
#' @param .Object a RNAmod object. 
#' @param number a number defining an experiment. 
#' @param modifications name of modification to be used for analysis. 
#'
#' @return
#' @export
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
    
    # get experiment data and create outout folders
    experiment <- getExperimentData(.Object,number)
    
    
    message("Analyzing modifications for Sample '",
            unique(experiment["SampleName"]),
            "'...")
    
    # add the default modification class, which is just used for handling read 
    # positions (5'-end), but not modification detection.
    modifications <- append("default",modifications)
    modClasses <- .load_mod_classes(modifications)
    message("Detecting modification types: '",
            paste(names(modClasses[names(modClasses) != "default"]), 
                  collapse = "', '"),
            "'")
    
    # extract gene boundaries
    gff <- .Object@.dataGFF
    gff <- gff[S4Vectors::mcols(gff)$type %in% RNAMOD_MOD_CONTAINING_FEATURES,]
    # combine path and file name
    files <- paste0(getInputFolder(.Object),experiment$BamFile)
    
    # assemble param for scanBam
    param <- .assemble_scanBamParam(gff, 
                                    .Object@.mapQuality,
                                    .get_acceptable_chrom_ident(files))
    
    # detect modifications in each file
    data <- lapply(files,
                   FUN = .get_positions,
                   gff,
                   param,
                   modClasses)
    
    # Merge data from all replicates and detect modifications:
    mod_positions <- vector(mode = "list", length = length(modClasses))
    names(mod_positions) <- names(modClasses)
    
    for(i in seq_along(modClasses)){
      mod_positions[[i]] <- parseMod(modClasses[[i]],
                                     gff,
                                     .Object@.dataFasta,
                                     data)
    }
    mod_positions <- mod_positions[!is.na(mod_positions)]
    
    positions <- vector(mode = "list", length = length(modClasses))
    names(positions) <- names(modClasses)
    for(i in seq_along(modClasses)){
      positions[[i]] <- mergePositionsOfReplicates(modClasses[[i]],
                                                   gff,
                                                   .Object@.dataFasta,
                                                   data)
    }
    positions <- positions[!is.na(positions)]
    
    # Construct DataFrame from found modifications
    df <- .construct_DataFrame_from_mod_result(mod_positions,gff)
    
    # Save found modifications as gff file
    setGff(.Object,
           GRanges(df),
           number,
           modifications)
    message("Saved detected modifications as gff3 file.")
    
    # Save found modifications asSummarizedExperiment
    se <- .construct_SE_from_mod_result(experiment,df,gff,positions)
    setSummarizedExperiment(.Object,
                            se,
                            number,
                            modifications)
    message("Saved detected modifications as SummarizedExperiment.")
  }
)

# load classes for modification analysis
.load_mod_classes <- function(modifications){
  
  modClasses <- vector(mode = "list", length = length(modifications))
  for(i in seq_along(modifications)){
    className <- paste0("mod_",modifications[[i]])
    
    # try to create modification detection classes
    tryCatch(
      class <- new(className),
      error = function(e) stop("Class for detecting ",
                               modifications[[1]],
                               " does not exist (",className,").",
                               call. = FALSE)
    )
    if( !existsMethod("convertReadsToPositions",signature(class(class),
                                                          "numeric",
                                                          "GRanges",
                                                          "DataFrame") ) )
      stop("Function convertReadsToPositions() not defined for ",class(class))
    if( !existsMethod("parseMod",signature(class(class),
                                           "GRanges",
                                           "FaFile",
                                           "list") ) )
      stop("Function parseMod() not defined for ",class(class))
    if( !existsMethod("mergePositionsOfReplicates",signature(class(class),
                                                             "GRanges",
                                                             "FaFile",
                                                             "list") ) )
      stop("Function mergePositionsOfReplicates() not defined for ",class(class))
    modClasses[[i]] <- class
  }
  names(modClasses) <- modifications
  
  return(modClasses)
}


# detect modifications in each file
.get_positions <- function(bamFile,
                           gff,
                           param,
                           modClasses){
  # Construct Dataframe from scanBam data
  bamData <- Rsamtools::scanBam(bamFile, param=param)
  bamData <- .convert_bam_to_DataFrame(bamData, param=param)
  
  # Total counts
  totalCounts <- Rsamtools::idxstatsBam(bamFile, param=param)
  totalCounts <- sum(totalCounts$mapped)
  
  # process result
  bamData <- S4Vectors::split(bamData,
                              bamData$ID)
  
  # for testing
  # bamData <- bamData[names(bamData) %in% c("RDN18")]
  bamData <- bamData[names(bamData) %in% c("tS(CGA)C")]
  # bamData <- bamData[names(bamData) %in% c("RDN18","tS(CGA)C")]
  
  positions <- lapply(bamData,
  # res <- BiocParallel::bplapply(bamData,
                      FUN = .get_positions_in_transcript,
                      totalCounts,
                      gff,
                      modClasses)
  names(positions) <- names(bamData)
  return(positions)
}

# For each transcript get positional data
# This can be individually done for different modification types
.get_positions_in_transcript <- function(data,totalCounts,gff,mods){
  # get ID and GRanges
  ID <- unique(data$ID)
  gff <- gff[(is.na(S4Vectors::mcols(gff)$ID) & 
                S4Vectors::mcols(gff)$Name == ID) |
               (!is.na(S4Vectors::mcols(gff)$ID) & 
                  S4Vectors::mcols(gff)$ID == ID),]
  
  # Parse reads for all modifications based on sequence
  resData <- vector(mode="list",length = length(mods))
  names(resData) <- names(mods)
  for(i in seq_along(mods)){
    resData[[i]] <- convertReadsToPositions(mods[[i]],totalCounts,gff,data)
  }
  resData <- resData[!is.na(resData)]
  return(resData)
}

# Construct DataFrame from RNAmod results per modification from individual
# gene results
.construct_DataFrame_from_mod_result <- function(data,gff){
  modTypes <- names(data)
  for(j in seq_along(modTypes)){
    modType <- data[[j]]
    genes <- names(modType)
    genesDf <- lapply(modType,
                      .get_dataframe_per_gene)
    
    g <- gff[S4Vectors::mcols(gff)$ID %in% genes |
               S4Vectors::mcols(gff)$Name %in% genes]
    strand <- as.character(strand(g))
    chrom <- as.character(seqnames(g))
    
    for(i in seq_along(genes)){
      browser()
      genesDf[[i]]$strand <- strand[[i]]
      genesDf[[i]]$chrom <- chrom[[i]]
      genesDf[[i]]$Parent <- genes[[i]]
      genesDf[[i]]$source <- rep("RNAmodR",nrow(genesDf[[i]]))
      genesDf[[i]]$type <- rep("RNAMOD",nrow(genesDf[[i]]))
      genesDf[[i]]$RNAmod_type <- rep(modTypes[j],nrow(genesDf[[i]]))
      genesDf[[i]]$score <- genesDf[[i]]$RNAmod_signal
      genesDf[[i]]$ID <- paste0(genesDf[[i]]$Parent,"_",genesDf[[i]]$ID)
    }
  }
  
  df <- do.call(rbind,genesDf)
  df <- df[,c("chrom","start","end","strand","source","type","score","ID",
              "Parent","RNAmod_type","RNAmod_signal","RNAmod_signal_sd",
              "RNAmod_p.value","RNAmod_nbReplicates")]
  
  df$source <- factor(df$source)
  df$type <- factor(df$type)
  return(df)
}

# Construct DataFrame from RNAmod results per gene
.get_dataframe_per_gene <- function(gene){
  df <- S4Vectors::DataFrame(start = vapply(gene,"[[",numeric(1),"location"),
                  end = vapply(gene,"[[",numeric(1),"location"),
                  ID = names(gene),
                  RNAmod_signal = vapply(gene,"[[",numeric(1),"signal"),
                  RNAmod_signal_sd = vapply(gene,"[[",numeric(1),"signal.sd"),
                  RNAmod_p.value = vapply(gene,"[[",numeric(1),"p.value"),
                  RNAmod_nbReplicates = vapply(gene,"[[",numeric(1),"nbsamples"))
  return(df)
}


.construct_SE_from_mod_result <- function(experiment,data,gff,positions){
  # col- and rownames
  modTypes <- unique(data$RNAmod_type)
  rownames <- S4Vectors::mcols(gff)$ID
  rownames[is.na(rownames)] <- S4Vectors::mcols(gff[is.na(S4Vectors::mcols(gff)$ID),])$Name
  
  # dims
  nrow <- length(gff)
  ncol <- 1
  
  # assays
  geneMods <- lapply(rownames,
                     function(name){
                       return(data[data$Parent == name,])
                     })
  assays <-  lapply(modTypes,
                    function(modType){
                      mods <- lapply(geneMods, 
                                     function(gene){
                                       gene[gene$RNAmod_type == modType,]
                                     })
                      
                      data <- vapply(mods,nrow,numeric(1))
                      return(matrix(data,
                                    nrow,
                                    dimnames = list(rownames,
                                                    modType)))
                    })         
  names(assays) <- modTypes
  assays <- S4Vectors::SimpleList(assays)
  
  # col- and rowdata
  colData <- experiment[, c(colnames(experiment)[1:3])]
  colData <- colData[!duplicated(colData),]
  rowData <- gff
  
  # construct SE
  se <- SummarizedExperiment::SummarizedExperiment(assays = assays,
                                                   colData = colData,
                                                   rowData = rowData)
  # Modify SE to include DataFrame on modified positions
  SummarizedExperiment::rowData(se)$mods <- geneMods
  
  # Save read position data in SE
  posTypes <- names(positions)
  genePos <- lapply(rownames,
                     function(name){
                       res <- lapply(posTypes, 
                                     function(posType){
                                       positions[[posType]][[name]]
                       })
                       names(res) <- posTypes
                       return(res)
                     })
  names(genePos) <- rownames
  SummarizedExperiment::rowData(se)$positions <- genePos
  
  return(se)
}
