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
    
    # get experiment data and extract bam file path
    experiment <- getExperimentData(.Object,number)
    fileName <- paste0(getOutputFolder(.Object),
                       "RNAmod_",
                       unique(experiment["SampleName"]),
                       "_",
                       paste(modifications,collapse = "-"),
                       ".gff3")
    message("Analyzing modifications for Sample '",
            unique(experiment["SampleName"]),
            "'...")
    
    modClasses <- .load_mod_classes(modifications)
    message("Detecting modification types: '",
            paste(names(modClasses), 
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
                   FUN = .detect_mod,
                   gff,
                   .Object@.dataFasta,
                   param,
                   modClasses)
    
    # Merge data from all replicates
    # Modifications:
    mods <- lapply(modClasses,
    # mods <- BiocParallel::bplapply(modClasses,
                                   FUN = mergeModsOfReplicates,
                                   gff,
                                   .Object@.dataFasta,
                                   lapply(data,"[[","mods"))
    # Positions:
    positions <- lapply(modClasses,
    # positions <- BiocParallel::bplapply(modClasses,
                   FUN = mergePositionsOfReplicates,
                   gff,
                   .Object@.dataFasta,
                   lapply(data,"[[","positions"))
    
    
    # Construct DataFrame from found modifications
    df <- .construct_DataFrame_from_mod_result(mods,gff)
    
    # Save found modifications as gff file
    rtracklayer::export.gff3(GRanges(df), con = fileName)
    message("Saved detected modifications as gff3 file in: ",
            fileName)
    
    # Save found modifications asSummarizedExperiment
    se <- .construct_SE_from_mod_result(experiment,df,gff,positions)
    save(se, file = gsub(".gff",".RData",fileName))
    message("Saved detected modifications as SummarizedExperiment: ",
            gsub(".gff3",".RData",fileName))
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
    if( !existsMethod("parseMod",signature(class(class),
                                           "numeric",
                                           "GRanges",
                                           "DNAString",
                                           "DataFrame") ) )
      stop("Function parseMod() not defined for ",class(class))
    if( !existsMethod("mergeModsOfReplicates",signature(class(class),
                                                        "GRanges",
                                                        "FaFile",
                                                        "list") ) )
      stop("Function mergeModsOfReplicates() not defined for ",class(class))
    
    if( !existsMethod("convertReadsToPositions",signature(class(class),
                                                          "numeric",
                                                          "GRanges",
                                                          "DNAString",
                                                          "DataFrame") ) )
      stop("Function convertReadsToPositions() not defined for ",class(class))
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
.detect_mod <- function(bamFile,
                        gff,
                        fasta,
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
  bamData <- bamData[grepl("RDN18",names(bamData))]
  
  mods <- lapply(bamData,
                 # res <- BiocParallel::bplapply(bamData,
                 FUN = .detect_mod_in_transcript,
                 totalCounts,
                 gff,
                 fasta,
                 modClasses)
  names(mods) <- names(bamData)
  
  positions <- lapply(bamData,
                 # res <- BiocParallel::bplapply(bamData,
                 FUN = .get_positions_in_transcript,
                 totalCounts,
                 gff,
                 fasta,
                 modClasses)
  names(positions) <- names(bamData)
  return(list(mods = mods,
              positions = positions))
}

# For each transcript check for all modifications 
.detect_mod_in_transcript <- function(data,totalCounts,gff,fasta,mods){
  # get ID, GRanges and sequence
  ID <- unique(data$ID)
  gff <- gff[(is.na(S4Vectors::mcols(gff)$ID) & 
                S4Vectors::mcols(gff)$Name == ID) |
               (!is.na(S4Vectors::mcols(gff)$ID) & 
                  S4Vectors::mcols(gff)$ID == ID),]
  seq <- getSeq(fasta,gff)[[1]]
  
  # Parse reads for all modifications based on sequence
  resData <- vector(mode="list",length = length(mods))
  for(i in seq_along(mods)){
    resData[[i]] <- parseMod(mods[[i]],totalCounts,gff,seq,data)
  }
  names(resData) <- names(mods)
  return(resData)
}

# For each transcript get positional data
# This can be individually done for different modification types
.get_positions_in_transcript <- function(data,totalCounts,gff,fasta,mods){
  # get ID and GRanges
  ID <- unique(data$ID)
  gff <- gff[(is.na(S4Vectors::mcols(gff)$ID) & 
                S4Vectors::mcols(gff)$Name == ID) |
               (!is.na(S4Vectors::mcols(gff)$ID) & 
                  S4Vectors::mcols(gff)$ID == ID),]
  seq <- getSeq(fasta,gff)[[1]]
  
  # Parse reads for all modifications based on sequence
  resData <- vector(mode="list",length = length(mods))
  for(i in seq_along(mods)){
    resData[[i]] <- convertReadsToPositions(mods[[i]],totalCounts,gff,seq,data)
  }
  names(resData) <- names(mods)
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
              "RNAmod_p.value","RNAmod_p.value.sd","RNAmod_nReplicates",
              "RNAmod_inNReplicate")]
  
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
                  RNAmod_p.value.sd = vapply(gene,"[[",numeric(1),"p.value.sd"),
                  RNAmod_nReplicates = vapply(gene,"[[",numeric(1),"replicates"),
                  RNAmod_inNReplicate = (vapply(gene,"[[",numeric(1),"replicates")-vapply(gene,"[[",numeric(1),"undetectedInNReplicate")))
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
  genePos <- lapply(rownames,
                     function(name){
                       res <- lapply(modTypes, 
                                     function(modType){
                                       positions[[modType]][[name]]
                       })
                       names(res) <- modTypes
                       return(res)
                     })
  names(genePos) <- rownames
  SummarizedExperiment::rowData(se)$positions <- genePos
  
  return(se)
}
