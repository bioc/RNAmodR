#' @include RNAmod.R
NULL


#' @name parseForModifications
#' 
#' @title parseForModifications
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
  f = "parseForModifications", 
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
    
    
    message("Searching for modifications in sample '",
            unique(experiment["SampleName"]),
            "'...")
    
    # add the default modification class, which is just used for handling read 
    # positions (5'-end), but not modification detection.
    modClasses <- .load_mod_classes(append("default",modifications))
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
    data <- data[!is.null(data)]
    if(length(data) == 0){
      stop("No reads detected in any bam file",
           call. = FALSE)
    }
    
    # Merge data from all replicates and detect modifications:
    mod_positions <- vector(mode = "list", length = length(modClasses))
    names(mod_positions) <- names(modClasses)
    for(i in seq_along(modClasses)){
      mod_positions[[i]] <- parseMod(modClasses[[i]],
                                     gff,
                                     force(.Object@.dataFasta),
                                     data)
    }
    # General cleanup
    mod_positions <- mod_positions[!is.na(mod_positions)]
    # Specific cleanup
    mod_positions <- lapply(mod_positions, function(modPerTypes){
      modPerTypes[!vapply(modPerTypes,is.null,logical(1))]
    })
    # check if any modification could be detected
    nMods <- lapply(mod_positions, function(modPerTypes){
      sum(unlist(lapply(modPerTypes,length)))
    })
    if( sum(unlist(nMods)) == 0){
      stop("No modifications detected. Aborting...",
           call. = FALSE)
    }
    
    # Merge position data
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
    se <- .construct_SE_from_mod_result(experiment,gff,df,positions)
    setSummarizedExperiment(.Object,
                            se,
                            number,
                            modifications)
    message("Saved detected modifications as SummarizedExperiment.")
  }
)


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
  
  # browser()
  # for testing
  # m7G 
  # bamData <- bamData[names(bamData) %in% c("RDN18-1")]
  # bamData <- bamData[names(bamData) %in% c("tC(GCA)B")]
  
  # D 
  # bamData <- bamData[names(bamData) %in% c("tH(GUG)E1")]
  bamData <- bamData[names(bamData) %in% c("tI(AAU)B")]
  # bamData <- bamData[names(bamData) %in% c("tD(GUC)B")]
  
  # m3C
  # bamData <- bamData[names(bamData) %in% c("tS(CGA)C")]
  # bamData <- bamData[names(bamData) %in% c("RDN25-1")]
  
  # combination
  # bamData <- bamData[names(bamData) %in% c("RDN18-1",
  #                                          "tS(CGA)C",
  #                                          "tC(GCA)B")]
  # if( getOption("RNAmod_debug") ){
  #   bamData <- bamData[names(bamData) %in% getOption("RNAmod_debug_transcripts")]
  # }
  
  if(length(bamData) == 0){
    warning("No reads detected in bam file '",
            bamFile,
            "'")
    return(NULL)
  }
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
.get_positions_in_transcript <- function(data,totalCounts,gff,modClasses){
  # get ID and GRanges
  ID <- unique(data$ID)
  gff <- .subset_gff_for_unique_transcript(gff, ID)
  
  # Parse reads for all modifications based on sequence
  resData <- vector(mode="list",length = length(modClasses))
  names(resData) <- names(modClasses)
  for(i in seq_along(modClasses)){
    resData[[i]] <- convertReadsToPositions(modClasses[[i]],totalCounts,gff,data)
  }
  resData <- resData[!is.na(resData)]
  return(resData)
}

# Construct DataFrame from RNAmod results per modification from individual
# gene results
.construct_DataFrame_from_mod_result <- function(data,gff){
  modTypes <- names(data)
  l <- lapply(seq_along(modTypes), function(j){
    modType <- data[[j]]
    genes <- names(modType)
    genesDf <- lapply(modType,
                      .get_dataframe_per_gene)
    nMods <- unlist(lapply(genesDf, function(df){
      nrow(df)
    }))
    if( sum(nMods) == 0 ) return(NULL)
    
    # This way os selection matches to one construting the initial read
    # DataFrame
    g <- gff[as.character(S4Vectors::mcols(gff)$ID) %in% genes |
               (as.character(S4Vectors::mcols(gff)$Name) %in% genes &
                  is.na(as.character(S4Vectors::mcols(gff)$ID)))]
    S4Vectors::mcols(gff)[is.na(as.character(S4Vectors::mcols(gff)$ID))]$ID <- 
      S4Vectors::mcols(gff)[is.na(as.character(S4Vectors::mcols(gff)$ID))]$Name
    S4Vectors::mcols(g)$ID <- factor(S4Vectors::mcols(g)$ID, levels = genes)
    g <- g[order(S4Vectors::mcols(g)$ID)]
    
    strand <- as.character(BiocGenerics::strand(g))
    chrom <- as.character(GenomeInfoDb::seqnames(g))
    genesDf <- lapply(seq_along(genes), function(i){
      if(nrow(genesDf[[i]]) == 0) return(NULL)
      df <- genesDf[[i]]
      df$strand <- strand[[i]]
      df$chrom <- chrom[[i]]
      df$Parent <- genes[[i]]
      df$source <- rep("RNAmodR",nrow(df))
      df$type <- rep("RNAMOD",nrow(df))
      df$RNAmod_type <- rep(modTypes[j],nrow(df))
      df$score <- df$RNAmod_signal
      df$ID <- paste0(df$Parent,"_",df$ID)
      return(df)
    })
    genesDf <- genesDf[!vapply(genesDf, is.null, logical(1))]
    if(length(genesDf) == 0) return(NULL)
    res <- do.call(rbind,genesDf)
    res <- res[order(res$start),]
    return(res)
  })
  l <- l[!vapply(l, is.null, logical(1))]
  if(length(l) == 0){
    stop("No modifications detected. Aborting...",
         call. = FALSE)
  }
  df <- do.call(rbind,l)
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


.construct_SE_from_mod_result <- function(experiment,gff,data,positions){
  # rownames and colnames(modTypes)
  modTypes <- unique(data$RNAmod_type)
  rownames <- S4Vectors::mcols(gff)$ID
  rownames[is.na(rownames)] <- 
    S4Vectors::mcols(gff[is.na(S4Vectors::mcols(gff)$ID),])$Name
  
  # generate col- and rowdata
  colData <- experiment[, c(colnames(experiment)[1:3])]
  colData <- colData[1,]
  rowData <- gff
  
  # generate additional row data
  geneMods <- .get_mod_data_per_row(rownames, data)
  genePos <- .get_position_data_per_row(rownames, positions)
  
  # generate assay data
  assays <- .get_assay_data(modTypes, geneMods, rownames)
    
  # construct SE
  se <- SummarizedExperiment::SummarizedExperiment(assays = assays,
                                                   colData = colData,
                                                   rowData = rowData)
  # Modify SE to include DataFrame on modified positions
  SummarizedExperiment::rowData(se)$mods <- geneMods
  # Save read position data in SE
  SummarizedExperiment::rowData(se)$positions <- genePos
  
  return(se)
}

# return a modification result per rowname. can be an empty result
.get_mod_data_per_row <- function(rownames, data){
  mod <- lapply(rownames,
                function(name){
                  return(data[data$Parent == name,])
                })
  names(mod) <- rownames
  return(mod)
}

# return a position result per rowname. can be an empty result
.get_position_data_per_row <- function(rownames, positions){
  posTypes <- names(positions)
  pos <- lapply(rownames,
                function(name){
                  res <- lapply(posTypes, 
                                function(posType){
                                  positions[[posType]][[name]]
                                })
                  names(res) <- posTypes
                  return(res)
                })
  names(pos) <- rownames
  return(pos)
}

# generate a list of named lists per modification
.get_assay_data <- function(modTypes, geneMods, rownames){
  nrow <- length(rownames)
  assays <-  lapply(modTypes,
                    function(modType){
                      mods <- lapply(geneMods, 
                                     function(gene){
                                       gene[gene$RNAmod_type == modType,]
                                     })
                      
                      data <- vapply(mods,nrow,numeric(1))
                      return(matrix(data,
                                    nrow,
                                    dimnames = list(rownames)))
                    })         
  names(assays) <- modTypes
  assays <- S4Vectors::SimpleList(assays)
  return(assays)
}