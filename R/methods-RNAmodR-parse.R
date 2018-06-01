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
    if(is.null(res)) return(NULL)
    # Save found modifications as gff file
    setGffResult(x,
                 res$gr,
                 number)
    message("Saved detected modifications as gff3 file.")
    setSummarizedExperiment(x,
                            res$se,
                            number)
    message("Saved detected modifications as SummarizedExperiment.")
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
    param <- .assemble_scanBamParam(gff_subset, 
                                    mapQuality,
                                    .get_acceptable_chrom_ident(files))
    # retrieve the quantifier and identifier needed
    quantifier <- loadQuantifier(x)
    identifier <- loadIdentifier(x)
    # get transcript quantification
    data <- sapply(quantifier, function(quant){
      quantifyReadData(quant,
                       x,
                       gff,
                       param)
    }, simplify = FALSE, USE.NAMES = TRUE)
    data <- sapply(identifier, function(quant){
      identifyModification(ident,
                           x,
                           subsetData(ident,
                                      data),
                           gff,
                           param)
    }, simplify = FALSE, USE.NAMES = TRUE)
    
    
    
    # retrieve the analysis types need
    modClasses <- .load_mod_classes(modifications)
    analysisTypes <- stats::setNames(vapply(modClasses, 
                                            getDataType, 
                                            character(1)),
                                     modifications)
    # group the analysis types
    analysisGroups <- stats::setNames(lapply(unique(analysisTypes), 
                                             function(x){
                                               names(analysisTypes[analysisTypes == x])
                                             }),unique(analysisTypes))
    # load the analysis classes
    analysisClasses <- .load_analysis_classes(names(analysisGroups))
   
    
    message(Sys.time(), ": Getting read data...")
    # load data into each analysis class
    analysisClasses <- sapply(names(analysisGroups), function(className){
      convertReadsToData(analysisClasses[[className]],
                         files,
                         conditions,
                         gff,
                         param)
    }, simplify = FALSE, USE.NAMES = TRUE)
    # parse data in each analysis class for the subset of modifications
    # this merges data from all replicates for the analysis
    message(Sys.time(), ": Parsing for modifications...")
    analysisClasses <- sapply(names(analysisGroups), function(className){
      modClassesSubset <- modClasses[vapply(modClasses, 
                                            getDataType, 
                                            character(1)) == className]
      parseMod(analysisClasses[[className]],
               gff,
               fasta,
               modClassesSubset)
    }, simplify = FALSE, USE.NAMES = TRUE)
    # check if any modification could be detected
    nMods <- vapply(names(analysisGroups), function(className){
      length(getModifications(analysisClasses[[className]]))
    }, numeric(1))
    if( sum(nMods) == 0){
      message("No modifications detected. Aborting...")
      return(NULL)
    }
    # Merge position data
    message(Sys.time(), ": Summarizing data...")
    analysisClasses <- sapply(names(analysisGroups), function(className){
      mergeDataOfReplicates(analysisClasses[[className]])
    }, simplify = FALSE, USE.NAMES = TRUE)
    # Retrieve position data
    positions <- sapply(names(analysisGroups), function(className){
      getPositions(analysisClasses[[className]])
    }, simplify = FALSE)
    # Retrieve modification data
    mod_positions <- sapply(names(analysisGroups), function(className){
      getModifications(analysisClasses[[className]])
    }, simplify = FALSE)
    # Construct DataFrame from found modifications
    df <- .construct_DataFrame_from_mod_result(mod_positions,
                                               gff_subset)
    gr <- GRanges(df)
    # Save found modifications asSummarizedExperiment
    se <- .construct_SE_from_mod_result(name,
                                        files,
                                        gff_subset,
                                        df,
                                        positions)
    return(list(gr = gr,
                se = se))
  }
)

# Construct DataFrame from RNAmod results per modification from individual
# gene results
.construct_DataFrame_from_mod_result <- function(data,
                                                 gff){
  l <- lapply(seq_along(data), function(i){
    analysisType <- names(data[i])
    analysisTypeData <- data[[i]]
    
    genes <- names(analysisTypeData)
    genesDf <- lapply(analysisTypeData,
                      .get_dataframe_per_gene)
    nMods <- unlist(lapply(genesDf, function(df){
      nrow(df)
    }))
    if( sum(nMods) == 0 ) return(NULL)
    # browser()
    # This way os selection matches to one construting the initial read
    # DataFrame
    g <- gff[as.character(S4Vectors::mcols(gff)$ID) %in% genes |
               (as.character(S4Vectors::mcols(gff)$Name) %in% genes &
                  is.na(as.character(S4Vectors::mcols(gff)$ID)))]
    S4Vectors::mcols(gff)[is.na(as.character(S4Vectors::mcols(gff)$ID)),]$ID <- 
      S4Vectors::mcols(gff)[is.na(as.character(S4Vectors::mcols(gff)$ID)),]$Name
    S4Vectors::mcols(g)$ID <- factor(S4Vectors::mcols(g)$ID, levels = genes)
    g <- g[order(S4Vectors::mcols(g)$ID),]
    
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
      df$score <- df$RNAmodR_signal
      return(df)
    })
    genesDf <- genesDf[!vapply(genesDf, is.null, logical(1))]
    if(length(genesDf) == 0) return(NULL)
    res <- do.call(rbind,genesDf)
    res <- res[order(res$chrom,
                     res$start,
                     res$end),]
    return(res)
  })
  
  l <- l[!vapply(l, is.null, logical(1))]
  if(length(l) == 0){
    stop("No modifications detected. Aborting...",
         call. = FALSE)
  }
  df <- do.call(rbind,l)
  df <- df[,c("chrom","start","end","strand","source","type","score","ID",
              "Parent","RNAmodR_type","RNAmodR_signal","RNAmodR_signal_sd",
              "RNAmodR_z","RNAmodR_nbReplicates")]
  
  df$source <- factor(df$source)
  df$type <- factor(df$type)
  return(df)
}

# Construct DataFrame from RNAmod results per gene
.get_dataframe_per_gene <- function(gene){
  df <- S4Vectors::DataFrame(start = vapply(gene,"[[",numeric(1),"location"),
              end = vapply(gene,"[[",numeric(1),"location"),
              ID = names(gene),
              RNAmodR_type = vapply(gene,"[[",character(1),"type"),
              RNAmodR_signal = vapply(gene,"[[",numeric(1),"signal"),
              RNAmodR_signal_sd = vapply(gene,"[[",numeric(1),"signal.sd"),
              RNAmodR_z = vapply(gene,"[[",numeric(1),"z"),
              RNAmodR_nbReplicates = vapply(gene,"[[",numeric(1),"nbsamples"))
  return(df)
}

# conctructs a SummarizedExperiment from modifications and position data
.construct_SE_from_mod_result <- function(name,
                                          files,
                                          gff,
                                          data,
                                          positions){
  # rownames and colnames(modTypes)
  modTypes <- unique(data$RNAmodR_type)
  rownames <- S4Vectors::mcols(gff)$ID
  rownames[is.na(rownames)] <- 
    S4Vectors::mcols(gff[is.na(S4Vectors::mcols(gff)$ID),])$Name
  
  # generate col- and rowdata
  colData <- S4Vectors::DataFrame(name = name, 
                                  files = S4Vectors::SimpleList(files))
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
                                       gene[gene$RNAmodR_type == modType,]
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