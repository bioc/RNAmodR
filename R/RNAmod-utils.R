#' @include RNAmod.R
NULL


#' @rdname setupWorkEnvir
#'
#' @param experimentName a name for the experiment, e.g. "xyz"
#'
#' @return TRUE
#' @export
#'
#' @examples
#' experimentName <- "test"
#' setupWorkEnvir(experimentName)
setMethod(
  f = "setupWorkEnvir", 
  signature = signature(experimentName = "character"),
  definition = function(experimentName = "experiment") {
    assertive::assert_is_scalar(experimentName)
    assertive::assert_is_a_non_missing_nor_empty_string(experimentName)
    
    createFolder = function(folder){
      if (!dir.exists(folder)) {
        dir.create(folder, recursive = TRUE)
      }
    }
    
    dataFolder <- paste0(experimentName, "/data/")
    resultFolder <- paste0(experimentName, "/results/")
    
    createFolder(dataFolder)
    createFolder(resultFolder)
    
    file.copy(from = system.file("extdata", 
                                 "example_experiment_layout.csv",
                                 package = "RNAmod"),
              to = paste0(dataFolder,
                          "experiment_layout.csv") )
    
    message(paste0("work enviroment '",experimentName,"' created."))
    
    return(invisible(TRUE))
  }
)

# reading BAM input ------------------------------------------------------------

#' converts results from scanBam to DataFrame
#' returns a simplified DataFrame for read data read out from a BAM
#' file
#' 
#' @importFrom stats setNames
#' @importClassesFrom S4Vectors DataFrame
.convert_bam_to_DataFrame <- function(dataList, 
                                      param){
  .unlist <- function (x){
    ## do.call(c, ...) coerces factor to integer, which is undesired
    x1 <- x[[1L]]
    if (is.factor(x1)) {
      structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
      do.call(c, x)
    }
  } 
  ids <- .get_IDs_from_scanBamParam(dataList,param)
  dataList <- unname(dataList) # names not valid for unlisting of result
  elts <- stats::setNames(Rsamtools::bamWhat(param), 
                          Rsamtools::bamWhat(param))
  
  lst <- lapply(elts, function(elt) .unlist(lapply(dataList, "[[", elt)))
  lst$ID <- ids
  lst <- lst[c("ID"="ID",elts)]
  
  df <- do.call(S4Vectors::DataFrame, lst)
  return(df)
}

.get_IDs_from_scanBamParam <- function(dataList,param){
  geneNames <- vector(mode = "list", 
                      length=length(names(Rsamtools::bamWhich(param))))
  # construct start-end part
  # format = start-end
  geneNames <- lapply(Rsamtools::bamWhich(param), function(x){
    pos <- paste0(IRanges::start(x),"-",IRanges::end(x))
    # If ID column is NA use Name column
    res <- S4Vectors::mcols(x)$ID
    res[is.na(res)] <- S4Vectors::mcols(x[is.na(S4Vectors::mcols(x)$ID),])$Name
    
    names(res) <- paste0(IRanges::start(x),"-",IRanges::end(x))
    return(res)
  })
  # combine chromosome with start-end part
  # format = chr:start-end
  names(geneNames) <- names(Rsamtools::bamWhich(param))
  for(i  in seq_along(geneNames)){
    names(geneNames[[i]]) <- paste0(names(geneNames[i]),
                                    ":",
                                    names(geneNames[[i]]))
  }
  geneNames <- setNames(unlist(geneNames),unlist(lapply(geneNames, names)))
  
  # Repeat every "chr:start-end" = "ID" pair for number of reads mapped to this
  # part
  geneNames2 <- vector(mode="list",length(geneNames))
  for(i in seq_along(geneNames)){
    geneNames2[[i]] <- list(rep(geneNames[i],
                                length(dataList[[names(geneNames[i])]]$seq)))
  }
  geneNames <- setNames(unlist(geneNames2),unlist(lapply(geneNames2, names)))
  return(geneNames)
}

#' scans BAM files and converts them to DataFrame each 
#' 
#' @import Biostrings
.get_bam_data <- function(files, 
                          param){
  
  res <- c()
  files <- files[grepl(".bam",files)]
  if( length(files) > 0){
    
    FUN <-  function(x, param){
      requireNamespace("Rsamtools", quietly = TRUE)
      
      y <- scanBam(x, param=param)
      y <- .convert_bam_to_DataFrame(y, param)
      return(y)
    }
    
    res <- bplapply(files,
                    FUN,
                    param = param)
    names(res) <- files
  }
  return(res)
  
}

# Returns the number of reads contained in a bam file for parameters set in
# param
.get_bam_read_count <- function(files, quality){
  
  res <- c()
  files <- files[grepl(".bam",files)]
  if( length(files) > 0){
    for(i in 1:length(files)) {
      res[[i]] <- countBam(files[[i]], 
                           param = ScanBamParam(mapqFilter = quality))$records
    }
  }
  return(res)
}


#' return parameters to be used by scanBam
#' gRangeInput a GRanges object containing the ranges for search in the BAM file
#' quality quality argument used for scanBamParam
#' 
#' @importFrom GenomeInfoDb seqnames
#' @importFrom Rsamtools ScanBamParam
.assemble_scanBamParam <- function(gRangeInput,
                                   quality,
                                   acceptableChromIdent){
  # ScanBamParam expects GRangesList each members matching a chromosome
  gRangeList <- GenomicRanges::GRangesList()
  listNames <- rep("",length(unique(GenomeInfoDb::seqnames(gRangeInput))))
  for(i in seq_along(unique(GenomeInfoDb::seqnames(gRangeInput)))){
    ident <- as.character(unique(GenomeInfoDb::seqnames(gRangeInput))[i])
    if( ident %in% acceptableChromIdent ){
      gRangeList <- append(gRangeList, GenomicRanges::GRangesList(gRangeInput[GenomeInfoDb::seqnames(gRangeInput) == ident,]) )
      listNames[[i]] <- ident
    } else {
      warning("Not matching chromosome identifier in gff and bam file. ",
              "Skipping data for chromosome '",ident,"'",
              call. = FALSE)
    }
  }
  names(gRangeList) <- listNames[listNames != ""]
  
  which <- gRangeList
  what <- c("rname", "strand", "pos", "qwidth", "seq", "mapq")
  param <- Rsamtools::ScanBamParam(which=which, 
                                   what=what
                                   # since tRNA and
                                   # mapqFilter = quality
                                   )
  return(param)
}

# Extracts sequence names aka. chromosome identifier from list of bam files
.get_acceptable_chrom_ident <- function(bamFiles){
  seqnames <- lapply(bamFiles, function(file){
    res <- Rsamtools::idxstatsBam(file)
    return(as.character(res$seqnames))
  })
  return(unique(unlist(seqnames)))
}

# Check if name exists in any parent env apart from global env
.where <- function(name, n = 1) {
  env = parent.frame( n = n)
  while(!identical(env, sys.frame(which = 0))) {
    if (exists(name, envir = env, inherits = FALSE)) {
      # success case
      return(env)
    }
    # inspect parent
    n <- n + 1
    env <- parent.frame(n = n)
  }
  stop("Can't find ", name, call. = FALSE)
}

# subsets a gRanges object for entries concerning a single ID
# - Name/ID for parent entry
# - al childs based on Parent name equals ID
.subset_gff_for_unique_transcript <- function(gff, ID){
  res <- gff[(is.na(S4Vectors::mcols(gff)$ID) & 
                S4Vectors::mcols(gff)$Name == ID) |
               (!is.na(S4Vectors::mcols(gff)$ID) & 
                  S4Vectors::mcols(gff)$ID == ID) |
               (!is.na(as.character(S4Vectors::mcols(gff)$Parent)) & 
                  as.character(S4Vectors::mcols(gff)$Parent) == ID),]
  if(length(res) != 1){
    stop("Error gff annotation. Non unique ID detected: ",
         ID,
         .call = FALSE)
  }
  return(res)
}



#' @rdname matchFastaToGff
#'
#' @param inputFasta fasta file to be used for renaming
#' @param inputGFF gff file containing gene, mRNA and CDS annotation data
#'
#' @return TRUE
#' @export
#' 
#' @import Biostrings
#' @import rtracklayer
setMethod(
  f = "matchFastaToGff", 
  signature = signature(inputFasta = "character", 
                        inputGFF = "character"),
  definition = function(inputFasta , 
                        inputGFF){
    requireNamespace("Biostrings", quietly = TRUE)
    requireNamespace("rtracklayer", quietly = TRUE)
    assertive::assert_all_are_existing_files(c(inputFasta,inputGFF))
    
    fsa <- readDNAStringSet(inputFasta)
    gff <- import.gff3(inputGFF)
    
    # number of chromosomes and sequences have to match
    if( !assertive::are_same_length(names(fsa), unique(rtracklayer::chrom(gff)))){
      stop("Fasta and GFF file don't have matching number of sequences/",
           "sequence identifiers. Please make sure that the number of ",
           "fasta sequences matches the number of unique chromosome ",
           "identifiers in the GFF file. The names of the fasta ",
           "sequences will be overridden with the chromosomal ",
           "identifier from the GFF file.",
           call. = FALSE)
    }
    
    names(fsa) <- unique(chrom(gff))
    
    fileName <- gsub(".fsa","_fixed.fsa",inputFasta)
    writeXStringSet(fsa, fileName)
    
    message(paste0("Fasta file ", fileName, " written."))
    
    return(invisible(TRUE))
  }
)