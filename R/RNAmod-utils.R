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
  ids <- .get_IDs_from_scanBamParam_all_reads(param,
                                              dataList)
  dataList <- unname(dataList) # names not valid for unlisting of result
  elts <- stats::setNames(Rsamtools::bamWhat(param), 
                          Rsamtools::bamWhat(param))
  
  lst <- lapply(elts, function(elt) .unlist(lapply(dataList, "[[", elt)))
  lst$ID <- ids
  lst <- lst[c("ID"="ID",elts)]
  
  df <- do.call(S4Vectors::DataFrame, lst)
  return(df)
}

# repeats names in the format chromosome:start-end according to the number of 
# reads for each in element bamWhich(param)
.get_IDs_from_scanBamParam_all_reads <- function(param,
                                                 dataList){
  geneNames <- .get_IDs_from_scanBamParam(param)
  # Repeat every "chr:start-end" = "ID" pair for number of reads mapped to this
  # part
  geneNames2 <- lapply(seq_along(geneNames), function(i){
    rep(geneNames[i],
             length(dataList[[names(geneNames[i])]]$seq))
  })
  geneNames <- setNames(unlist(geneNames2),
                        unlist(lapply(geneNames2, names)))
  return(geneNames)
}

# assembles a name in the format of chromosome:start-end for each element in 
# bamWhich(param)
.get_IDs_from_scanBamParam <- function(param){
  # construct start-end part
  # format = start-end
  geneNamesPerChrom <- lapply(Rsamtools::bamWhich(param), function(x){
    # If ID column is NA use Name column
    res <- S4Vectors::mcols(x)$ID
    res[is.na(res)] <- S4Vectors::mcols(x[is.na(S4Vectors::mcols(x)$ID),])$Name
    names(res) <- paste0(IRanges::start(x),"-",IRanges::end(x))
    return(res)
  })
  # combine chromosome with start-end part
  # format = chr:start-end
  names(geneNamesPerChrom) <- names(Rsamtools::bamWhich(param))
  geneNamesPerChrom <- lapply(seq_along(geneNamesPerChrom), function(i){
    setNames(geneNamesPerChrom[[i]],
             paste0(names(geneNamesPerChrom[i]),
                    ":",
                    names(geneNamesPerChrom[[i]])))
  })
  names(geneNamesPerChrom) <- names(Rsamtools::bamWhich(param))
  geneNames <- setNames(unlist(geneNamesPerChrom),
                        unlist(lapply(geneNamesPerChrom, names)))
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
    
    res <- BiocParallel::bplapply(files,
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
    res 
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
      gRangeList <- append(gRangeList, 
                           GenomicRanges::GRangesList(
                             gRangeInput[GenomeInfoDb::seqnames(gRangeInput) 
                                         == ident,]) )
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
  flags <- scanBamFlag(isSecondaryAlignment = FALSE)
  param <- Rsamtools::ScanBamParam(flag = flags,
                                   which=which, 
                                   what=what,
                                   mapqFilter = .get_map_quality())
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

# GRanges helper functions -----------------------------------------------------

# subsets a gRanges object for entries concerning a single ID
# - Name/ID for parent entry
# - all childs based on Parent name equals ID
.subset_gff_for_unique_transcript <- function(gff, ID){
  browser()
  res <- gff[(is.na(S4Vectors::mcols(gff)$ID) & 
                S4Vectors::mcols(gff)$Name == ID) |
               (!is.na(S4Vectors::mcols(gff)$ID) & 
                  S4Vectors::mcols(gff)$ID == ID) |
               (!is.na(as.character(S4Vectors::mcols(gff)$Parent)) & 
                  as.character(S4Vectors::mcols(gff)$Parent) == ID &
                  S4Vectors::mcols(gff)$type %in% RNAMOD_MOD_SEQ_FEATURES),]
  .order_GRanges(res)
}

# returns an order based on the named of the chromosome and start and end 
# coordinates
.order_GRanges <- function(gr){
  gr[order(GenomeInfoDb::seqnames(gr),
           BiocGenerics::start(gr),
           BiocGenerics::end(gr))]
}

# returns all identifiers whcih can be used as a identifier in the Parent
# column
.get_parent_identifiers <- function(gr){
  parents <- as.character(gr$ID)
  parents[is.na(parents)] <- as.character(gr[is.na(parents)]$Name)
  parents[is.na(parents)] <- as.character(gr[is.na(parents)]$gene)
  parents
}

# sequence handling ------------------------------------------------------------

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
    if( !assertive::are_same_length(names(fsa), 
                                    unique(rtracklayer::chrom(gff)))){
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

# retrieve concats DNAStrings based on the strand
.get_seq_for_unique_transcript <- function(gff,fafile,ID){
  if(length(unique(as.character(BiocGenerics::strand(gff)))) != 1) {
    stop("Ambigeous type of GRanges given. Expects only one strand.",
         call. = FALSE)
  }
  seq <- getSeq(fafile,gff)
  if(as.character(BiocGenerics::strand(gff)) == "-"){
    return(unlist(rev(seq)))
  }
  return(unlist(seq))
}


# debug helper -----------------------------------------------------------------

.print_location_info <- function(location, locs){
  message("location: ",location," //locs: ",paste(locs,collapse = ","))
}
.print_transcript_info <- function(ID, iterationN){
  message(ID, " - iteration: ",iterationN)
}
