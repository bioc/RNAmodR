#' @include RNAmodR.R
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
                                 package = "RNAmodR"),
              to = paste0(dataFolder,
                          "experiment_layout.csv") )
    
    message(paste0("work enviroment '",experimentName,"' created."))
    
    return(invisible(TRUE))
  }
)

# reading BAM input ------------------------------------------------------------

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
  # what <- c("rname", "strand", "pos", "qwidth", "seq", "mapq")
  what <- c("seq", "mapq")
  flags <- scanBamFlag(isSecondaryAlignment = FALSE)
  param <- Rsamtools::ScanBamParam(flag = flags,
                                   which = which, 
                                   what = what,
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

# # Check if name exists in any parent env apart from global env
# .where <- function(name, n = 1) {
#   env = parent.frame( n = n)
#   while(!identical(env, sys.frame(which = 0))) {
#     if (exists(name, envir = env, inherits = FALSE)) {
#       # success case
#       return(env)
#     }
#     # inspect parent
#     n <- n + 1
#     env <- parent.frame(n = n)
#   }
#   stop("Can't find ", name, call. = FALSE)
# }

# BiocGeneric helper functions -------------------------------------------------

# strand related functions
.get_strand <- function(x){
  as.character(BiocGenerics::strand(x))
}
.get_unique_strand <- function(x){
  unique(.get_strand(x))
}
.is_minus_strand <- function(x) {
  all(as.logical(.get_strand(x) == "-"))
}
.is_on_minus_strand <- function(x) {
  all(.is_on_correct_strand(x,"-"))
}
# is on the minus strand?
.is_on_correct_strand <- function(x, strand) {
  as.logical(.get_strand(x) == strand) 
}
.is_on_correct_strand2 <- function(x, gr) {
  as.logical(.get_strand(x) == .get_unique_strand(gr)) 
}

# GRanges helper functions -----------------------------------------------------

# get aggregated identifiers
.get_gr_ids <- function(gr,
                        na.rm = TRUE){
  IDs <- unique(c(gr$ID,gr$Name))
  if(na.rm){
    IDs <- IDs[!is.na(IDs)]
  }
  IDs
}

# get unique ids
# returns all identifiers which can be used as a identifier in the Parent
# column aggrgated by precident: ID > Name > gene
.get_unique_identifiers <- function(gr){
  ids <- as.character(gr$ID)
  ids[is.na(ids)] <- as.character(gr[is.na(ids)]$Name)
  ids[is.na(ids)] <- as.character(gr[is.na(ids)]$gene)
  ids
}

# subset to types present in RNAMODR_MOD_SEQ_FEATURES
.subset_rnamod_transcript_features <- function(gr){
  gr[S4Vectors::mcols(gr)$type %in% RNAMODR_MOD_TRANSCRIPT_FEATURES,]
}
# subset to types present in RNAMODR_MOD_SEQ_FEATURES
.subset_rnamod_containing_features <- function(gr){
  gr[S4Vectors::mcols(gr)$type %in% RNAMODR_MOD_CONTAINING_FEATURES,]
}

# subsets a gRanges object for entries concerning a single ID
# - Name/ID for parent entry
# - all childs based on Parent name equals ID
.subset_gff_for_unique_transcript <- function(gff, 
                                              ID,
                                              wo.childs = TRUE){
  gr <- gff[(is.na(S4Vectors::mcols(gff)$ID) & 
                S4Vectors::mcols(gff)$Name == ID) |
               (!is.na(S4Vectors::mcols(gff)$ID) & 
                  S4Vectors::mcols(gff)$ID == ID),]
  if(!wo.childs){
    gr <- append(gr, .get_childs(gff,ID))
  }
  .order_GRanges(gr)
}
# recursive search for all childs of an ID
.get_childs <- function(gff, ID){
  gr <- gff[!is.na(as.character(gff$Parent)) & 
              as.character(gff$Parent) == ID,]
  if(length(gr) == 0) return(GRanges())
  # Use ID and Name for child search
  IDs <- .get_gr_ids(gr)
  IDs <- IDs[IDs != ID]
  # Walk up the parent chain
  grl <- append(list(gr),lapply(IDs, function(x){.get_childs(gff,x)}))
  grl <- grl[!vapply(grl,is.null,logical(1))]
  return(unlist(GRangesList(grl)))
}

# subset GRanges for highest ranking parent
.get_parent_annotations <- function(gr,
                                    forceSingle = FALSE,
                                    doRecursiveSearch = FALSE,
                                    IDs){
  if(!doRecursiveSearch){
    res <- gr[is.na(as.character(gr$Parent)),]
    if(length(res) > 1 && forceSingle) 
      stop("No distinct single parent annotation detectable.")
    return(res)
  }
  if(missing(IDs)) IDs <- .get_gr_ids(gr)
  l <- lapply(IDs, .get_parent, gr = gr)
  l <- unique(unlist(l[!vapply(l,is.null,logical(1))]))
  .order_GRanges(gr[(!is.na(gr$ID) & gr$ID %in% l) |
                      (!is.na(gr$Name) & gr$Name %in% l),])
}
# quasi recursive search for all possible parents
.get_parent <- function(ID,
                        gr){
  parent <- as.character(gr[(!is.na(gr$ID) & gr$ID == ID) |
                              (!is.na(gr$Name) & gr$Name == ID),]$Parent)
  parent <- parent[!is.na(parent)]
  if(length(parent) == 0) return(ID)
  return(unlist(lapply(parent, .get_parent, gr)))
}

# returns an order based on the named of the chromosome and start and end 
# coordinates
.order_GRanges <- function(gr){
  gr[order(GenomeInfoDb::seqnames(gr),
           BiocGenerics::start(gr),
           BiocGenerics::end(gr))]
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
