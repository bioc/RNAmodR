#' @include RNAmod.R
NULL

RNAMOD_SEP_SEQNAMES <- "---"

#' @name convertGeneToChrom
#' 
#' @aliases convertGeneToChrom convertGeneTypeToChrom
#' 
#' @title convert gene sequences to individual chromosomes
#'
#' @param gfffile 
#' @param fafile 
#' @param geneTypes 
#' @param ident 
#' @param saveTrimmedFile 
#' @param appendToOriginal 
#'
#' @return
#' @export
#' 
#' @importFrom S4Vectors mcols Rle
#' @importFrom rtracklayer import.gff3 export.gff3 chrom
#' @importFrom Rsamtools FaFile indexFa getSeq
#' @importFrom GenomeInfoDb seqnames
#' @importFrom stringr str_split
#' @importFrom BiocGenerics which unlist width start end strand
#' @import Biostrings
#' @import tRNAscan2GRanges
#'
#' @examples
#' \donttest{
#' gfffile <-
#' fafile <-
#' convertGeneTypeToChrom(gfffile,
#'                        fafile,
#'                        "tRNA_gene",
#'                        "tRNA")
#' convertGeneToChrom(gfffile,
#'                    fafile,
#'                    "RDN18-1",
#'                    "rRNA")
#' }
convertGeneTypeToChrom <- function(gfffile,
                                   fafile,
                                   geneTypes,
                                   ident,
                                   saveTrimmedFile = FALSE,
                                   appendToOriginal = FALSE){
  # Check inputs
  assertive::assert_all_are_existing_files(gfffile)
  assertive::assert_all_are_existing_files(fafile)
  assertive::assert_is_a_non_empty_string(ident)
  assertive::assert_is_a_bool(saveTrimmedFile)
  assertive::assert_is_a_bool(appendToOriginal)
  # Load files
  gff <- rtracklayer::import.gff3(gfffile)
  # check that gene types are in types column of gff file
  checkType <- unique(S4Vectors::mcols(gff)$type)
  assertive::assert_is_subset(geneTypes, checkType)
  # get individual ID of genes with gene type
  gff_sub <- gff[S4Vectors::mcols(gff)$type %in% geneTypes,]
  IDs <- S4Vectors::mcols(gff_sub)$ID
  # hand of the individual IDs
  convertGeneToChrom(gfffile,fafile,IDs,ident,appendToOriginal,saveTrimmedFile)
}

#' @rdname convertGeneToChrom
#'
#' @param genes 
#' 
#' @export
convertGeneToChrom <- function(gfffile,
                               fafile,
                               genes,
                               ident,
                               saveTrimmedFile = FALSE,
                               appendToOriginal = FALSE){
  # Check inputs
  assertive::assert_all_are_existing_files(gfffile)
  assertive::assert_all_are_existing_files(fafile)
  assertive::assert_is_a_non_empty_string(ident)
  assertive::assert_is_a_bool(saveTrimmedFile)
  assertive::assert_is_a_bool(appendToOriginal)
  # get the aggregated file names
  fileNames <- .get_file_names(gfffile,
                               fafile,
                               ident)
  # Load files
  gff <- rtracklayer::import.gff3(gfffile)
  fa <- Rsamtools::FaFile(fafile)
  # check that genes are in ID or Name column of gff file
  checkGenes <- c(unique(S4Vectors::mcols(gff)$ID),
                  unique(S4Vectors::mcols(gff)$Name))
  checkGenes <- unique(checkGenes)
  assertive::assert_is_non_empty(genes)
  assertive::assert_is_subset(genes, checkGenes)
  # get the matchin gff entries
  gff_subset <- .get_gff_subset(gff, genes)
  # Remove all entries from the original gff, which are part of the subset
  gff_trimmed <- .get_trimmed_gff(gff, 
                                  gff_subset,
                                  saveTrimmedFile,
                                  fileNames[["trim_gfffile"]])
  # extract the coding sequences of all subset features
  seqs_subset <- .extract_exon_sequences(gff_subset,fa)
  # remove indistinguishable sequences (applies usually to tRNA, rRNA and such)
  # apply by using the columns id-name-gene for the search
  gff_subset_clustered <- .cluster_results_gff(gff_subset,seqs_subset)
  # cluster by sequence
  seqs_subset_clustered <- .cluster_results_seqs(seqs_subset)
  # match gff and fasta sequences and create a new gff annotation and fasta 
  # sequences
  gff_subset_clustered <- .update_gff_subset(gff_subset_clustered,
                                             seqs_subset_clustered)
  seqs_subset_clustered <- .update_seqs(seqs_subset_clustered)
  # write general output to files
  .write_output(gff_subset_clustered,
                seqs_subset_clustered,
                gff_subset,
                fileNames)
  # if appendToOriginal == TRUE saved appended version of the original files
  # with the clustered subset files
  if(saveTrimmedFile && appendToOriginal){
    fileNames[["out_fafile"]] <- 
      appendFastaFiles(fileNames[["in_fafile"]],
                       fileNames[["convert_fafile"]],
                       ident = ident,
                       msg = "Saved appended fasta file: ")
    fileNames[["out_gfffile"]] <- 
      appendGFF(fileNames[["trim_gfffile"]],
                fileNames[["convert_gfffile"]],
                ident = ident,
                msg = "Saved appended gff3 file: ")
  }
  # return file names
  return(invisible(fileNames))
}

# get a subset off gff annotation from ID, name or parent
.get_gff_subset <-  function(gff, genes){
  # create gff subset
  gff_sub <- gff[(is.na(S4Vectors::mcols(gff)$ID) & 
                    S4Vectors::mcols(gff)$Name %in% genes) |
                   (!is.na(S4Vectors::mcols(gff)$ID) & 
                      S4Vectors::mcols(gff)$ID %in% genes) |
                   (!is.na(S4Vectors::mcols(gff)$gene) & 
                      S4Vectors::mcols(gff)$gene %in% genes) |
                   (!is.na(as.character(S4Vectors::mcols(gff)$Parent)) & 
                      as.character(S4Vectors::mcols(gff)$Parent) %in% genes) |
                   (is.na(S4Vectors::mcols(gff)$ID) & 
                      S4Vectors::mcols(gff)$Name %in% paste0(genes,"_mRNA")) |
                   (!is.na(S4Vectors::mcols(gff)$ID) & 
                      S4Vectors::mcols(gff)$ID %in% paste0(genes,"_mRNA")) |
                   (!is.na(as.character(S4Vectors::mcols(gff)$Parent)) & 
                      as.character(S4Vectors::mcols(gff)$Parent) %in% 
                      paste0(genes,"_mRNA")),]
  if(length(gff_sub) == 0){
    stop("No genes found in supplied gff file with the given genes names.",
         call. = FALSE)
  }
  gff_sub <- gff_sub[order(rtracklayer::chrom(gff_sub), 
                           BiocGenerics::start(gff_sub))]
  gff_sub
}
# get the trimmed gff annotation
.get_trimmed_gff <- function(gff, 
                             gff_subset,
                             saveTrimmedFile,
                             trim_gfffile){
  gff_trimmed <- IRanges::subsetByOverlaps(gff, 
                                           gff_subset,
                                           type = "equal")
  gff_trimmed <- IRanges::subsetByOverlaps(gff, 
                                           gff_trimmed,
                                           type = "equal",
                                           invert = TRUE)
  if(saveTrimmedFile){
    # save the reduced gff annotation file, which can be used for further steps
    .write_gff(gff_trimmed, 
               trim_gfffile,
               "Saved trimmed gff3 file: ")
  }
  gff_trimmed
}


#' @rdname convertGeneToChrom
#'
#' @param genes 
#' 
#' @export
converttRNAscanToChrom <- function(gfffile,
                                   fafile,
                                   tRNAscanfile,
                                   ident,
                                   saveTrimmedFile = FALSE,
                                   appendToOriginal = FALSE){
  # Check inputs
  assertive::assert_all_are_existing_files(gfffile)
  assertive::assert_all_are_existing_files(fafile)
  assertive::assert_all_are_existing_files(tRNAscanfile)
  assertive::assert_is_a_non_empty_string(ident)
  assertive::assert_is_a_bool(saveTrimmedFile)
  assertive::assert_is_a_bool(appendToOriginal)
  # get the aggregated file names
  fileNames <- .get_file_names(gfffile,
                               fafile,
                               ident)
  # Load files
  gff <- rtracklayer::import.gff3(gfffile)
  # load tRNAscan data
  tRNAscan <- tRNAscan2GRanges::tRNAscan2GFF(tRNAscanfile)
  # generate a trimmed version and save if necessary
  gff_subset <- .get_trimmed_tRNAscan(tRNAscan,
                                      gff,
                                      saveTrimmedFile,
                                      fileNames[["trim_gfffile"]])
  # unknown tRNA should also be masked. therefore exchange gff_sub for tRNAscan
  gff_subset <- tRNAscan[,colnames(S4Vectors::mcols(tRNAscan)) %in% 
                           c("source","type","score","phase","ID","Name","gene")]
  # extract the coding sequences of tRNA
  seqs_subset <- DNAStringSet(tRNAscan$tRNA_seq)
  # Add CCA to end, if the CCA end is not encoded
  CCAseq <- DNAStringSet(lapply(S4Vectors::mcols(tRNAscan)$tRNA_CCA.end, 
                                function(bool){
    if(!bool) return(DNAString("CCA"))
    DNAString("")
  }))
  seqs_subset <- DNAStringSet(lapply(seq_len(length(seqs_subset)), function(i){
    d <- DNAStringSet(list(seqs_subset[[i]],CCAseq[[i]]))
    unlist(d)
  }))
  # set new names which become chromosome names
  names(seqs_subset) <- tRNAscan$ID
  # remove indistinguishable sequences from tRNAscan
  # get clustered subset
  tRNA_subset_clustered <- tRNAscan[,!(colnames(S4Vectors::mcols(tRNAscan)) %in%
                                         c("tRNA_seq",
                                           "tRNA_str"))]
  seqs_subset_clustered <- .cluster_results_seqs(seqs_subset)
  # create new gff annotation
  tRNA_subset_clustered <- 
    tRNA_subset_clustered[which(tRNA_subset_clustered$ID 
                                %in% names(seqs_subset_clustered)),]
  tRNA_subset_clustered <- .reset_coord(tRNA_subset_clustered, 
                                        seqs_subset_clustered)
  if( is.null(tRNA_subset_clustered$Name)) tRNA_subset_clustered$Name <- NA
  if( is.null(tRNA_subset_clustered$gene)) tRNA_subset_clustered$gene <- NA
  tRNA_subset_clustered <-
    tRNA_subset_clustered[, colnames(S4Vectors::mcols(tRNA_subset_clustered))
                          %in% c("source", "type", "score",
                                 "phase", "ID", "Name", "gene")]
  tRNA_subset_clustered <- GenomicRanges::GRanges(
    S4Vectors::Rle(tRNA_subset_clustered$ID),
    ranges = IRanges::ranges(tRNA_subset_clustered),
    strand = as.character(BiocGenerics::strand(tRNA_subset_clustered)),
    S4Vectors::mcols(tRNA_subset_clustered))
  # write output
  .write_output(tRNA_subset_clustered,
                seqs_subset_clustered,
                gff_subset,
                fileNames)
  if(saveTrimmedFile && appendToOriginal){
    fileNames[["out_fafile"]] <- appendFastaFiles(fileNames[["in_fafile"]],
                                          fileNames[["convert_fafile"]],
                                          ident = ident,
                                          msg = "Saved appended fasta file: ")
    fileNames[["out_gfffile"]] <- appendGFF(fileNames[["trim_gfffile"]],
                                          fileNames[["convert_gfffile"]],
                                          ident = ident,
                                          msg = "Saved appended gff3 file: ")
  }
  return(invisible(fileNames))
}

.get_trimmed_tRNAscan <- function(tRNAscan,gff,saveTrimmedFile,trim_gfffile){
  # test if chromosome names of tRNAscan need to be fixed
  chromTest <- intersect(unique(as.character(GenomeInfoDb::seqnames(tRNAscan))),
                         unique(as.character(GenomeInfoDb::seqnames(gff))))
  if(!identical(unique(as.character(GenomeInfoDb::seqnames(tRNAscan))), 
                chromTest)){
    # if yes
    # construct a translation library
    chrN <- length(unique(as.character(GenomeInfoDb::seqnames(tRNAscan))))
    chrom_names <- unique(as.character(GenomeInfoDb::seqnames(gff)))[seq_len(chrN)]
    names(chrom_names) <- unique(as.character(GenomeInfoDb::seqnames(tRNAscan)))
    # get new chromosome names
    newChromNames <- unlist(lapply(GenomeInfoDb::seqnames(tRNAscan), 
                                   function(x){
                                     chrom_names[names(chrom_names) == x]
                                   }))
    # create updated tRNAscan object for subsetting
    tRNAscan_new <- GenomicRanges::GRanges(
      S4Vectors::Rle(newChromNames),
      ranges = IRanges::ranges(tRNAscan),
      strand = as.character(BiocGenerics::strand(tRNAscan)),
      S4Vectors::mcols(tRNAscan))
  } else {
    tRNAscan_new <- tRNAscan
  }
  # Because tRNA His does not exactly overlap, enlarge the boundaries each by
  # 3 nt
  start(tRNAscan_new) <- start(tRNAscan_new) - 3
  end(tRNAscan_new) <- end(tRNAscan_new) + 3
  # Remove all entries from the original gff, which are part of the tRNAscan
  gff_subset <- suppressWarnings(IRanges::subsetByOverlaps(gff, 
                                                           tRNAscan_new,
                                                           type = "within"))
  if(length(gff_subset) == 0){
    stop("No overlap found in supplied gff and tRNAscan file.",
         call. = FALSE)
  }
  gff_subset <- gff_subset[order(rtracklayer::chrom(gff_subset), 
                                 BiocGenerics::start(gff_subset))]
  # Remove all entries from the original gff, which are part of the subset
  gff_trimmed <- .get_trimmed_gff(gff, 
                                  gff_subset,
                                  saveTrimmedFile,
                                  trim_gfffile)
  return(gff_subset)
}

# identifing annotations and manipulations identifiers -------------------------
# merge seq_identifier
.get_unique_seqnames <- function(gff){
  paste0(S4Vectors::mcols(gff)$ID,
         RNAMOD_SEP_SEQNAMES,
         S4Vectors::mcols(gff)$Name,
         RNAMOD_SEP_SEQNAMES,
         S4Vectors::mcols(gff)$gene)
}

# split seq_identifier
.reverse_unique_seqnames <- function(names){
  res <- str_split(names, RNAMOD_SEP_SEQNAMES)
  ID <- vapply(res,function(x){x[1]},character(1))
  Name <- vapply(res,function(x){x[2]},character(1))
  gene <- vapply(res,function(x){x[3]},character(1))
  ID[ID == "NA"] <- NA
  Name[Name == "NA"] <- NA
  gene[gene == "NA"] <- NA
  return(list(ID = ID,
              Name = Name,
              gene = gene))
}

# remove special character from chromosome names
.condense_chrom_names <- function(names){
  # remove unnecessary additions
  res <- lapply(names, function(x){
    x <- gsub("_CDS","",x)
    x <- gsub("_noncoding_exon","",x)
    x[is.na(x)] <- ""
    x
  })
  names(res) <- names(names)
  # remove duplicated naming elements
  for(i in seq_along(res[[1]])){
    x <- unique(c(res$ID[[i]],res$Name[[i]],res$gene[[i]]))
    x <- x[x != ""]
    res$ID[[i]] <- x[1]
    res$Name[[i]] <- x[2]
    res$gene[[i]] <- x[3]
  }
  res <- lapply(res, function(x){
    x[is.na(x)] <- ""
    x
  })
  # condense naming elements intro string
  res <- paste0(res$ID,
                "_",
                res$Name,
                "_",
                res$gene)
  res <- gsub("_$","",res)
  res <- gsub("^_","",res)
  res <- gsub("_$","",res)
  res <- gsub("^_","",res)
  # escape unsupported characters in chromnames
  return(.fix_chrom_names(res))
}
# escape unsupported characters in chromnames
.fix_chrom_names <- function(x){
  gsub("(?![a-zA-Z0-9.:^*$@!+_?-|]).","-",
       x, 
       perl = TRUE)
}

# getting sequences ------------------------------------------------------------
# extract the exons sequences and merges them into one per gene
.extract_exon_sequences <- function(gff_sub, fa){
  # subset to exons only
  exons <- gff_sub[S4Vectors::mcols(gff_sub)$type %in% c("CDS",
                                                         "noncoding_exon",
                                                         "exon")]
  # get sequences for exons
  exons_seq <- Biostrings::getSeq(fa,exons)
  names(exons_seq) <- .get_unique_seqnames(exons)
  # get unique names and iterate and merge split sequences
  uniq_gene_names <- unique(names(exons_seq))
  seqs <- Biostrings::DNAStringSet(lapply(uniq_gene_names, function(name){
    seq <- exons_seq[BiocGenerics::which(names(exons_seq) == name)]
    strand <- unique(as.character(BiocGenerics::strand(exons[.get_unique_seqnames(exons) == name])))
    # If strand == "-" intron sequences are in reverse order
    if(strand == "-"){
      return(BiocGenerics::unlist(rev(seq)))
    }
    return(BiocGenerics::unlist(seq))
  }))
  # rename sequences to match parent sequence
  names(seqs) <- uniq_gene_names
  # names(seqs) <- gsub("_CDS","",names(seqs))
  # names(seqs) <- gsub("_noncoding_exon","",names(seqs))
  return(seqs)
}

# clustering of annotations and sequneces --------------------------------------
# removed duplicated sequences
.cluster_results_gff <- function(gff_sub,seqs){
  # get duplicated sequences must match reverse of first function call in
  # .combine_gff_subset_with_seqs
  duplicated_seq_names <- names(seqs[duplicated(seqs)])
  duplicated_seq_names <- append(duplicated_seq_names,
                                 gsub("_CDS","",duplicated_seq_names))
  duplicated_seq_names <- append(duplicated_seq_names,
                                 gsub("_noncoding_exon","",duplicated_seq_names))
  duplicated_seq_names <- append(duplicated_seq_names,
                                 gsub("_exon","",duplicated_seq_names))
  duplicated_seq_names <- unique(duplicated_seq_names)
  gff_sub <- gff_sub[!(.get_unique_seqnames(gff_sub) %in% duplicated_seq_names)]
  return(gff_sub)
}
.cluster_results_seqs <- function(seqs){
  return(seqs[!duplicated(seqs)])
}

# matching annotations and sequences -------------------------------------------
# create new gff annotations
.update_gff_subset <- function(gff_sub, seqs){
  # subset gff to parents and childs
  gff_children <- gff_sub[.get_unique_seqnames(gff_sub) %in% names(seqs)]
  gff_children <- gff_children[!duplicated(.get_unique_seqnames(gff_children))]
  gff_parents <- gff_sub[(!is.na(as.character(S4Vectors::mcols(gff_sub)$ID)) &
                         as.character(S4Vectors::mcols(gff_sub)$ID) %in% 
                          as.character(S4Vectors::mcols(gff_children)$Parent)) |
                         (!is.na(as.character(S4Vectors::mcols(gff_sub)$Name)) &
                         as.character(S4Vectors::mcols(gff_sub)$Name) %in% 
                           as.character(S4Vectors::mcols(gff_children)$Parent))]
  if( length(gff_children) != length(gff_parents) |
      length(seqs) != length(gff_children)){
    warnings("Mismatching number of annotations and sequences",
             call. = FALSE)
  }
  # setup coordinates for genes based on seqs as chromosomes
  gff_parents <- .reset_coord(gff_parents, seqs)
  # chromosome names
  chrom_names <- .reverse_unique_seqnames(.get_unique_seqnames(gff_children))
  chrom_names <- .condense_chrom_names(chrom_names)
  # update/create and combine GRanges
  gff_children <- GenomicRanges::GRanges(S4Vectors::Rle(chrom_names),
                                        ranges = IRanges::ranges(gff_parents),
                                        strand = as.character(BiocGenerics::strand(gff_parents)),
                                        S4Vectors::mcols(gff_children))
  gff_parents <- GenomicRanges::GRanges(S4Vectors::Rle(chrom_names),
                                        ranges = IRanges::ranges(gff_parents),
                                        strand = as.character(BiocGenerics::strand(gff_parents)),
                                        S4Vectors::mcols(gff_parents))
  gff_c <- c(gff_parents,gff_children)
  gff_c <- sort(gff_c)
  return(gff_c)
}
# create new fasta file
.update_seqs <- function(seqs){
  # remove brackets from chrom names
  chrom_names <- .reverse_unique_seqnames(names(seqs))
  chrom_names <- .condense_chrom_names(chrom_names)
  names(seqs) <- chrom_names
  return(seqs)
}
# reset coordinates for new chromosomes
.reset_coord <- function(gff,seqs){
  BiocGenerics::start(gff) <- 1
  BiocGenerics::end(gff) <- BiocGenerics::width(seqs)
  BiocGenerics::strand(gff) <- "+"
  gff
}

# writing files of files -------------------------------------------------------
# create filenames
.get_file_names <- function(gfffile,
                            fafile,
                            ident){
  fileEndingFa <- ".fa"
  fileEndingGFF <- ".gff"
  fileNames <- list(in_fafile = fafile,
                    in_gfffile = gfffile)
  fileNames[["convert_gfffile"]] <- .get_file_name(gfffile, 
                                                   fileEndingGFF, 
                                                   ident)
  
  fileNames[["convert_fafile"]] <- .get_file_name(fafile, 
                                                  fileEndingFa, 
                                                  ident)
  fileNames[["mask_gfffile"]] <- .get_file_name(gfffile, 
                                                fileEndingGFF, 
                                                ident,
                                                "mask")
  fileNames[["trim_gfffile"]] <- .get_file_name(gfffile, 
                                                fileEndingGFF, 
                                                ident,
                                                "trim")
  fileNames
}
.get_file_name <- function(file, fileending, ident, modString = ""){
  if(modString != "") modString <- paste0(modString,"_")
  paste0(dirname(file),"/",
         unlist(stringr::str_split(basename(file),fileending))[1],
         "_",
         modString,
         ident,
         fileending,
         unlist(stringr::str_split(basename(file),fileending))[2])
}

# write output
.write_output <- function(gff_subset_clustered,
                          seqs_subset_clustered,
                          gff_subset,
                          fileNames){
  # save subset as single gff and fasta file
  .write_gff(gff_subset_clustered, 
             fileNames[["convert_gfffile"]],
             "Saved extracted gff3 file: ")
  .write_fa(seqs_subset_clustered, 
            fileNames[["convert_fafile"]],
            "Saved extracted fasta file: ")
  # save the masking gff annotation, which can be used for masking sequences 
  # with bedtools
  .write_gff(gff_subset, 
             fileNames[["mask_gfffile"]],
             "Saved mask gff3 file: ")
}

# write gff file
.write_gff <- function(gff_save, fileName, message){
  if(missing(message)){
    message <- ""
  }
  # remove brackets from chrom names
  chrom_names <- .fix_chrom_names(
    as.character(GenomeInfoDb::seqnames(gff_save)))
  gff_save <- GenomicRanges::GRanges(S4Vectors::Rle(chrom_names),
                          ranges = IRanges::ranges(gff_save),
                          strand = as.character(BiocGenerics::strand(gff_save)),
                          S4Vectors::mcols(gff_save))
  rtracklayer::export.gff3(gff_save,con = fileName)
  message(message,fileName)
  return(invisible(fileName))
}

# write fasta file
.write_fa <- function(seqs_save, fileName, message){
  if(missing(message)){
    message <- ""
  }
  # remove brackets from chrom names
  names(seqs_save) <- .fix_chrom_names(names(seqs_save))
  Biostrings::writeXStringSet(seqs_save,filepath = fileName)
  message(message,fileName)
  return(invisible(fileName))
}

# appending two or more fasta/gff files ----------------------------------------

#' @name appendFastaFiles
#' 
#' @title appendFastaFiles
#' 
#' @aliases appendFastaFiles appendGFF
#' 
#' @description
#' append one or more fasta or gff files
#'
#' @param ... 
#' @param ident 
#' @param msg 
#'
#' @return file name of the saved file
#' @export
#'
#' @examples
#' \donttest{
#' appendFastaFiles(fafile, 
#'                  fafile2, 
#'                  ident = "test", 
#'                  msg = "Test fasta file saved: ")
#' appendGFF(gfffile, 
#'           gfffile2, 
#'           ident = "test", 
#'           msg = "Test gff3 file saved: ")
#' }
appendFastaFiles <- function(..., ident, msg){
  files <- as.character(unlist(list(...)))
  assertive::assert_all_are_existing_files(files)
  seqs <- lapply(files, readDNAStringSet)
  seq <- do.call(c, seqs)
  return(.write_fa(seq,
                   .get_file_name(files[[1]],
                                  ".fa",
                                  ident,
                                  "append"),
                   msg))
}

#' @rdname appendFastaFiles
#' 
#' @export
#' 
#' @examples
appendGFF <- function(..., ident, msg){
  files <- as.character(unlist(list(...)))
  assertive::assert_all_are_existing_files(files)
  gffs <- lapply(files, import.gff3)
  
  # merge chrom names
  chrom_names <- lapply(gffs,function(gff){
    as.character(GenomeInfoDb::seqnames(gff))
  })
  chrom_names <- c(unlist(chrom_names))
  
  # merge ranges
  ranges <- lapply(gffs,function(gff){
    IRanges::ranges(gff)
  })
  ranges <- do.call(c, ranges)
  
  # merge strands
  strand <- lapply(gffs,function(gff){
    as.character(BiocGenerics::strand(gff))
  })
  strand <- c(unlist(strand))
  browser()
  # merge mcols
  DF <- lapply(gffs,function(gff){
    S4Vectors::mcols(gff)
  })
  # Force ID, Name, gene and Parent columns to be set even if empty
  colnames <- Reduce(intersect, lapply(DF,colnames))
  colnames <- unique(c(colnames,"ID","Name","gene","Parent"))
  DF <- lapply(DF,function(d){
    diff <- setdiff(colnames,colnames(d))
    if(length(diff) > 0){
      d_add <- S4Vectors::DataFrame(
        setNames(lapply(diff, function(names){
          rep(NA_character_,nrow(d))
        }),
        diff))
      # cannot combine character with CharacterList
      if(length(diff[diff == "Parent"]) != 0){
        d_add[,"Parent"] <- IRanges::CharacterList(rep(list(character(0)),
                                                       nrow(d)))
      }
      d <- cbind(d,d_add)
    }
    d[,colnames]
  })
  DF <- do.call(rbind, DF)
  # construct new gff file
  gff <- GenomicRanges::GRanges(S4Vectors::Rle(chrom_names),
                                ranges = ranges,
                                strand = strand,
                                DF)
  return(.write_gff(gff,
                    .get_file_name(files[[1]],
                                   ".gff",
                                   ident,
                                   "append"),
                    msg))
}