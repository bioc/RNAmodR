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
#' @param appendToOriginal 
#'
#' @return
#' @export
#' 
#' @importFrom S4Vectors mcols Rle
#' @importFrom rtracklayer import.gff3 export.gff3
#' @importFrom Rsamtools FaFile indexFa getSeq
#' @importFrom GenomeInfoDb seqnames
#' @importFrom stringr str_split
#' @importFrom BiocGenerics which unlist width start end strand
#' @import Biostrings
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
                                   appendToOriginal = FALSE){
  # Check inputs
  assertive::assert_all_are_existing_files(gfffile)
  assertive::assert_all_are_existing_files(fafile)
  assertive::assert_is_a_non_empty_string(ident)
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
  convertGeneToChrom(gfffile,fafile,IDs,ident,appendToOriginal)
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
                               appendToOriginal = FALSE){
  # Check inputs
  assertive::assert_all_are_existing_files(gfffile)
  assertive::assert_all_are_existing_files(fafile)
  assertive::assert_is_a_non_empty_string(ident)
  assertive::assert_is_a_bool(appendToOriginal)
  
  # Load files
  gff <- rtracklayer::import.gff3(gfffile)
  fa <- Rsamtools::FaFile(fafile)
  
  # check that genes are in ID or Name column of gff file
  checkGenes <- c(unique(S4Vectors::mcols(gff)$ID),
                  unique(S4Vectors::mcols(gff)$Name))
  checkGenes <- unique(checkGenes)
  assertive::assert_is_subset(genes, checkGenes)
  # browser()
  # create gff subset
  gff_sub <- gff[(is.na(S4Vectors::mcols(gff)$ID) & 
                    S4Vectors::mcols(gff)$Name %in% genes) |
                   (!is.na(S4Vectors::mcols(gff)$ID) & 
                      S4Vectors::mcols(gff)$ID %in% genes) |
                   (!is.na(as.character(S4Vectors::mcols(gff)$Parent)) & 
                      as.character(S4Vectors::mcols(gff)$Parent) %in% genes) |
                   (is.na(S4Vectors::mcols(gff)$ID) & 
                      S4Vectors::mcols(gff)$Name %in% paste0(genes,"_mRNA")) |
                   (!is.na(S4Vectors::mcols(gff)$ID) & 
                      S4Vectors::mcols(gff)$ID %in% paste0(genes,"_mRNA")) |
                   (!is.na(as.character(S4Vectors::mcols(gff)$Parent)) & 
                      as.character(S4Vectors::mcols(gff)$Parent) %in% paste0(genes,"_mRNA")),]
  # Remove all entries part of the subset
  gff_updated <- gff[(!is.na(S4Vectors::mcols(gff)$Name) &
                        !(S4Vectors::mcols(gff)$Name %in% S4Vectors::mcols(gff_sub)$Name)) |
                       (!is.na(S4Vectors::mcols(gff)$ID) &
                          !(S4Vectors::mcols(gff)$ID %in% S4Vectors::mcols(gff_sub)$ID)),]
  
  # extract the coding sequences of all features
  seqs <- .extract_exon_sequences(gff_sub,fa)
  
  # remove indistinguishable sequences (applies usually to tRNA, rRNA and such)
  gff_exon <- .cluster_results_gff(gff_sub,seqs)
  fa_seq_exon <- .cluster_results_seqs(seqs)
  
  # create new gff annotation and fasta sequences
  gff_exon <- .update_gff_subset(gff_exon,fa_seq_exon)
  fa_seq_exon <- .update_seqs(fa_seq_exon)
  
  # save input file names
  resultFiles <- c()
  resultFiles[["in_gfffile"]] <- gfffile
  resultFiles[["in_fafile"]] <- fafile
  
  # save subset as single gff and fasta file
  resultFiles[["convert_gfffile"]] <- .write_gff(gff_exon, 
                                                 ident, 
                                                 gfffile,
                                                 "Saved extracted gff3 file: ")
  resultFiles[["mask_gfffile"]] <- .write_gff(gff_sub, 
                                              paste0("mask_",
                                                     ident), 
                                              gfffile,
                                              "Saved mask gff3 file: ")
  resultFiles[["convert_fafile"]] <- .write_fa(fa_seq_exon, 
                                               ident, 
                                               fafile,
                                               "Saved extracted fasta file: ")
  
  if(appendToOriginal){
    resultFiles[["out_fafile"]] <- appendFastaFiles(resultFiles[["in_fafile"]],
                                                    resultFiles[["convert_fafile"]],
                                                    ident = ident,
                                                    msg = "Saved appended fasta file: ")
    resultFiles[["out_gfffile"]] <- appendGFF(resultFiles[["in_gfffile"]],
                                              resultFiles[["convert_gfffile"]],
                                              ident = ident,
                                              msg = "Saved appended gff3 file: ")
  }
  return(invisible(resultFiles))
}

# merge seq_identifier
.get_unique_seqnames <- function(gff){
  res <- lapply(gff, function(g){
    as.character(paste0(S4Vectors::mcols(g)$ID,
                        RNAMOD_SEP_SEQNAMES,
                        S4Vectors::mcols(g)$Name,
                        RNAMOD_SEP_SEQNAMES,
                        S4Vectors::mcols(g)$gene))
  })
  return(unlist(res))
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

# extract the exons sequences and merges them into one per gene
.extract_exon_sequences <- function(gff_sub, fa){
  # subset to exons only
  exons <- gff_sub[S4Vectors::mcols(gff_sub)$type %in% c("CDS",
                                                         "noncoding_exon",
                                                         "exon")]
  # get sequences for exons
  exons_seq <- Biostrings::getSeq(fa,exons)
  names(exons_seq) <- .get_unique_seqnames(exons)
  
  # get unqiue names and iterate and merge split sequences
  uniq_gene_names <- unique(names(exons_seq))
  seqs <- Biostrings::DNAStringSet(lapply(uniq_gene_names, function(name){
    seq <- exons_seq[BiocGenerics::which(names(exons_seq) == name)]
    return(BiocGenerics::unlist(seq))
  }))
  # rename sequences to match parent sequence
  names(seqs) <- uniq_gene_names
  # names(seqs) <- gsub("_CDS","",names(seqs))
  # names(seqs) <- gsub("_noncoding_exon","",names(seqs))
  return(seqs)
}

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



# create new gff annotations
.update_gff_subset <- function(gff_sub, seqs){
  
  
  # subset gff to parents and childs
  gff_children <- gff_sub[.get_unique_seqnames(gff_sub) %in% names(seqs)]
  gff_children <- gff_children[!duplicated(as.character(S4Vectors::mcols(gff_children)$Name))]
  gff_parents <- gff_sub[(!is.na(as.character(S4Vectors::mcols(gff_sub)$ID)) &
                           as.character(S4Vectors::mcols(gff_sub)$ID) %in% as.character(S4Vectors::mcols(gff_children)$Parent)) |
                           (!is.na(as.character(S4Vectors::mcols(gff_sub)$Name)) &
                           as.character(S4Vectors::mcols(gff_sub)$Name) %in% as.character(S4Vectors::mcols(gff_children)$Parent))]
  
  if( length(gff_children) != length(gff_parents) |
      length(seqs) != length(gff_children)){
    warnings("Mismatching number of annotations and sequences",
             .call = FALSe)
  }
  
  # setup coordinates for genes based on seqs as chromosomes
  BiocGenerics::start(gff_parents) <- 1
  BiocGenerics::end(gff_parents) <- BiocGenerics::width(seqs)
  BiocGenerics::strand(gff_parents) <- "+"
  
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
  res <- gsub("(?![a-zA-Z0-9.:^*$@!+_?-|]).","-",res, perl = TRUE)
  return(res)
}

# write gff file
.write_gff <- function(gff, ident, gfffile, message){
  fileEnding <- ".gff3"
  fileName <- paste0(dirname(gfffile),"/",
                     unlist(stringr::str_split(basename(gfffile),fileEnding))[1],
                     "_",
                     ident,
                     fileEnding)
  
  # remove brackets from chrom names
  chrom_names <- as.character(GenomeInfoDb::seqnames(gff))
  gff <- GenomicRanges::GRanges(S4Vectors::Rle(chrom_names),
                                ranges = IRanges::ranges(gff),
                                strand = as.character(BiocGenerics::strand(gff)),
                                S4Vectors::mcols(gff))
  
  rtracklayer::export.gff3(gff,con = fileName)
  message(message,fileName)
  return(fileName)
}

# write fasta file
.write_fa <- function(seqs, ident, fafile, message){
  fileEnding <- ".fa"
  fileName <- paste0(dirname(fafile),"/",
                     unlist(stringr::str_split(basename(fafile),fileEnding))[1],
                     "_",
                     ident,
                     fileEnding,
                     unlist(stringr::str_split(basename(fafile),fileEnding))[2])
  Biostrings::writeXStringSet(seqs,filepath = fileName)
  message(message,fileName)
  return(fileName)
}


#' @name appendFastaFiles
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
  fas <- lapply(files, readDNAStringSet)
  fa <- do.call(append, fas)
  
  ident <- paste0("merge_",ident)
  
  return(.write_fa(fa,ident,files[[1]],msg))
}

#' @rdname appendFastaFiles
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
  
  # merge mcols
  DF <- lapply(gffs,function(gff){
    S4Vectors::mcols(gff)
  })
  colnames <- do.call(intersect, lapply(DF,colnames))
  DF <- lapply(DF,function(d){
    d[,colnames]
  })
  DF <- do.call(rbind, DF)
  
  # construct new gff file
  gff <- GenomicRanges::GRanges(S4Vectors::Rle(chrom_names),
                                ranges = ranges,
                                strand = strand,
                                DF)
  
  ident <- paste0("merge_",ident)
  
  return(.write_gff(gff,ident,files[[1]],msg))
}

maskGene <- function(gfffile,
                     fafile,
                     gff_updated,
                     gff_removed,
                     ident){
  # Check inputs
  assertive::assert_all_are_existing_files(gfffile)
  assertive::assert_all_are_existing_files(fafile)
  
  
  out_fafile
  
  res <- list(out_gfffile = .write_gff(gff_updated, 
                                       paste0("masked_",ident), 
                                       gfffile,
                                       "Saved masked gff3 file: "),
              out_fafile = out_fafile)
}