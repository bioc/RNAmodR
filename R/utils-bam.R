#' @include RNAmodR.R
NULL

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
    stats::setNames(geneNamesPerChrom[[i]],
                    paste0(names(geneNamesPerChrom[i]),
                           ":",
                           names(geneNamesPerChrom[[i]])))
  })
  names(geneNamesPerChrom) <- names(Rsamtools::bamWhich(param))
  geneNames <- stats::setNames(unlist(geneNamesPerChrom),
                               unlist(lapply(geneNamesPerChrom, names)))
  return(geneNames)
}


# return parameters to be used by scanBam
# ranges a GRanges object containing the ranges for search in the BAM file
# quality quality argument used for scanBamParam
#' @importFrom GenomeInfoDb seqnames
.assemble_scanBamParam <- function(ranges,
                                   quality,
                                   chromosomes){
  # "subset" to parent ranges aka the transcripts
  ranges <- .get_parent_annotations(ranges)
  seqnames <- unique(as.character(GenomeInfoDb::seqnames(ranges)))
  f <- seqnames %in% chromosomes
  if(length(seqnames[!f]) > 0){
    warning("Not matching chromosome identifier in gff and bam file. ",
            "Skipping data for chromosome '",
            paste(length(seqnames[!f]), collapse = "', '"),
            "'.",
            call. = FALSE)
  }
  # keep only ranges for which sequences are available
  ranges <- ranges[as.character(GenomeInfoDb::seqnames(ranges)) %in% seqnames[f],]
  rangesList <- split(ranges,
                      GenomeInfoDb::seqnames(ranges))
  # only retrieve primary alignments
  flags <- scanBamFlag(isSecondaryAlignment = FALSE)
  param <- Rsamtools::ScanBamParam(flag = flags,
                                   which = rangesList,
                                   mapqFilter = quality)
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
