#' @include RNAmodR.R
NULL

# construct GRanges object from found modifications ----------------------------

#' @rdname RNAmodR-development
setMethod(f = "constructModRanges",
          signature = c(range = "GRanges", data = "DataFrame"),
          function(range, data, modType, scoreFun, source, type) {
            if(nrow(data) == 0L) {
              return(GenomicRanges::GRanges()) 
            }
            assertive::assert_is_a_non_empty_string(modType)
            assertive::assert_is_closure_function(scoreFun)
            assertive::assert_is_a_non_empty_string(source)
            assertive::assert_is_a_non_empty_string(type)
            if(!.is_valid_modType(modType)){
              stop("Modification '",modType,"' not found in the short name ",
                   "alphabet from the Modstrings package. ",
                   "'shortName(ModRNAString())'",call. = FALSE)
            }
            positions <- as.integer(rownames(data))
            scores <- do.call(scoreFun,
                              list(data))
            if(!is.list(scores)){
              stop("result of 'scoreFun' must be a list.")
            }
            data_length <- unique(vapply(scores, length, integer(1)))
            if(any(length(positions) != data_length)){
              stop("Number of positions and scores do not match.")
            }
            seqnames_unique <- unique(GenomeInfoDb::seqnames(range))
            ranges <- IRanges::IRanges(start = positions,
                                       width = 1L)
            strand <- unique(BiocGenerics::strand(range))
            mranges <- 
              GenomicRanges::GRanges(seqnames = seqnames_unique,
                                     ranges = ranges,
                                     strand = strand,
                                     seqinfo = GenomeInfoDb::seqinfo(range),
                                     mod = rep(modType, data_length),
                                     source = source,
                                     type = type,
                                     scores)
            mranges
          }
)
