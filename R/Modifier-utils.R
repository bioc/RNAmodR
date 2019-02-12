#' @include RNAmodR.R
NULL

# error propagation ------------------------------------------------------------

################################################################################
# modified from Lee Pang, 2015,
# (https://oddhypothesis.blogspot.com/2015/01/easy-error-propagation-in-r.html)
# (https://www.r-bloggers.com/easy-error-propagation-in-r/)
################################################################################

#' @importFrom dplyr mutate_
#' @importFrom stats D
.mutate_with_error <- function(.data, f){
  exprs = list(
    # expression to compute new variable values
    deparse(f[[3]]),
    # expression to compute new variable errors
    sprintf('sqrt(%s)',
            paste(vapply(all.vars(f[[3]]),
                         function(v) {
                           dfdp = deparse(stats::D(f[[3]], 
                                                   v))
                           sprintf('(d%s*(%s))^2', 
                                   v,
                                   dfdp)
                         },
                         character(1)),
                  collapse = '+'))
  )
  names(exprs) = c(
    deparse(f[[2]]),
    sprintf('d%s', deparse(f[[2]]))
  )
  if( nrow(.data) > 0 ){
    # the standard evaluation alternative of mutate()
    .data <- dplyr::mutate_(.data, .dots = exprs)
  }
  return(.data)
}

# construct GRanges object from found modifications ----------------------------

#' @rdname RNAmodR-internals
setMethod(f = ".constructModRanges",
          signature = c(range = "GRanges", data = "DataFrame"),
          function(range, data, modType, scoreFun, source, type) {
            if(nrow(data) == 0L) {
              return(GenomicRanges::GRanges()) 
            }
            positions <- as.integer(rownames(data))
            if(as.character(GenomicRanges::strand(range)) == "-"){
              positions <- GenomicRanges::end(range) - positions + 1L
            } else {
              positions <- GenomicRanges::start(range) + positions - 1L
            }
            mranges <- do.call(
              GenomicRanges::GRanges,
              c(list(seqnames = rep(as.character(GenomicRanges::seqnames(range)),
                                    nrow(data)),
                     ranges = IRanges::IRanges(start = positions,
                                               width = 1L),
                     strand = GenomicRanges::strand(range),
                     seqinfo = GenomeInfoDb::seqinfo(range),
                     mod = rep(modType,nrow(data))),
                source = source,
                type = type,
                do.call(scoreFun,
                        list(data))))
            mranges
          }
)
