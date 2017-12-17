#' @include RNAmod.R
NULL

RNAMOD_D_P_THRESHOLD <- 0.001
RNAMOD_D_SIGMA_THRESHOLD <- 5


#' @rdname parseMod
#' 
#' @description 
#' \code{mod_D}
#'
#' @return
#' @export
#'
#' @examples
setMethod(
  f = "parseMod",
  signature = signature(object = "mod_D",
                        counts = "numeric",
                        gff = "GRanges",
                        seq = "DNAString",
                        data = "DataFrame"),
  definition = function(object,
                        counts,
                        gff,
                        seq,
                        data) {
    
    
    print("D")
    
    
  }
)
