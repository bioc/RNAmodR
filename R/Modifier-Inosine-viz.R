#' @include RNAmodR.R
#' @include Modifier-Inosine-class.R
NULL

RNAMODR_NUCLEOTIDE_COLOUR <- 
  c("." = biovizBase::getBioColor("RNA_BASES_N")[["N"]],
    "G" = biovizBase::getBioColor("RNA_BASES_N")[["G"]],
    "A" = biovizBase::getBioColor("RNA_BASES_N")[["A"]],
    "U" = biovizBase::getBioColor("RNA_BASES_N")[["U"]],
    "C" = biovizBase::getBioColor("RNA_BASES_N")[["C"]])

#' @rdname ModInosine
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "ModInosine",
                        coord = "GRanges"),
  definition = function(x,
                        coord,
                        type = NA,
                        window.size = 15L,
                        ...) {
    callNextMethod(x = x,
                   coord = coord,
                   type = "score",
                   window.size = window.size,
                   ...)
  }
)

setMethod(
  f = ".dataTracksByCoord",
  signature = signature(x = "ModInosine",
                        data = "GRanges"),
  definition = function(x,
                        data,
                        args) {
    requireNamespace("Gviz")
    browser()


    NULL
  }
)
