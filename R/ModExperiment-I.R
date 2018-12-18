#' @include RNAmodR.R
#' @include ModExperiment-class.R
NULL

#' @name ModInosine
#' @title ModInosine
#' @description 
#' title
NULL

setClass("ModInosine",
         contains = c("ModExperiment"),
         prototype = list(mod = "I"))

setMethod(
  f = "initialize", 
  signature = signature(.Object = "ModInosine"),
  definition = function(.Object,
                        bamfiles,
                        fasta,
                        gff) {
    .Object <- callNextMethod(.Object,
                              bamfiles,
                              fasta,
                              gff)
    return(.Object)
  }
)

.find_inosine <- function(mod,
                          data){
  browser()
  letters <- CharacterList(strsplit(as.character(sequences(data)),""))
  modifications <- mapply(
    function(m,d,l,r){
      browser()
      rownames(m) <- seq_len(width(r))
      rownames(d) <- seq_len(width(r))
      # m <- m[l == "A",]
      m <- m[l == "A" &
               m$means.G > m$means.A &
               !is.na(m$means.G) &
               !is.na(m$means.A),]
      if(nrow(m) == 0L) return(NULL)
      m
    },
    mod,
    data,
    letters,
    split(.get_parent_annotations(ranges(data)),
          seq_along(ranges(data))))
  modifications
}


#' @rdname ModInosine
#' @export
ModInosine <- function(bamfiles,
                       fasta,
                       gff){
  browser()
  ans <- new("ModInosine",
             bamfiles,
             fasta,
             gff)
  modName <- fullName(ModRNAString())[
    which(ans@mod == shortName(ModRNAString()))]
  message("Starting to search for '",
          tools::toTitleCase(modName),
          "'...")
  ans@data <- do.call("PileupPosData",
                      list(bamfiles = ans@bamfiles,
                           fasta = ans@fasta,
                           gff = ans@gff))
  ans@data@sequences <- RNAStringSet(ans@data@sequences)
  mod <- .aggregate_pile_up(ans@data)
  mod <- .find_inosine(mod,
                       ans@data)
  ans@modifications <- mod
  ans
}

