#' @include RNAmodR.R
#' @include Modifier-class.R
NULL

#' @name ModInosine
#' @aliases Inosine ModifierInosine
#' 
#' @title ModInosine
#' @description 
#' title
NULL

#' @rdname ModInosine
#' @export
setClass("ModInosine",
         contains = c("Modifier"),
         prototype = list(mod = "I",
                          dataType = "PileupPosData"))
POS_DATA_TYPE_INOSINE <- "PileupPosData"


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

.ModInosineFromCharacter <- function(bamfiles,
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
  ans@data <- do.call(ans@dataType,
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

.ModInosineFromPosData <- function(x){
  ans <- new("ModInosine",
             x@bamfiles,
             x@fasta,
             x@gff)
  # check data type
  ans@data <- .norm_data_type(ans,x)
  ans
}

setGeneric( 
  name = "ModInosine",
  def = function(x,
                 ...) standardGeneric("ModInosine")
)
# Create Modifier class from file character, fasta and gff file
#' @rdname ModInosine
#' @export
setMethod("ModInosine",
          signature = c(x = "character"),
          function(x,
                   fasta,
                   gff,
                   modifications = NULL){
            ans <- .ModInosineFromCharacter(x,
                                            fasta,
                                            gff)
            ans@modifications <- .norm_modifications(ans,modifications)
            ans
          }
)

# Create Modifier class from bamfiles, fasta and gff file
#' @rdname ModInosine
#' @export
setMethod("ModInosine",
          signature = c(x = "BamFileList"),
          function(x,
                   fasta,
                   gff,
                   modifications = NULL){
            ans <- .ModInosineFromCharacter(x,
                                            fasta,
                                            gff)
            ans@modifications <- .norm_modifications(ans,modifications)
            ans
          }
)

# Create Modifier class from existing PosData
#' @rdname ModInosine
#' @export
setMethod("ModInosine",
          signature = c(x = "PosData"),
          function(x,
                   modifications = NULL){
            ans <- .ModInosineFromPosData(x)
            ans@modifications <- .norm_modifications(ans,modifications)
            ans
          }
)