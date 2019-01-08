#' @include RNAmodR.R
#' @include Modifier-class.R
NULL

#' @name ModifierSet
#' 
#' @title ModifierSet
#' @description 
#' title
NULL

#' @rdname ModifierSet
#' @export
setClass("ModifierSet",
         contains = c("VIRTUAL",
                      "SimpleList"),
         prototype = list(elementType = "Modifier"))


setMethod("pcompareRecursively", "ModifierSet", function(x) FALSE)

setMethod("seqtype", "ModifierSet",
          function(x) modtype(new(elementType(x)))
)

setMethod(
  f = "initialize", 
  signature = signature(.Object = "ModifierSet"),
  definition = function(.Object,
                        ...) {
    callNextMethod(.Object,
                   ...)
  }
)

# not supported functions ------------------------------------------------------

setMethod(f = "relistToClass",
          signature = c(x = "ModifierSet"),
          function(x) {
            browser()
            NULL
          })

# contructor -------------------------------------------------------------------

.contains_only_Modifier <- function(x){
  classNames <- unique(vapply(x, function(z){class(z)[[1]]},character(1)))
  if(length(classNames) != 1L){
    return(FALSE)
  }
  x <- try(.norm_modifiertype(classNames), silent = TRUE)
  if (is(x, "try-error")){
    return(FALSE)
  }
  TRUE
}

.contains_only_bamfiles <- function(x){
  x <- unname(x)
  classNames <- vapply(x, function(z){class(z)[[1]]},character(1))
  if(!all(unique(classNames) %in% c("BamFileList","character","list"))){
    return(FALSE)
  }
  namedRequired <- x[classNames %in% c("character","list")]
  names <- unique(names(unlist(namedRequired)))
  if(is.null(names) || 
     !all(tolower(names) %in% c("treated","control"))){
    return(FALSE)
  }
  x <- unlist(x)
  x[classNames %in% c("BamFileList","BamFile")] <- 
    lapply(x[classNames %in% c("BamFileList","BamFile")],
           path)
  if(!all(file.exists(unlist(x)))){
    return(FALSE)
  }
  TRUE
}

.bamfiles_to_ModifierSet <- function(modifiertype,
                                     x,
                                     fasta,
                                     gff,
                                     modifications,
                                     ...){
  # check and normalize input
  fasta <- .norm_fasta(fasta, modifiertype)
  gff <- .norm_gff(gff, modifiertype)
  modifiertype <- .norm_modifiertype(modifiertype)
  if(!is.list(x)){
    x <- list(x)
  }
  names <- names(x)
  if(is.null(names)){
    names <- vector(mode = "list", length = length(x))
  }
  x <- lapply(x,
              .norm_bamfiles,
              modifiertype)
  ni <- seq_along(x)
  # do analysis by calling the Modifier classes
  x <- BiocParallel::bpmapply(
    function(i,z,n){
      if(!is.null(n)){
        message(i,". ",modifiertype," analysis '",n,"':")
      } else {
        message(i,". ",modifiertype," analysis:")
      }
      # do not parallelize further
      BiocParallel::register(SerialParam())
      #
      do.call(modifiertype,list(z,
                                fasta = fasta,
                                gff = gff,
                                modifications = modifications))
    },
    ni,
    x,
    names,
    SIMPLIFY = FALSE)
  names(x) <- names
  # pass results to ModifierSet object
  new2(.get_class_name_for_set_from_modifier_type(modifiertype),
       listData = x)
}

.modifer_to_ModifierSet <- function(x, ...){
  modifiertype <- modifiertype(x[[1]])
  new2(.get_class_name_for_set_from_modifier_type(modifiertype),
       listData = x)
}

setMethod(f = "ModifierSet",
          signature = c(x = "list"),
          function(modifiertype, x, ...) {
            if(.contains_only_Modifier(x)){
              return(.modifer_to_ModifierSet(x, ...))
            }
            if(.contains_only_bamfiles(x)){
              args <- list(...)
              return(.bamfiles_to_ModifierSet(modifiertype, 
                                              x, 
                                              args[["fasta"]], 
                                              args[["gff"]], 
                                              ...))
            }
            stop("'x' must be a list containining only elements of the same ",
                 "type\nof 'Modifer' or elements of type ('BamFileList', ",
                 "'character', 'list') which are coercible\nto a named ",
                 "BamFileList. In the latter case, elements must contain named ",
                 "vectors or lists('treated' or 'control')\nand the files ",
                 "referenced must exist. Please note, that the list a",
                 call. = FALSE)
          })

setMethod(f = "ModifierSet",
          signature = c(x = "character"),
          function(modifiertype, x, fasta, gff, ...) {
            browser()
            .bamfiles_to_ModifierSet(modifiertype, x, fasta, gff, ...)
          })

setMethod(f = "ModifierSet",
          signature = c(x = "BamFileList"),
          function(modifiertype, x, fasta, gff, ...) {
            browser()
            .bamfiles_to_ModifierSet(modifiertype, x, fasta, gff, ...)
          })

setMethod(f = "ModifierSet",
          signature = c(x = "Modifier"),
          function(modifiertype, x, ...) {
            browser()
            .modifer_to_ModifierSet(x, ...)
          })
