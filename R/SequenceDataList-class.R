#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

#' @name SequenceDataList-class
#' @aliases SequenceDataList
#' 
#' @title The SequenceDataList class
#' 
#' @description 
#' The \code{SequenceDataList} class is used to hold \code{SequenceData} objects
#' as its elements. It is derived from the 
#' \code{\link[S4Vectors:List-class]{List}}.
#' 
#' @param ... The elements to be included in the \code{SequenceDataList}.
#' 
#' @return a \code{SequenceDataList}
#' 
#' @examples
#' data(psd,package="RNAmodR")
#' data(e5sd,package="RNAmodR")
#' sdl <- SequenceDataList(psd,e5sd)
NULL

#' @rdname SequenceDataList-class
#' @export
setClass("SequenceDataList",
         contains = c("List"),
         slots = c(listData = "list"),
         prototype = list(elementType = "SequenceData"))

# show method ------------------------------------------------------------------
#' @rdname SequenceData-functions
setMethod("show", "SequenceDataList",
          function(object)
          {
            lo <- length(object)
            cat(classNameForDisplay(object), " of length ", lo,
                "\n", sep = "")
            if (!is.null(names(object)))
              cat(S4Vectors:::labeledLine("names", names(object)))
            ranges_mcols <- mcols(object@listData[[1]]@ranges@unlistData,
                                  use.names = FALSE)
            nhead <- S4Vectors::get_showHeadLines()
            ntail <- S4Vectors::get_showTailLines()
            nc <- if (is.null(ranges_mcols)) 0L else ncol(ranges_mcols)
            nr <- if (is.null(ranges_mcols)) 0L else nrow(ranges_mcols)
            nms <- rownames(ranges_mcols)
            if (nr <= (nhead + ntail + 1L)) {
              out <-
                as.matrix(format(as.data.frame(
                  lapply(ranges_mcols, showAsCell),
                  optional = TRUE)))
            } else {
              out <-
                rbind(as.matrix(format(as.data.frame(
                  lapply(ranges_mcols, function(x)
                    showAsCell(head(x, nhead))),
                  optional = TRUE))),
                  rbind(rep.int("...", nc)),
                  as.matrix(format(as.data.frame(
                    lapply(ranges_mcols, function(x) 
                      showAsCell(tail(x, ntail))),
                    optional = TRUE))))
              rownames(out) <- S4Vectors:::.rownames(nms, nr, nhead, ntail) 
            }
            classinfo <-
              matrix(unlist(lapply(ranges_mcols, function(x) {
                paste0("<", classNameForDisplay(x)[1],
                       ">")
              }), use.names = FALSE), nrow = 1,
              dimnames = list("", colnames(out)))
            out <- rbind(classinfo, out)
            cat("- Ranges metadata columns:\n")
            print(out, quote = FALSE, right = TRUE)
          })

# parallelSlotNames ------------------------------------------------------------
#' @rdname RNAmodR-internals
setMethod("parallelSlotNames", "SequenceDataList",
          function(x) c("listData", callNextMethod())
)

# accessors --------------------------------------------------------------------
#' @rdname SequenceData-functions
setMethod("names", "SequenceDataList", function(x) names(as.list(x)))
#' @rdname SequenceData-functions
setReplaceMethod("names", "SequenceDataList",
                 function(x, value) {
                   names(x@listData) <- value
                   x
                 })

# constructors -----------------------------------------------------------------
.compare_element_metadata <- function(input,FUN){
  input <- lapply(input,FUN)
  first_input <- input[[1]]
  ans <- vapply(input[seq.int(2,length(input))],
                function(i){
                  all(all(first_input == i))
                },
                logical(1))
  all(ans)
}

# not exported. Only used internally
new_SequenceDataList_from_list <- function(Class, x, ..., mcols){
  if (!extends(Class, "SequenceDataList")){
    stop("class ", Class, " must extend SequenceDataList")
  }
  if (!is.list(x)){
    stop("'x' must be a list")
  }
  if (is.array(x)) { # drop any unwanted dimensions
    tmp_names <- names(x)
    dim(x) <- NULL # clears the names
    names(x) <- tmp_names
  }
  class(x) <- "list"
  proto <- new(Class)
  ans_elementType <- elementType(proto)
  if (is(S4Vectors::mcols(proto, use.names = FALSE), "DataFrame")){
    mcols <- S4Vectors:::make_zero_col_DataFrame(length(x))
  }
  extends_elementType <- vapply(x,
                                function(xi){
                                  extends(class(xi), ans_elementType)
                                },
                                logical(1))
  if (!all(extends_elementType)){
    stop("all elements in 'x' must be ", ans_elementType, " objects")
  }
  # check that all bamfiles, sequences and annotation information are the same
  if(!.compare_element_metadata(x,"ranges")){
    stop("Annotation data is not equal.",stop = FALSE)
  }
  if(!.compare_element_metadata(x,"sequences")){
    stop("Sequence data is not equal.",stop = FALSE)
  }
  # class name as default names
  if(is.null(names(x))){
    names(x) <- vapply(x,class,character(1))
  }
  #
  if (missing(mcols)){
    return(new2(Class, listData = x, ..., check = FALSE))
  }
  new2(Class, listData = x, ..., elementMetadata = mcols, check = FALSE)
}

#' @rdname SequenceDataList-class
#' @export
SequenceDataList <- function(...){
  args <- list(...)
  if (length(args) == 1L && extends(class(args[[1L]]), "SequenceData")){
    args <- args[[1L]]
  }
  new_SequenceDataList_from_list("SequenceDataList", args)
}

# Validity ---------------------------------------------------------------------
.valid.SequenceDataList.listData <- function(x){
  elementTypeX <- elementType(x)
  if (!all(vapply(as.list(x),
                  function(xi) extends(class(xi), elementTypeX),
                  logical(1)))){
    return(paste("the 'listData' slot must be a list containing",
                 elementTypeX, "objects"))
  }
  if(!.compare_element_metadata(x,"ranges")){
    return("Annotation data is not equal.")
  }
  if(!.compare_element_metadata(x,"sequences")){
    return("Sequence data is not equal.")
  }
  NULL
}
.valid.SequenceDataList <- function(x){
  c(.valid.SequenceDataList.listData(x),
    unlist(lapply(x,validObject)))
}
S4Vectors::setValidity2("SequenceDataList", .valid.SequenceDataList)


# classNameForDisplay ----------------------------------------------------------
setMethod("classNameForDisplay", "SequenceDataList",
          function(x) "SequenceDataList")


# Subsetting -------------------------------------------------------------------
#' @rdname RNAmodR-internals
setMethod("getListElement", "SequenceDataList",
          function(x, i, exact = TRUE)
            getListElement(x@listData, i, exact = exact)
)

# looping ----------------------------------------------------------------------
setMethod("lapply", "SequenceDataList",
          function(X, FUN, ...) lapply(as.list(X), match.fun(FUN), ...)
)

# coercion ---------------------------------------------------------------------
# forth and ...
setAs("SequenceDataList", "list", function(from) as.list(from))
.as.list.SequenceDataList <- function(x, use.names = TRUE)
{
  if (!isTRUEorFALSE(use.names)){
    stop("'use.names' must be TRUE or FALSE")
  }
  ans <- x@listData
  if (!use.names){
    names(ans) <- NULL
  }
  ans
}
setMethod("as.list", "SequenceDataList", .as.list.SequenceDataList)

# ... back
setAs("list", "SequenceDataList",
      function(from) new_SequenceDataList_from_list("SequenceDataList", from))
setAs("ANY", "SequenceDataList", function(from) {
  new_SequenceDataList_from_list("SequenceDataList", as.list(from))
})

# additional accessors ---------------------------------------------------------
.subaccessors <- function(x,FUN,ans_Class){
  x_not_NULL <- !vapply(x,is.null,logical(1))
  do.call(ans_Class,
          lapply(x[x_not_NULL],aggregate))
}

#' @rdname SequenceData-functions
#' @export
setMethod(f = "seqinfo", 
          signature = signature(x = "SequenceDataList"),
          definition = function(x){
            seqinfo(x[[1]])
          })

#' @rdname SequenceData-functions
#' @export
setMethod(f = "names", 
          signature = signature(x = "SequenceDataList"),
          definition = function(x){
            names(x[[1]])
          })

#' @rdname SequenceData-functions
#' @export
setMethod(f = "sequences", 
          signature = signature(x = "SequenceDataList"),
          definition = function(x){
            sequences(x[[1]])
          })

#' @rdname SequenceData-functions
#' @export
setMethod(f = "ranges", 
          signature = signature(x = "SequenceDataList"),
          definition = function(x){
            ranges(x[[1]])
          })

#' @rdname SequenceData-functions
#' @export
setMethod(f = "bamfiles", 
          signature = signature(x = "SequenceDataList"),
          definition = function(x){
            bamfiles(x[[1]])
          })

#' @rdname aggregate
#' @export
setMethod("aggregate",
          signature = c(x = "SequenceDataList"),
          function(x, condition = "Treated"){
            ans <- do.call(S4Vectors::SimpleList,
                           lapply(x,aggregate, condition = condition))
            names(ans) <- names(x@listData)
            ans
          })
