#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

#' @name SequenceDataList
#' 
#' @title SequenceDataList
#' 
#' @description 
#' title
#' 
NULL

setClassUnion("SequenceData_OR_NULL", c("SequenceData","NULL"))

#' @rdname SequenceDataList
#' @export
setClass("SequenceDataList",
         contains = c("Vector"),
         slots = c(listData = "list",
                   elementType = "SequenceData_OR_NULL"))

# show method ------------------------------------------------------------------

setMethod("show", "SequenceDataList",
          function(object)
          {
            lo <- length(object)
            cat(classNameForDisplay(object), " of length ", lo,
                "\n", sep = "")
            if (!is.null(names(object)))
              cat(labeledLine("names", names(object)))
          })

# parallelSlotNames ------------------------------------------------------------

setMethod("parallelSlotNames", "SequenceDataList",
          function(x) c("listData", callNextMethod())
)

# accessors --------------------------------------------------------------------
setMethod("elementType", "SequenceDataList", function(x) x@elementType)
setMethod("elementNROWS", "SequenceDataList",
          function(x)
          {
            y <- as.list(x)
            if (length(y) == 0L) {
              ans <- integer(0)
              ## We must return a named integer(0) if 'x' is named
              names(ans) <- names(x)
              return(ans)
            }
            if (length(dim(y[[1L]])) < 2L)
              return(elementNROWS(y))
            return(sapply(y, NROW))
          }
)
setMethod("isEmpty", "List", function(x) all(elementNROWS(x) == 0L))
setMethod("names", "SequenceDataList", function(x) names(as.list(x)))
setReplaceMethod("names", "SequenceDataList",
                 function(x, value) {
                   names(x@listData) <- value
                   x
                 })

# constructors -----------------------------------------------------------------

# not exported. Only used internally
new_SequenceDataList_from_list <- function(Class, x, ..., mcols)
{
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
    mcols <- S4Vectors::make_zero_col_DataFrame(length(x))
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
  .check_FUN <- function(input,FUN){
    input <- lapply(input,FUN)
    first_input[[1]]
    ans <- lapply(input[seq.int(2,length(input))],
                  function(i){
                    all(first_input == i)
                  })
    all(ans)
  }
  if(.check_FUN(x,"ranges")){
    stop("Annotation data is not equal.",stop = FALSE)
  }
  if(.check_FUN(x,"sequences")){
    stop("Sequence data is not equal.",stop = FALSE)
  }
  #
  if (missing(mcols)){
    return(new2(Class, listData = x, ..., check = FALSE))
  }
  new2(Class, listData = x, ..., elementMetadata = mcols, check = FALSE)
}

#' @rdname SequenceDataList
#' @export
SequenceDataList <- function(...){
  args <- list(...)
  if (length(args) == 1L && extends(class(args[[1L]]), "SequenceData")){
    args <- args[[1L]]
  }
  new2("SequenceDataList",
       listData = args,
       check = FALSE)
}

# Validity ---------------------------------------------------------------------

.valid.SequenceDataList.listData <- function(x)
{
  elementTypeX <- elementType(x)
  if (!all(sapply(as.list(x),
                  function(xi) extends(class(xi), elementTypeX))))
    return(paste("the 'listData' slot must be a list containing",
                 elementTypeX, "objects"))
  NULL
}
.valid.SequenceDataList <- function(x)
{
  c(.valid.SequenceDataList.listData(x))
}
S4Vectors::setValidity2("SequenceDataList", .valid.SequenceDataList)


# classNameForDisplay ----------------------------------------------------------

setMethod("classNameForDisplay", "SequenceDataList",
          function(x) "SequenceDataList")


# Subsetting -------------------------------------------------------------------

# ### Assume 'x' and 'i' are parallel List objects (i.e. same length).
# ### Returns TRUE iff 'i' contains non-NA positive values that are compatible
# ### with the shape of 'x'.
# .is_valid_NL_subscript <- function(i, x)
# {
#   unlisted_i <- unlist(i, use.names=FALSE)
#   if (!is.integer(unlisted_i))
#     unlisted_i <- as.integer(unlisted_i)
#   if (anyMissingOrOutside(unlisted_i, lower=1L))
#     return(FALSE)
#   x_eltNROWS <- elementNROWS(x)
#   i_eltNROWS <- elementNROWS(i)
#   if (any(unlisted_i > rep.int(x_eltNROWS, i_eltNROWS)))
#     return(FALSE)
#   return(TRUE)
# }
# 
# ### Assume 'x' and 'i' are parallel List objects (i.e. same length).
# ### Returns the name of one of the 3 supported fast paths ("LL", "NL", "RL")
# ### or NA if no fast path can be used.
# .select_fast_path <- function(i, x)
# {
#   ## LEPType (List Element Pseudo-Type): same as "elementType" except for
#   ## RleList objects.
#   if (is(i, "RleList")) {
#     i_runvals <- runValue(i)
#     i_LEPType <- elementType(i_runvals)
#   } else {
#     i_LEPType <- elementType(i)
#   }
#   if (extends(i_LEPType, "logical")) {
#     ## 'i' is a List of logical vectors or logical-Rle objects.
#     ## We select the "LL" fast path ("Logical List").
#     return("LL")
#   }
#   if (extends(i_LEPType, "numeric")) {
#     ## 'i' is a List of numeric vectors or numeric-Rle objects.
#     if (is(i, "RleList")) {
#       i2 <- i_runvals
#     } else {
#       i2 <- i
#     }
#     if (.is_valid_NL_subscript(i2, x)) {
#       ## We select the "NL" fast path ("Number List").
#       return("NL")
#     }
#   }
#   if (extends(i_LEPType, "IntegerRanges")) {
#     ## 'i' is a List of IntegerRanges objects.
#     ## We select the "RL" fast path ("IntegerRanges List").
#     return("RL")
#   }
#   return(NA_character_)
# }
# 
# ### Assume 'x' and 'i' are parallel List objects (i.e. same length).
# ### Truncate or recycle each list element of 'i' to the length of the
# ### corresponding element in 'x'.
# .adjust_elt_lengths <- function(i, x)
# {
#   x_eltNROWS <- unname(elementNROWS(x))
#   i_eltNROWS <- unname(elementNROWS(i))
#   idx <- which(x_eltNROWS != i_eltNROWS)
#   ## FIXME: This is rough and doesn't follow exactly the truncate-or-recycle
#   ## semantic of normalizeSingleBracketSubscript() on a logical vector or
#   ## logical-Rle object.
#   for (k in idx)
#     i[[k]] <- rep(i[[k]], length.out=x_eltNROWS[k])
#   return(i)
# }
# 
# ### Assume 'x' and 'i' are parallel List objects (i.e. same length),
# ### and 'i' is a List of logical vectors or logical-Rle objects.
# .unlist_LL_subscript <- function(i, x)
# {
#   i <- .adjust_elt_lengths(i, x)
#   unlist(i, use.names=FALSE)
# }
# 
# ### Assume 'x' and 'i' are parallel List objects (i.e. same length),
# ### and 'i' is a List of numeric vectors or numeric-Rle objects.
# .unlist_NL_subscript <- function(i, x)
# {
#   offsets <- c(0L, end(IRanges::PartitioningByEnd(x))[-length(x)])
#   i <- i + offsets
#   unlist(i, use.names=FALSE)
# }
# 
# ### Assume 'x' and 'i' are parallel List objects (i.e. same length),
# ### and 'i' is a List of IntegerRanges objects.
# .unlist_RL_subscript <- function(i, x)
# {
#   unlisted_i <- unlist(i, use.names=FALSE)
#   offsets <- c(0L, end(IRanges::PartitioningByEnd(x))[-length(x)])
#   IRanges::shift(unlisted_i, shift=rep.int(offsets, elementNROWS(i)))
# }
# 
# ### Fast subset by List of logical vectors or logical-Rle objects.
# ### Assume 'x' and 'i' are parallel List objects (i.e. same length).
# ### Propagate 'names(x)' only. Caller is responsible for propagating 'mcols(x)'
# ### and 'metadata(x)'.
# .fast_subset_List_by_LL <- function(x, i)
# {
#   ## Unlist 'x' and 'i'.
#   unlisted_x <- unlist(x, use.names=FALSE)
#   unlisted_i <- .unlist_LL_subscript(i, x)
# 
#   ## Subset.
#   unlisted_ans <- extractROWS(unlisted_x, unlisted_i)
# 
#   ## Relist.
#   group <- rep.int(seq_along(x), elementNROWS(x))
#   group <- extractROWS(group, unlisted_i)
#   ans_partitioning <- IRanges::PartitioningByEnd(group, NG=length(x),
#                                                  names=names(x))
#   relist(unlisted_ans, ans_partitioning)
# }
# 
# ### Fast subset by List of numeric vectors or numeric-Rle objects.
# ### Assume 'x' and 'i' are parallel List objects (i.e. same length).
# ### Propagate 'names(x)' only. Caller is responsible for propagating 'mcols(x)'
# ### and 'metadata(x)'.
# .fast_subset_List_by_NL <- function(x, i)
# {
#   ## Unlist 'x' and 'i'.
#   unlisted_x <- unlist(x, use.names=FALSE)
#   unlisted_i <- .unlist_NL_subscript(i, x)
# 
#   ## Subset.
#   unlisted_ans <- extractROWS(unlisted_x, unlisted_i)
# 
#   ## Relist.
#   ans_breakpoints <- cumsum(unname(elementNROWS(i)))
#   ans_partitioning <- IRanges::PartitioningByEnd(ans_breakpoints,
#                                                  names=names(x))
#   relist(unlisted_ans, ans_partitioning)
# }
# 
# ### Fast subset by List of IntegerRanges objects.
# ### Assume 'x' and 'i' are parallel List objects (i.e. same length).
# ### Propagate 'names(x)' only. Caller is responsible for propagating 'mcols(x)'
# ### and 'metadata(x)'.
# .fast_subset_List_by_RL <- function(x, i)
# {
#   i_eltNROWS <- elementNROWS(i)
#   if (all(i_eltNROWS == 1L)) {
#     unlisted_i <- unlist(i, use.names=FALSE)
#     return(IRanges::windows(x, unlisted_i))
#   }
# 
#   ## Unlist 'x' and 'i'.
#   unlisted_x <- unlist(x, use.names=FALSE)
#   unlisted_i <- .unlist_RL_subscript(i, x)
# 
#   ## Subset.
#   unlisted_ans <- extractROWS(unlisted_x, unlisted_i)
# 
#   ## Relist.
#   ans_breakpoints <- cumsum(unlist(sum(width(i)), use.names=FALSE))
#   ans_partitioning <- IRanges::PartitioningByEnd(ans_breakpoints,
#                                                  names=names(x))
#   relist(unlisted_ans, ans_partitioning)
# }
# 
# ### Subset a List object by a list-like subscript.
# subset_List_by_List <- function(x, i)
# {
#   li <- length(i)
#   if (is.null(names(i))) {
#     lx <- length(x)
#     if (li > lx)
#       stop("list-like subscript is longer than ",
#            "list-like object to subset")
#     if (li < lx)
#       x <- x[seq_len(li)]
#   } else {
#     if (is.null(names(x)))
#       stop("cannot subscript an unnamed list-like object ",
#            "by a named list-like object")
#     if (!identical(names(i), names(x))) {
#       i2x <- match(names(i), names(x))
#       if (anyMissing(i2x))
#         stop("list-like subscript has names not in ",
#              "list-like object to subset")
#       x <- x[i2x]
#     }
#   }
#   ## From here, 'x' and 'i' are guaranteed to have the same length.
#   if (li == 0L)
#     return(x)
#   if (!is(x, "SimpleList")) {
#     ## We'll try to take a fast path.
#     if (is(i, "List")) {
#       fast_path <- .select_fast_path(i, x)
#     } else {
#       i2 <- as(i, "List")
#       i2_elttype <- elementType(i2)
#       if (length(i2) == li && all(sapply(i, is, i2_elttype))) {
#         fast_path <- .select_fast_path(i2, x)
#         if (!is.na(fast_path))
#           i <- i2
#       } else {
#         fast_path <- NA_character_
#       }
#     }
#     if (!is.na(fast_path)) {
#       fast_path_FUN <- switch(fast_path,
#                               LL=.fast_subset_List_by_LL,
#                               NL=.fast_subset_List_by_NL,
#                               RL=.fast_subset_List_by_RL)
#       ans <- as(fast_path_FUN(x, i), class(x))  # fast path
#       ## Propagate 'metadata(x)' and 'mcols(x)'.
#       metadata(ans) <- metadata(x)
#       mcols(ans) <- mcols(x, use.names=FALSE)
#       return(ans)
#     }
#   }
#   ## Slow path (loops over the list elements of 'x').
#   for (k in seq_len(li))
#     x[[k]] <- extractROWS(x[[k]], i[[k]])
#   return(x)
# }
# 
# .adjust_value_length <- function(value, i_len)
# {
#   value_len <- length(value)
#   if (value_len == i_len)
#     return(value)
#   if (i_len %% value_len != 0L)
#     warning("number of values supplied is not a sub-multiple ",
#             "of the number of values to be replaced")
#   rep(value, length.out=i_len)
# }
# 
# ### Assume 'x' and 'i' are parallel List objects (i.e. same length).
# .fast_lsubset_List_by_List <- function(x, i, value)
# {
#   ## Unlist 'x', 'i', and 'value'.
#   unlisted_x <- unlist(x, use.names=FALSE)
#   fast_path <- .select_fast_path(i, x)
#   unlist_subscript_FUN <- switch(fast_path,
#                                  LL=.unlist_LL_subscript,
#                                  NL=.unlist_NL_subscript,
#                                  RL=.unlist_RL_subscript)
#   unlisted_i <- unlist_subscript_FUN(i, x)
#   if (length(value) != 1L) {
#     value <- .adjust_value_length(value, length(i))
#     value <- .adjust_elt_lengths(value, i)
#   }
#   unlisted_value <- unlist(value, use.names=FALSE)
# 
#   ## Subset.
#   unlisted_ans <- replaceROWS(unlisted_x, unlisted_i, unlisted_value)
# 
#   ## Relist.
#   ans <- as(relist(unlisted_ans, x), class(x))
#   metadata(ans) <- metadata(x)
#   ans
# }
# 
# lsubset_List_by_List <- function(x, i, value)
# {
#   lx <- length(x)
#   li <- length(i)
#   if (li == 0L) {
#     ## Surprisingly, in that case, `[<-` on standard vectors does not
#     ## even look at 'value'. So neither do we...
#     return(x)
#   }
#   lv <- length(value)
#   if (lv == 0L)
#     stop("replacement has length zero")
#   value <- normalizeSingleBracketReplacementValue(value, x)
#   if (is.null(names(i))) {
#     if (li != lx)
#       stop("when list-like subscript is unnamed, it must have the ",
#            "length of list-like object to subset")
#     if (!is(x, "SimpleList")) {
#       ## We'll try to take a fast path.
#       if (is(i, "List")) {
#         fast_path <- .select_fast_path(i, x)
#       } else {
#         i2 <- as(i, "List")
#         i2_elttype <- elementType(i2)
#         if (length(i2) == li && all(sapply(i, is, i2_elttype))) {
#           fast_path <- .select_fast_path(i2, x)
#           if (!is.na(fast_path))
#             i <- i2
#         } else {
#           fast_path <- NA_character_
#         }
#       }
#       if (!is.na(fast_path))
#         return(.fast_lsubset_List_by_List(x, i, value))  # fast path
#     }
#     i2x <- seq_len(li)
#   } else {
#     if (is.null(names(x)))
#       stop("cannot subset an unnamed list-like object ",
#            "by a named list-like subscript")
#     i2x <- match(names(i), names(x))
#     if (anyMissing(i2x))
#       stop("list-like subscript has names not in ",
#            "list-like object to subset")
#     if (anyDuplicated(i2x))
#       stop("list-like subscript has duplicated names")
#   }
#   value <- .adjust_value_length(value, li)
#   ## Slow path (loops over the list elements of 'x').
#   for (k in seq_len(li))
#     x[[i2x[k]]] <- replaceROWS(x[[i2x[k]]], i[[k]], value[[k]])
#   return(x)
# }

setMethod("[", "SequenceDataList",
          function(x, i, j, ..., drop = TRUE){
            if (length(list(...)) > 0L){
              stop("invalid subsetting")
            }
            if (missing(i) || !is(i, "list_OR_List") || is(i, "IntegerRanges")) {
              ans <- S4Vectors:::subset_along_ROWS(x, i, drop = drop)
            } else {
              ans <- S4Vectors:::subset_List_by_List(x, i)
            }
            if (!missing(j)){
              mcols(ans) <- mcols(ans, use.names = FALSE)[ , j, drop = FALSE]
            }
            ans
          }
)
setReplaceMethod("[", "SequenceDataList",
                 function(x, i, j, ..., value){
                   if (!missing(j) || length(list(...)) > 0L){
                     stop("invalid subsetting")
                   }
                   if (!missing(i) && is(i, "list_OR_List") && !is(i, "IntegerRanges")){
                     return(S4Vectors:::lsubset_List_by_List(x, i, value))
                   }
                   callNextMethod(x, i, value = value)
                 }
)
setMethod("[[", "SequenceDataList",
          function(x, i, j, ...){
            dotArgs <- list(...)
            if (length(dotArgs) > 0L){
              dotArgs <- dotArgs[names(dotArgs) != "exact"]
            }
            if (!missing(j) || length(dotArgs) > 0L){
              stop("incorrect number of subscripts")
            }
            ## '...' is either empty or contains only the 'exact' arg.
            getListElement(x, i, ...)
          }
)
setMethod("$", "SequenceDataList", 
          function(x, name) x[[name, exact = FALSE]])
setReplaceMethod("[[", "SequenceDataList",
                 function(x, i, j, ..., value){
                   if (!missing(j) || length(list(...)) > 0)
                     stop("invalid replacement")
                   setListElement(x, i, value)
                 }
)
setReplaceMethod("$", "SequenceDataList",
                 function(x, name, value) {
                   x[[name]] <- value
                   x
                 })
setMethod("setListElement", "SequenceDataList",
          S4Vectors:::setListElement_default)
setMethod("getListElement", "SequenceDataList",
          function(x, i, exact = TRUE)
            getListElement(x@listData, i, exact = exact)
)

# looping ----------------------------------------------------------------------
setMethod("lapply", "SequenceDataList",
          function(X, FUN, ...) lapply(as.list(X), match.fun(FUN), ...)
)

# coercion ---------------------------------------------------------------------

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
setAs("list", "SequenceDataList",
      function(from) new_SequenceDataList_from_list("SequenceDataList", from))

# additional accessors ---------------------------------------------------------

.subaccessors <- function(x,FUN,ans_Class){
  x_not_NULL <- !vapply(x,is.null,logical(1))
  do.call(ans_Class,
          lapply(x[x_not_NULL],aggregate))
}

#' @name SequenceDataList
#' @export
setMethod(f = "gff", 
          signature = signature(x = "SequenceDataList"),
          definition = function(x){
            fasta(x[[1]])
          })

#' @name SequenceDataList
#' @export
setMethod(f = "fasta", 
          signature = signature(x = "SequenceDataList"),
          definition = function(x){
            fasta(x[[1]])
          })

#' @name SequenceDataList
#' @export
setMethod(f = "sequences", 
          signature = signature(x = "SequenceDataList"),
          definition = function(x){
            sequences(x[[1]])
          })

#' @name SequenceDataList
#' @export
setMethod(f = "ranges", 
          signature = signature(x = "SequenceDataList"),
          definition = function(x){
            ranges(x[[1]])
          })

#' @name SequenceDataList
#' @export
setMethod(f = "bamfiles", 
          signature = signature(x = "SequenceDataList"),
          definition = function(x){
            bamfiles(x[[1]])
          })

#' @name SequenceDataList
#' @export
setMethod("aggregate",
          signature = c(x = "SequenceDataList"),
          function(x){
            .subaccessors(x,"aggregate",S4Vectors::SimpleList)
          })