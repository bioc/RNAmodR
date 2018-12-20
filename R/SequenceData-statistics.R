#' @include RNAmodR.R
NULL


# classes for representing statistical information during search for modifications

# - derive from GPos object since it will be tied to single sequence positions
# - different types of class for t-test, MWU and such, depends on what is needed
# - derive class
# - add stat function as slot
# - add stat function depending on type of class
# - keep only data and calculate statistic on demand?

# classes for storing stational data in final result

# - same classes probably, since they do no differ


# setMethod("ttest",
#           signature = "ModData",
#           function(x){
#             
#           })
