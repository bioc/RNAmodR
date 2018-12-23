#' @include RNAmodR.R
NULL


# BiocGeneric helper functions -------------------------------------------------

# strand related functions
.get_strand <- function(x){
  as.character(BiocGenerics::strand(x))
}
.get_unique_strand <- function(x){
  unique(.get_strand(x))
}
.is_minus_strand <- function(x) {
  all(as.logical(.get_strand(x) == "-"))
}
.is_on_minus_strand <- function(x) {
  all(.is_on_correct_strand(x,"-"))
}
# is on the minus strand?
.is_on_correct_strand <- function(x, strand) {
  as.logical(.get_strand(x) == strand) 
}
.is_on_correct_strand2 <- function(x, gr) {
  as.logical(.get_strand(x) == .get_unique_strand(gr)) 
}

# GRanges helper functions -----------------------------------------------------

.sanitize_to_seqnames <- function(x){
  gsub("\\(","-",gsub("\\)","-",x))
}

# subset to types present in MOD_TRANSCRIPT_FEATURES
.subset_mod_transcript_features <- function(gr){
  gr[S4Vectors::mcols(gr)$type %in% MOD_TRANSCRIPT_FEATURES,]
}
# subset to types present in MOD_CONTAINING_FEATURES
.subset_mod_containing_features <- function(gr){
  gr[S4Vectors::mcols(gr)$type %in% MOD_CONTAINING_FEATURES,]
}

# returns an order based on the named of the chromosome and start and end 
# coordinates
.order_GRanges <- function(gr){
  gr[order(GenomeInfoDb::seqnames(gr),
           BiocGenerics::start(gr),
           BiocGenerics::end(gr))]
}

.get_children_factor <- function(ranges){
  if(is.null(ranges$Parent)){
    ids <- ranges$ID
    return(factor(ids))
  }
  parentPos <- which(is.na(as.character(ranges$Parent)))
  parents <- ranges[parentPos]$ID
  childrenPos <- lapply(parents, function(x){.get_children(ranges,
                                                           x,
                                                           x)})
  names(parentPos) <- parents
  pos <- c(parentPos,unlist(childrenPos))
  pos <- pos[order(pos)]
  factor(names(pos))
}

.get_children <- function(ranges,
                          Parent,
                          ID){
  childPos <- which(!is.na(as.character(ranges$Parent)) &
                      as.character(ranges$Parent) == ID)
  children <- ranges[childPos]$ID
  if(length(children) == 0L) return(NULL)
  names(childPos) <- rep(Parent,length(childPos))
  childPos <- c(childPos,
                lapply(children, function(x){.get_children(ranges,
                                                           Parent,
                                                           x)}))
  childPos
}

# gets the parent of each element of the GRangesList
.get_parent_annotations <- function(grl){
  meta_parents <- metadata(grl)[["parents"]]
  if(!is.null(meta_parents)){
    return(meta_parents)
  }
  if(is(grl,"CompressedList")){
    if(!is.null(grl@unlistData$Parent)){
      ans <- grl@unlistData[is.na(as.character(grl@unlistData$Parent))]
    } else {
      ans <- grl@unlistData
    }
  } else {
    browser()
    ans <- lapply(grl,
                  function(gr){
                    if(!is.null(gr$Parent)){
                      ans <- gr[is.na(gr$Parent)]
                    } else {
                      ans <- gr
                    }
                    ans
                  })
    ans <- do.call(c,ans)
  }
  if(length(ans) != length(grl)){
    stop("Only one parent per GRanges object in GRangesList allowed.")
  }
  if(any(ans$ID != names(grl))){
    stop("Names of GRangesList and parent IDs do not match.")
  }
  return(ans)
}

# # old stuff
# 
# # get aggregated identifiers
# .get_gr_ids <- function(gr,
#                         na.rm = TRUE){
#   IDs <- unique(c(gr$ID,gr$Name,gr$gene))
#   if(na.rm){
#     IDs <- IDs[!is.na(IDs)]
#   }
#   IDs
# }
# 
# # get unique ids
# # returns all identifiers which can be used as a identifier in the Parent
# # column aggrgated by precident: ID > Name > gene
# .get_unique_identifiers <- function(gr){
#   ids <- as.character(gr$ID)
#   ids[is.na(ids)] <- as.character(gr[is.na(ids)]$Name)
#   ids[is.na(ids)] <- as.character(gr[is.na(ids)]$gene)
#   ids
# }
# 
# 
# 
# # recursive search for all childs of an ID
# .get_childs <- function(ranges,
#                         chrom,
#                         ID){
#   gr <- ranges[as.character(GenomeInfoDb::seqnames(ranges)) == chrom &
#                  !is.na(as.character(ranges$Parent)) &
#                  as.character(ranges$Parent) == ID,]
#   if(length(gr) == 0) return(GRanges())
#   # Use ID and Name for child search
#   IDs <- .get_gr_ids(gr)
#   IDs <- IDs[IDs != ID]
#   # Walk up the parent chain
#   grl <- append(list(gr),lapply(IDs, function(x){.get_childs(ranges,chrom,x)}))
#   grl <- grl[!vapply(grl,is.null,logical(1))]
#   return(unlist(GRangesList(grl)))
# }

# 
# 
# 
# # subset GRanges for highest ranking parent
# .get_parent_annotations <- function(gr,
#                                     forceSingle = FALSE,
#                                     doRecursiveSearch = FALSE,
#                                     IDs){
#   if(!doRecursiveSearch){
#     res <- gr[is.na(as.character(gr$Parent)),]
#     if(length(res) > 1 && forceSingle) 
#       stop("No distinct single parent annotation detectable.")
#     return(res)
#   }
#   if(missing(IDs)) IDs <- .get_gr_ids(gr)
#   l <- lapply(IDs, .get_parent, gr = gr)
#   l <- unique(unlist(l[!vapply(l,is.null,logical(1))]))
#   .order_GRanges(gr[(!is.na(gr$ID) & gr$ID %in% l) |
#                       (!is.na(gr$Name) & gr$Name %in% l),])
# }
# # quasi recursive search for all possible parents
# .get_parent <- function(ID,
#                         gr){
#   parent <- as.character(gr[(!is.na(gr$ID) & gr$ID == ID) |
#                               (!is.na(gr$Name) & gr$Name == ID),]$Parent)
#   parent <- parent[!is.na(parent)]
#   if(length(parent) == 0) return(ID)
#   return(unlist(lapply(parent, .get_parent, gr)))
# }
# 
# # subsets a gRanges object for entries concerning a single ID
# # - Name/ID for parent entry
# # - all childs based on Parent name equals ID
# .subset_ranges_for_unique_transcript <- function(ranges,
#                                                  ID,
#                                                  wo.childs = TRUE){
#   gr <- ranges[(is.na(S4Vectors::mcols(ranges)$ID) & 
#                   S4Vectors::mcols(ranges)$Name == ID) |
#                  (!is.na(S4Vectors::mcols(ranges)$ID) & 
#                     S4Vectors::mcols(ranges)$ID == ID),]
#   if(!wo.childs){
#     gr <- append(gr, .get_childs(ranges,
#                                  as.character(GenomeInfoDb::seqnames(gr)),
#                                  ID))
#   }
#   .order_GRanges(gr)
# }
# 
# 
# 
# # intron positions -------------------------------------------------------------
# 
# # aggregate the positions occupied by introns
# .get_intron_positions <- function(ranges, 
#                                   ID){
#   # get childs
#   childs <- .subset_ranges_for_unique_transcript(ranges,
#                                                  ID,
#                                                  wo.childs = FALSE)
#   lapply(childs[grepl("intron",childs$type),], 
#          function(intron){
#            Biostrings::start(intron):Biostrings::end(intron)
#          })
# }
# 
# # sequences --------------------------------------------------------------------
# 
# # get transcript sequence without removed
# .get_transcript_sequence <- function(ranges,ID,seq){
#   # get gr
#   gr <- .subset_ranges_for_unique_transcript(ranges, ID)
#   # get a list of introns and the position which are void
#   posToBeRemoved <- .get_intron_positions(ranges,
#                                           gr$ID)
#   # convert DNAString into character vector
#   seq <- as.character(split(seq, 1:length(seq)))
#   # set positions as names
#   if(.is_on_minus_strand(gr)){
#     names(seq) <- BiocGenerics::end(gr):BiocGenerics::start(gr)
#   } else{
#     names(seq) <- BiocGenerics::start(gr):BiocGenerics::end(gr)
#   }
#   # remove intron sequence
#   seq <- seq[!(names(seq) %in% unlist(posToBeRemoved))]
#   # set local position
#   names(seq) <- 1:length(seq)
#   return(seq)
# }
