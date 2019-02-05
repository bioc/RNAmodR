#' @include RNAmodR.R
#' @include Gviz-ModifiedSequenceTrack-class.R
NULL

.get_ModRNA_bio_color <- function(){
  alphabetNames <- alphabet(ModRNAString())
  alphabet <- rep("#33FF00",length(alphabetNames))
  names(alphabet) <- alphabetNames
  rna_color <- getBioColor(type="RNA_BASES_N")
  alphabet[match(names(rna_color),names(alphabet))] <- rna_color
  alphabet
}

.get_RNA_bio_color <- function(){
  alphabet <- getBioColor(type = "RNA_ALPHABET")
  names(alphabet)[names(alphabet) == "T"] <- "U"
  base_alphabet <- getBioColor(type = "RNA_BASES_N")
  alphabet[match(names(base_alphabet),names(alphabet))] <- base_alphabet
  alphabet
}


# additional import shortcuts --------------------------------------------------

.pxResolution <- Gviz:::.pxResolution
.dpOrDefault <- Gviz:::.dpOrDefault
.fontGp <- Gviz:::.fontGp
.chrName <- Gviz:::.chrName

# ##############################################################################
# reimplementation of Gviz functions since the bloody functions contain the
# DNA alphabet hardcoded

setMethod("start", "ModifiedSequenceTrack", function(x) NULL)
setMethod("end", "ModifiedSequenceTrack", function(x) NULL)
setMethod("width", "ModifiedSequenceTrack", function(x) NULL)
setMethod("chromosome", "ModifiedSequenceTrack",
          function(GdObject) GdObject@chromosome)
setReplaceMethod("chromosome",
                 signature = signature("ModifiedSequenceTrack"),
                 function(GdObject, value){
                   GdObject@chromosome <- .chrName(value[1])
                   return(GdObject)
                 })
setMethod("genome", "ModifiedSequenceTrack", function(x) x@genome)
#' @rdname RNAmodR-internals
setMethod("length", "ModifiedSequenceTrack",
          function(x){
            if(chromosome(x) %in% seqnames(x)){
              length(x@sequence[[chromosome(x)]]) 
            } else {
              0L
            }
          })
setMethod("consolidateTrack", 
          signature(GdObject = "ModifiedSequenceTrack"),
          function(GdObject, chromosome, ...) {
            if(!is.null(chromosome))
              chromosome(GdObject) <- chromosome
            GdObject <- callNextMethod(GdObject, ...)
            return(GdObject)
          })

#' @import grid
setMethod("drawGD",
          signature = signature("ModifiedSequenceTrack"),
          definition = function(GdObject,
                                minBase,
                                maxBase,
                                prepare = FALSE,
                                ...) {
            requireNamespace("grid")
            alphabet <- displayPars(GdObject)$fontcolor
            fcol <- .dpOrDefault(GdObject,
                                 "fontcolor",
                                 alphabet)
            cex <- getPar(GdObject, "cex")
            xscale <- if(!.dpOrDefault(GdObject, "reverseStrand", FALSE)){
              c(minBase, maxBase) 
            } else {
              c(maxBase, minBase)
            }
            pushViewport(viewport(xscale = xscale,
                                  clip = TRUE,
                                  gp = .fontGp(GdObject,cex = cex)))
            if(prepare){
              pres <- .pxResolution()
              nsp <-  max(as.numeric(convertHeight(stringHeight(
                stringWidth(names(alphabet))),
                "native")))
              nsp <- nsp/pres["y"]*2
              displayPars(GdObject) <- list("neededVerticalSpace" = nsp)
              popViewport(1)
              return(invisible(GdObject))
            }
            Gviz:::imageMap(GdObject) <- NULL
            delta <- maxBase - minBase
            if(delta == 0){
              return(invisible(GdObject))
            }
            lwidth <- min(max(as.numeric(convertUnit(
              stringWidth(names(alphabet)),"inches"))),0.15)
            perLetter <- Gviz:::vpLocation()$isize["width"] / 
              (maxBase - minBase + 1)
            diff <- .pxResolution(.dpOrDefault(GdObject, "min.width", 2),
                                  coord = "x")
            if(diff > 1 || (maxBase - minBase + 1) >= 10e6){
              grid.lines(x = unit(c(minBase, maxBase), "native"),
                         y = 0.5,
                         gp = gpar(col = .dpOrDefault(GdObject,
                                                      "col",
                                                      "darkgray"),
                                   lwd = .dpOrDefault(GdObject,
                                                      "lwd",
                                                      2)))
            } else {
              sequence <- as.character(as(subseq(GdObject,
                                                 start = minBase,
                                                 end = maxBase - 1),
                                          "Rle"))
              at <- seq((minBase + 0.5), maxBase - 1 + 0.5, by = 1)
              sequence[sequence == "-"] <- ""
              if(perLetter < 0.5 && .dpOrDefault(GdObject, "add53", FALSE)){
                sequence[c(1, length(sequence))] <- ""
              }
              col <- fcol[toupper(sequence)]
              if(lwidth < perLetter && 
                 !.dpOrDefault(GdObject, "noLetters", FALSE)){
                grid.text(x = unit(at, "native"),
                          y = 0.5,
                          label = sequence,
                          rot = .dpOrDefault(GdObject, "rotation", 0),
                          gp = gpar(col = col))
              } else {
                grid.rect(x = unit(at, "native"),
                          y = 0.05,
                          width = unit(1, "native"),
                          height = 0.9,
                          gp = gpar(fill = col, col = "white"),
                          just = c(0.5, 0))
              }
            }
            ## The direction indicators
            if(.dpOrDefault(GdObject, "add53", FALSE)){
              if(.dpOrDefault(GdObject, "complement", FALSE)){
                grid.text(label = expression("3'"),
                          x = unit(minBase + 0.1, "native"),
                          just = c(0, 0.5),
                          gp = gpar(col = "#808080", cex = 0.8))
                grid.text(label = expression("5'"),
                          x = unit(maxBase - 0.1, "native"),
                          just = c(1, 0.5),
                          gp = gpar(col = "#808080", cex = 0.8))
              } else {
                grid.text(label = expression("5'"),
                          x = unit(minBase + 0.1, "native"),
                          just = c(0, 0.5),
                          gp = gpar(col = "#808080",cex = 0.8))
                grid.text(label = expression("3'"),
                          x = unit(maxBase - 0.1, "native"),
                          just = c(1, 0.5),
                          gp = gpar(col = "#808080", cex = 0.8))
              }
            }
            popViewport(1)
            return(invisible(GdObject))
          }
)

setMethod("subseq",
          signature = signature("ModifiedSequenceTrack"),
          definition = function(x,
                                start = NA,
                                end = NA,
                                width = NA) {
            padding <- "-"
            if(!is.na(start[1] + end[1] + width[1])){
              warning("All 'start', 'stop' and 'width' are provided, ignoring ",
                      "'width'")
              width <- NA
            }
            ## We want start and end to be set if width is provided
            if(!is.na(width[1])){
              if(is.na(start) && is.na(end)){
                stop("Two out of the three in 'start', 'end' and 'width' have ",
                     "to be provided",
                     call. = FALSE)
              }
              if(is.na(start)){
                start <- end - width[1] + 1
              }
              if(is.na(end)){
                end <- start + width[1] - 1
              }
            }
            w <- length(x)
            if(is.na(start)){
              start <- 1
            }
            if(w > 0){
              if(is.na(end)){
                end <- w
              }
              rstart <- max(1, start[1], na.rm=TRUE)
              rend <- max(rstart, min(end[1], w, na.rm=TRUE))
            } else {
              if(is.na(end)){
                end <- start
              }
              rend <- end
              rstart <- start
            }
            if(rend < rstart){
              stop("'end' has to be bigger than 'start'",
                   call. = FALSE)
            }
            if((rend - rstart + 1) > 10e6){
              stop("Sequence is too big! Unable to extract",
                   call. = FALSE)
            }
            finalSeq <- rep(do.call(x@seqType,list(padding)), end-start+1)
            if(chromosome(x) %in% seqnames(x) && rend > rstart){
              chrSeq <- x@sequence[[chromosome(x)]]
              seq <- subseq(chrSeq, start=rstart, end=rend)
              subseq(finalSeq,
                     ifelse(start<1, abs(start)+2, 1),
                     width=rend-rstart+1) <- seq
            }
            if(.dpOrDefault(x, "complement", FALSE)){
              finalSeq <- complement(finalSeq)
            }
            return(finalSeq)
          }
)

setMethod("seqnames", 
          signature = signature("ModifiedSequenceTrack"),
          function(x){
            as.character(names(x@sequence))
          }
)
