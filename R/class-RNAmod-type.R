#' @include RNAmod.R
NULL



#' @name mod
#' 
#' @title mod
#'
#' @return
#' @export
#'
#' @examples
setClass("mod",
         contains = "VIRTUAL",
         slots = c(plotType = "character",
                   modType = "character"),
         prototype = list(plotType = "default")
)

setMethod(
  f = "show", 
  signature = signature(object = "mod"),
  definition = function(object) {
    
  }
)


#' @rdname mod-accessors
#'
#' @return character defining the plot type for this modification class
#' @export
#'
#' @examples
#' \donttest{
#' getPlotType(mod)
#' getModType(mod)
#' }
setMethod(
  f = "getPlotType", 
  signature = signature(object = "mod"),
  definition = function(object) {
    return(object@plotType)
  }
)

#' @rdname mod-accessors
#'
#' @return character defining the modification type for this modification
#' class
#' @export
setMethod(
  f = "getModType", 
  signature = signature(object = "mod"),
  definition = function(object) {
    return(object@modType)
  }
)




.plot_sample_data <- function(baseData,testData,name){
  requireNamespace("ggplot2", quietly = TRUE)
  
  # merge data
  df <- data.frame(x = c(rep(name,length(testData) + length(baseData))),
                   y = c(testData, baseData),
                   group = c(rep("pos",length(testData)), rep("base", length(baseData))))
  
  # plot data and create grid
  plot <- ggplot(df, aes_(x = ~x, y = ~y, colour = ~group)) +
    geom_violin(trim = FALSE) + 
    geom_jitter(width = 0.1) +
    scale_x_discrete(name = "position") + 
    scale_colour_brewer(palette = "Set1",
                        name = "position\ndata",
                        label = c("base" = "Base data",
                                  "pos" = "Position data"))
  plot2 <- plot + scale_y_log10(name = "counts per position", labels = scales::scientific)
  plot <- plot + scale_y_continuous(name = "counts per position", labels = scales::scientific)
  grid <- gridExtra::grid.arrange(plot, plot2, nrow = 1, ncol = 2)
  
  # Resurrect calling RNAmod object to get the output folder
  object <- get(".Object", envir = .where(".Object"))
  
  # create file path 
  folder <- paste0(getOutputFolder(object),
                   "sample/")
  if(!assertive::is_dir(folder)){
    dir.create(folder, recursive = TRUE)
  }
  fileName <- paste0(folder,
                     name,
                     ".pdf")
  
  # save plots as pdf
  ggsave(plot = grid,
         filename = fileName,
         width = 6,
         height = 4)
}

# converts the position on a transcript from start to end, to the correct
# genomic position
.convert_local_to_global_locations <- function(gff, loc){
  browser()
  strand <- unique(as.character(strand(gff)))
  if(strand == "-"){
    locations <- (end(gff) - loc) + 1
  } else {
    locations <- (start(gff) + loc) - 1
  }
  return(locations)
}