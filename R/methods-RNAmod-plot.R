#' @include RNAmod.R
NULL

RNAMOD_PLOT_POS_WIDTH <- 1
RNAMOD_PLOT_LAYER_HEIGHT <- 50
RNAMOD_PLOT_DATA_HEIGHT <- 100
RNAMOD_PLOT_FOCUS_WINDOW <- 50
RNAMOD_PLOT_MM_TO_INCH_F <- 0.03937008

#' @name getModPlot
#' 
#' @aliases getModPlot saveModPlots
#' 
#' @title visualization of modifications
#' 
#' @param .Object a RNAmod object. 
#' @param number a number defining an experiment. 
#' @param modifications name of modification to be used for analysis. 
#' @param gene a single gene name
#'
#' @return a grid of layer plot and any mod plot
#' @export
#' 
#' @import ggplot2
#'
#' @examples
#' \donttest{
#' getModPlot(mod,1,"m7G","RDN18-1")
#' saveModPlot(mod,1,"m7G",c("RDN18-1","tS\(CGA\)P"))
#' }
setMethod(
  f = "getModPlot", 
  signature = signature(.Object = "RNAmod",
                        se = "SummarizedExperiment",
                        modifications = "character",
                        gene = "character"),
  definition = function(.Object,
                        se,
                        modifications,
                        gene,
                        focus){
    # Check input
    .check_plot_input(se,
                      modifications,
                      gene,
                      focus,
                      "pdf")
    
    # get intersection of requested and available genes
    genesAvail <- intersect(rownames(se), gene)
    if(length(genesAvail) != length(gene)){
      if(length(genesAvail) > 0){
        warning("No information available for the following genes: '",
                paste0(setdiff(genes,genesAvail), collapse = "', '"),"'",
                call. = FALSE)
      } else {
        stop("No information available for any genes given: '",
             paste0(genes, collapse = "', '"),"'",
             call. = FALSE)
      }
    }
    genesAvail <- unlist(genesAvail)
    assertive::assert_is_a_non_empty_string(genesAvail)
    
    positions <- SummarizedExperiment::rowData(se[rownames(se) %in% genesAvail,])$positions[[1]]
    mods <- SummarizedExperiment::rowData(se[rownames(se) %in% genesAvail,])$mods[[1]]
    gff <- .Object@.dataGFF
    fasta <- .Object@.dataFasta
    
    modClasses <- .load_mod_classes(modifications)
    if(!focus){
      plot <- .plot_gene_with_modifications(genesAvail,
                                            positions,
                                            mods,
                                            gff,
                                            fasta,
                                            modClasses)
      plot <- plot$plot
    } else {
      plot <- .plot_per_modifications(genesAvail,
                                      positions,
                                      mods,
                                      gff,
                                      fasta,
                                      modClasses)
      plot <- plot$plot[[1]]
    }
    return(plot)
  }
)

.check_plot_input <- function(se,
                              modifications,
                              genes,
                              focus,
                              filetype){
  assertive::assert_all_are_non_empty_character(modifications)
  assertive::assert_all_are_non_empty_character(genes)
  assertive::assert_is_a_bool(focus)
  checkFileTypes <- c("pdf","png")
  if(!assertive::is_subset(filetype, checkFileTypes))
    stop("Unsupported file type '",filetype,"'",
         .call = FALSE)
}


#' @rdname getModPlot
#'
#' @return TRUE if successful
#' @export
#' 
#' @import ggplot2
setMethod(
  f = "saveModPlot", 
  signature = signature(.Object = "RNAmod",
                        se = "SummarizedExperiment",
                        modifications = "character",
                        genes = "character"),
  definition = function(.Object,
                        se,
                        modifications,
                        genes,
                        focus,
                        filetype){
    # Check input
    .check_plot_input(se,
                      modifications,
                      genes,
                      focus,
                      filetype)
    
    # create folder
    folder <- paste0(getOutputFolder(.Object),
                     "ModPlots/",
                     unique(SummarizedExperiment::colData(se)$SampleName),
                     "/")
    if(!assertive::is_dir(folder)){
      dir.create(folder, recursive = TRUE)
    }
    
    # get intersection of requested and available genes
    genesAvail <- intersect(rownames(se), genes)
    if(length(genesAvail) != length(genes)){
      if(length(genesAvail) > 0){
        warning("No information available for the following genes: '",
                paste0(setdiff(genes,genesAvail), collapse = "', '"),"'",
                call. = FALSE)
      } else {
        stop("No information available for any genes given: '",
             paste0(genes, collapse = "', '"),"'",
             call. = FALSE)
      }
    }
    
    positions <- SummarizedExperiment::rowData(se[rownames(se) %in% genesAvail,])$positions[[1]]
    mods <- SummarizedExperiment::rowData(se[rownames(se) %in% genesAvail,])$mods[[1]]
    gff <- .Object@.dataGFF
    fasta <- .Object@.dataFasta
    
    modClasses <- .load_mod_classes(modifications)
    
    if(!focus){
      for(i in seq_along(genesAvail)){
        plot <- .plot_gene_with_modifications(genesAvail[[i]],
                                              positions,
                                              mods,
                                              gff,
                                              fasta,
                                              modClasses)
        names <- paste0(paste(modifications, collapse = "-"),
                        "_",
                        genesAvail[[i]])
        .save_gene_mod_plot(filetype,
                            folder,
                            names,
                            list(plot$plot),
                            list(plot$width),
                            list(plot$height))
      }
    } else {
      for(i in seq_along(genesAvail)){
        plot <- .plot_per_modifications(genesAvail[[i]],
                                        positions,
                                        mods,
                                        gff,
                                        fasta,
                                        modClasses)
        names <- paste0(genesAvail[[i]],
                        "_",
                        plot$names)
        .save_gene_mod_plot(filetype,
                            folder,
                            names,
                            plot$plot,
                            plot$width,
                            plot$height)
      }
    }
    return(invisible(TRUE))
  }
)

# plotting a complete gene -----------------------------------------------------

.plot_gene_with_modifications <- function(geneName,
                                          positions,
                                          mods,
                                          gff,
                                          fasta,
                                          modClasses){
  gff_sub <- gff[(is.na(S4Vectors::mcols(gff)$ID) & 
                    S4Vectors::mcols(gff)$Name == geneName) |
                   (!is.na(S4Vectors::mcols(gff)$ID) & 
                      S4Vectors::mcols(gff)$ID == geneName),]
  seq <- Rsamtools::getSeq(fasta,gff_sub)[[1]]
  
  # create description layer plot
  layer <- .create_layer_data(gff,gff_sub)
  
  # get data for different plot types
  posData <- .aggregate_pos_data(gff_sub,positions,modClasses)
  modData <- .aggregate_mod_data(mods,modClasses)
  
  # to inhibit plotting by grid.arrange
  grDevices::pdf(file = NULL)
  
  dataPlots <- lapply(names(posData), function(type){
    pos <- posData[[type]]
    mods <- modData[[type]]
    
    plot <- .get_mod_plot(pos,mods)
    return(plot)
  })
  layerPlot <- .get_gene_plot(geneName, layer, posData[[1]], seq)
  #layerPlot <- list()
  
  nrow <- length(layerPlot) + length(dataPlots)
  width <- 40 + nrow(posData[[1]]) * RNAMOD_PLOT_POS_WIDTH
  height <- length(layerPlot) * RNAMOD_PLOT_LAYER_HEIGHT + 
    length(dataPlots) * RNAMOD_PLOT_DATA_HEIGHT
  
  # grid <- do.call(gridExtra::grid.arrange, append(dataPlots,
  #                                                 layerPlot,
  #                                                 list(nrow = nrow)))
  grid <- do.call(gridExtra::grid.arrange, append(dataPlots,
                                                  list(nrow = nrow)))
  
  # to inhibit plotting by grid.arrange
  grDevices::dev.off()
  
  return(list(plot = grid,
              width = width,
              height = height))
}

# plot per modification --------------------------------------------------------

.plot_per_modifications <- function(geneName,
                                    positions,
                                    mods,
                                    gff,
                                    fasta,
                                    modClasses){
  gff_sub <- gff[(is.na(S4Vectors::mcols(gff)$ID) & 
                    S4Vectors::mcols(gff)$Name == geneName) |
                   (!is.na(S4Vectors::mcols(gff)$ID) & 
                      S4Vectors::mcols(gff)$ID == geneName),]
  seq <- Rsamtools::getSeq(fasta,gff_sub)[[1]]
  
  # create description layer plot
  layer <- .create_layer_data(gff,gff_sub)
  
  # get data for different plot types
  posData <- .aggregate_pos_data(gff_sub,positions,modClasses)
  modData <- .aggregate_mod_data(mods,modClasses)
  modData <- lapply(names(modData), function(x){
    y <- modData[[x]]
    y$plotType <- x
    y
  })
  modData <- do.call(rbind,modData)
  
  plots <- lapply(1:nrow(modData), function(i){
    mod <- modData[i,]
    pos <- posData[[mod$plotType]]
    
    # focus on a modification
    localStart <- pos[pos$pos == mods$start,"localPos"] - RNAMOD_PLOT_FOCUS_WINDOW
    localEnd <- pos[pos$pos == mods$end,"localPos"] + RNAMOD_PLOT_FOCUS_WINDOW
    pos <- pos[pos$localPos > localStart & pos$localPos < localEnd,]
    
    # to inhibit plotting by grid.arrange
    grDevices::pdf(file = NULL)
    
    # get plots
    # layerPlot <- .get_gene_plot(geneName, layer, pos, seq)
    layerPlot <- list()
    plot <- .get_mod_plot(pos,mod)
    
    nrow <- length(layerPlot) + 1
    width <- 40 + nrow(pos) * RNAMOD_PLOT_POS_WIDTH
    height <- RNAMOD_PLOT_LAYER_HEIGHT + 
      RNAMOD_PLOT_DATA_HEIGHT
    
    # grid <- gridExtra::grid.arrange(plot, layerPlot, nrow = nrow)
    grid <- gridExtra::grid.arrange(plot, nrow = nrow)
    
    # to inhibit plotting by grid.arrange
    grDevices::dev.off()
    browser()
    return(list(plot = grid,
                width = width,
                height = height,
                names = rownames(mod)))
  })
  plot <- lapply(plots, "[[","plot")
  width <- lapply(plots, "[[","width")
  height <- lapply(plots, "[[","height")
  names <- lapply(plots, "[[","names")
  
  return(list(plot = plot,
              width = width,
              height = height,
              names = names))
}


# data aggregation function ----------------------------------------------------

# aggregates the position data for plotting of different plot types
.aggregate_pos_data <- function(gff,data, modClasses){
  plotTypes <- vapply(modClasses, getPlotType, character(1))
  plotTypes <- unique(plotTypes)
  data <- lapply(plotTypes, function(type, data){
    if(is.null(data[[type]])){
      stop("Positions data not available for type '",
           type,
           "'",
           .Call = FALSE)
    }
    res <- data[[type]]
    if( as.character(BiocGenerics::strand(gff)) == "-"){
      res <- res[order(res$pos, decreasing = TRUE),]
    }
    res$localPos <- 1:nrow(res)
    rownames(res) <- res$localPos
    return(res)
  }, data)
  names(data) <- plotTypes
  return(data)
}

# aggregates the modification data
.aggregate_mod_data <- function(data, modClasses){
  modTypes <- vapply(modClasses, getModType, character(1))
  modTypes <- unique(modTypes)
  data <- lapply(modTypes, function(type, data){
    if(nrow(data[data$RNAmod_type == type,]) == 0){
      return(NULL)
    }
    data[data$RNAmod_type == type,]
  }, data)
  names(data) <- vapply(modClasses, getPlotType, character(1))
  return(data)
}


# modification visualization ---------------------------------------------------

# returns a plot showing all modifications on one type of position data
.get_mod_plot <- function(pos,mods){
  requireNamespace("ggplot2", quietly = TRUE)
  
  break_FUN <- function(lim){
    if(length(lim[is.na(lim)])>0) return(lim)
    x <- 0
    y <- floor(abs(max(lim)/100))*100
    diff <- y - x
    n <- floor(diff / 100)
    as.numeric(c(x,unlist(lapply(1:n,function(z){100*z})),round(max(lim))))
  }
  # prepare mod data
  mods$localStart <- pos[pos$pos == mods$start,"localPos"]
  mods$localEnd <- pos[pos$pos == mods$end,"localPos"]
  mods$vStart <- max(pos[pos$localPos < mods$localStart+10 &
                           pos$localPos > mods$localStart-10,"mean"])*1.01
  modsPositions <- mods[mods$start == mods$end,]
  modsArea<- mods[mods$start != mods$end,]
  
  # initial plot setup
  plot <- ggplot(pos, aes_(x = ~localPos, y = ~mean)) +
    scale_x_continuous(name = "position of transcript [nt]") +
    scale_y_continuous(name = "mean(number of read ends)",
                       labels = scales::scientific,
                       limits = c(NA,max(pos$mean)*1.25)) +
    theme_bw()
  
  # plot modifications for an area
  if( nrow(modsArea) > 0){
    
  }
  
  # plot position data
  plot <- plot + geom_bar(stat = "identity")
  
  # plot modifications for a single position
  if( nrow(modsPositions) > 0){
    # setup p value text
    p_text <- unlist(lapply( modsPositions$RNAmod_p.value, function(p){
      if( p < 0.0001){
        return(paste0("< ",0.0001))
      } else {
        return(paste0(": ",round(p, digits = 4)))
      }
    }))
    # setup modification labels
    label <- paste0(modsPositions$RNAmod_type,
                    "\n",
                    rownames(modsPositions),
                    "\n \u03C3: ",
                    modsPositions$RNAmod_signal,
                    " (p ",
                    p_text,
                    ")")
    
    # plot modification marker
    plot <- plot + geom_segment(data = as.data.frame(modsPositions),
                                mapping = aes_(x = ~localStart,
                                               y = ~vStart,
                                               xend = ~localStart,
                                               yend = Inf,
                                               colour = ~RNAmod_type)) +
      scale_colour_brewer(name = "modification\ntype",
                          palette = "Set1") +
      annotate("label",
               x = modsPositions$localStart,
               y = Inf,
               label = label,
               vjust = "inward",
               size = 3)
  }
  
  return(plot)
}


# gene feature visualization ---------------------------------------------------

# returns a data.frame with the coordinates for plotting gene features
.create_layer_data <- function(gRanges,gRangeSelected){
  columns <- c("type","ID","Name","Parent")
  gRangeSelected <- gRangeSelected[,columns]
  gRanges <- gRanges[,columns]
  gRanges <- IRanges::subsetByOverlaps(gRanges, gRangeSelected)
  gRanges <- gRanges[as.character(S4Vectors::mcols(gRanges)$type) != "chromosome",]
  
  if( length(gRanges) == 0){
    return(NULL)
  }
  df <- S4Vectors::mcols(gRanges)[,c("type","ID","Name")]
  df$start <- BiocGenerics::start(gRanges)
  df$end <- BiocGenerics::end(gRanges)
  df$strand <- BiocGenerics::strand(gRanges)
  df <- as.data.frame(df)
  return(df)
}

# returns a plot showing all the gene annotation layers
.get_gene_plot <- function(geneName, 
                           layer,
                           pos,
                           seq){
  requireNamespace("ggplot2", quietly = TRUE)
  
  xlim <- c(min(pos$localPos),max(pos$localPos))
  ymin <- 0
  ymin2 <- 0
  browser()
  
  plot <- ggplot() +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    .get_read_line_colours() +
    scale_y_discrete(geneName) +
    labs(y = "none", x = "") 
  # plot chromosome
  plot <- plot + geom_segment(aes_(x = xlim[1],
                                   y = ~ymin2, 
                                   xend = xlim[2], 
                                   yend = ~ymin2, 
                                   colour = "chromosome"), 
                              size = 1)
  # # plot gene
  # plot <- plot + geom_segment(aes_(x = layer[layer$type == "gene","start"], 
  #                                  y = ~ymin2,
  #                                  xend = layer[layer$type == "gene","end"], 
  #                                  yend = ~ymin2, 
  #                                  colour = "gene"), 
  #                             size = 3)
  # # plot mRNA
  # plot <- plot + geom_segment(aes_(x = layer[layer$type == "mRNA","start"], 
  #                                  y = ~ymin2, 
  #                                  xend = layer[layer$type == "mRNA","end"], 
  #                                  yend = ~ymin2, 
  #                                  colour = "mRNA"), 
  #                             size = 1.5)
  # # plot CDS
  # plot <- plot + geom_segment(aes_(x = layer[layer$type == "CDS","start"], 
  #                                  y = ~ymin2, xend = layer[layer$type == "CDS","end"], 
  #                                  yend = ~ymin2, 
  #                                  colour = "CDS"), 
  #                             size = 2.0)
  # if( ("five_prime_UTR" %in% layer$type) ){
  #   # plot CDS
  #   plot <- plot + geom_segment(aes_(x = layer[layer$type == "five_prime_UTR","start"], 
  #                                    y = ~ymin2, 
  #                                    xend = layer[layer$type == "five_prime_UTR","end"], 
  #                                    yend = ~ymin2, 
  #                                    colour = "five_prime_UTR"), 
  #                               size = 1.0)
  # }
  # if( ("three_prime_UTR" %in% layer$type) ){
  #   # plot CDS
  #   plot <- plot + geom_segment(aes_(x = layer[layer$type == "three_prime_UTR","start"], 
  #                                    y = ~ymin2, 
  #                                    xend = layer[layer$type == "three_prime_UTR","end"], 
  #                                    yend = ~ymin2, 
  #                                    colour = "three_prime_UTR"), 
  #                               size = 1.0)
  # }
  # if( ("biological_region" %in% layer$type) ){
  #   uORFdata <- layer[layer$type == "biological_region",]
  #   for(i in seq_len(nrow(uORFdata))){
  #     # plot CDS
  #     plot <- plot + geom_segment(aes_(x = uORFdata[i,"start"], 
  #                                      y = ~ymin2, 
  #                                      xend = uORFdata[i,"end"], 
  #                                      yend = ~ymin2, 
  #                                      colour = "uORF"), 
  #                                 size = 2.0)
  #   }
  # }
  
  return(plot)
}

.get_read_line_colours <- function(){
  scale_colour_manual(values = c("Ribo" = "darkblue", 
                                 "RNA" = "darkgreen",  
                                 "gene" = "grey", 
                                 "mRNA" = "darkgrey", 
                                 "chromosome" = "black", 
                                 "CDS" = "orange",
                                 "intron" = "tan", 
                                 "uORF" = "red",
                                 "five_prime_UTR" = "violetred4",
                                 "three_prime_UTR" = "violetred"))
}
.get_read_fill_colours <- function(){
  scale_fill_manual(values = c("Ribo" = "darkblue", 
                               "RNA" = "darkgreen",  
                               "gene" = "grey", 
                               "mRNA" = "darkgrey", 
                               "chromosome" = "black", 
                               "CDS" = "orange", 
                               "intron" = "tan", 
                               "uORF" = "red",
                               "five_prime_UTR" = "violetred4",
                               "three_prime_UTR" = "violetred"))
}


# plot saving ------------------------------------------------------------------

#' save modification plots to file
#'
#' @param filetype "png" or "pdf"
#' @param folder vector of the folders to save to
#' @param names vector of the names to use to describe the file
#' @param plot vector of the plots to save
#' @param width vector of the width to use
#' @param height vector of the height to use
#'
#' @importFrom grDevices cairo_pdf png dev.off
#' @importFrom ggplot2 ggsave
#' @importFrom graphics plot
.save_gene_mod_plot <- function(filetype,
                                folder,
                                names,
                                plot,
                                width,
                                height){
  # check that for all plot all values are available
  x <- c(length(names),
         length(plot),
         length(width),
         length(height))
  if(!all( abs(x - mean(x)) < .Machine$double.eps ^ 0.5 ))
    stop("Not the same number of names, plots, height and width values ",
         "provided.",
         .call = TRUE)
  width <- unlist(width)
  height <- unlist(height)
  assertive::assert_all_are_non_empty_character(names)
  assertive::assert_all_are_whole_numbers(width)
  assertive::assert_all_are_whole_numbers(height)
  if(!assertive::is_a_number(getOption("RNAmod_dpi")))
    options("RNAmod_dpi",as.numeric(getOption("RNAmod_dpi")[[1]]))
  dpi <- getOption("RNAmod_dpi")
  
  fileNames <- paste0(folder,
                     "RNAmod_",
                     names,
                     ".",
                     filetype)
  for(i in 1:seq_along(length(fileNames))){
    if(filetype == "pdf"){
      if(assertive::r_has_cairo_capability() ){
        grDevices::cairo_pdf(fileNames[[i]],
                             width = (width[[i]]*RNAMOD_PLOT_MM_TO_INCH_F),
                             height = (height[[i]]*RNAMOD_PLOT_MM_TO_INCH_F))
        graphics::plot(plot[[i]])
        grDevices::dev.off()
      } else {
        ggplot2::ggsave(plot[[i]],
                        filename = fileNames[[i]],
                        units = "mm",
                        width = width[[i]],
                        height = height[[i]])
      }
    }
    if(filetype == "png"){
      if( assertive::r_has_cairo_capability() ){
        grDevices::png(fileNames[[i]],
                       units = "px",
                       width = (width[[i]]*
                                  RNAMOD_PLOT_MM_TO_INCH_F*
                                  dpi),
                       height = (height[[i]]*
                                   RNAMOD_PLOT_MM_TO_INCH_F*
                                   dpi),
                       type = "cairo-png")
        graphics::plot(plot[[i]])
        grDevices::dev.off()
      } else {
        ggplot2::ggsave(plot[[i]],
                        filename = fileNames[[i]],
                        units = "mm",
                        width = width[[i]],
                        height = height[[i]],
                        dpi = dpi)
      }
    }
    message("Saving mod plot for '",names[[i]],"'")
  }
}


