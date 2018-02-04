#' @include RNAmodR.R
NULL

RNAMODR_PLOT_POS_WIDTH <- 1
RNAMODR_PLOT_LAYER_HEIGHT <- 50
RNAMODR_PLOT_SEQ_SIZE <- 1
RNAMODR_PLOT_DATA_HEIGHT <- 100
RNAMODR_PLOT_FOCUS_WINDOW <- 50
RNAMODR_PLOT_MM_TO_INCH_F <- 0.03937008

#' @rdname getModPlot
#' @aliases getModPlot saveModPlots saveModPlotFromSummarizedExperiment
#' 
#' @title Visualization of modifications
#' 
#' @description 
#' \code{getModPlot} and \code{saveModPlot} plot the results by combining
#' position data and modifications with transcript layout data.
#' 
#' @param .Object a RNAmod object. 
#' @param number an experiemnt numnber. 
#' @param se a SummarizedExperiment containg the experimental data. 
#' @param gff a GRanges object with the genomic annotation data
#' @param fasta a FaFile object with the genomic sequences
#' @param modifications name of modification to be used for analysis. 
#' @param gene a single gene name used by \code{getModPlot}
#' @param genes gene names used by \code{saveModPlots}
#' @param folder folder for saving the output used by \code{saveModPlots}
#' @param focus optional: logical, whether to plot only the region next to the 
#' modification
#' @param filetype file type used for \code{saveModPlots}
#'
#' @return 
#' \code{getModPlot}: a \code{TableGrob} of layer plot and any mod plot
#' @export
#' 
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#'
#' @examples
#' \donttest{
#' getModPlot(mod,1,"m7G","RDN18-1")
#' saveModPlot(mod,1,"m7G",c("RDN18-1","tS\\(CGA\\)P"))
#' }
setMethod(
  f = "getModPlot", 
  signature = signature(.Object = "RNAmodR",
                        number = "numeric",
                        modifications = "character",
                        gene = "character"),
  definition = function(.Object,
                        number,
                        modifications,
                        gene,
                        focus){
    # get se, gff and fasta
    se <- getSummarizedExperiment(.Object,number,modifications)
    gff <- .Object@.dataGFF
    fasta <- .Object@.dataFasta
    # call plotting function
    getModPlotFromSummarizedExperiment(se = se,
                                       gff = gff,
                                       fasta = fasta,
                                       modifications = modifications,
                                       gene = gene,
                                       focus = focus)
  }
)
#' @rdname getModPlot
#' @export
setMethod(
  f = "getModPlotFromSummarizedExperiment", 
  signature = signature(se = "SummarizedExperiment",
                        gff = "GRanges",
                        fasta = "FaFile",
                        modifications = "character",
                        gene = "character"),
  definition = function(se,
                        gff,
                        fasta,
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
                paste0(setdiff(gene,genesAvail), collapse = "', '"),"'",
                call. = FALSE)
      } else {
        stop("No information available for any genes given: '",
             paste0(gene, collapse = "', '"),"'",
             call. = FALSE)
      }
    }
    genesAvail <- unlist(genesAvail)
    assertive::assert_is_a_non_empty_string(genesAvail)
    
    positions <- seGetPositions(se,genesAvail)[[1]]
    mods <- seGetModifications(se,genesAvail,modifications)[[1]]
    
    modClasses <- .load_mod_classes(modifications)
    if(!focus){
      plot <- .plot_gene_with_modifications(genesAvail,
                                            positions,
                                            mods,
                                            gff,
                                            fasta,
                                            modClasses)
    } else {
      plot <- .plot_per_modifications(genesAvail,
                                      positions,
                                      mods,
                                      gff,
                                      fasta,
                                      modClasses)
    }
    return(plot$plot[[1]])
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
         call. = FALSE)
}


#' @rdname getModPlot
#'
#' @return 
#' \code{saveModPlot}: TRUE if successful
#' 
#' @export
setMethod(
  f = "saveModPlot", 
  signature = signature(.Object = "RNAmodR",
                        number = "numeric",
                        modifications = "character",
                        genes = "character"),
  definition = function(.Object,
                        number,
                        modifications,
                        genes,
                        focus,
                        filetype){
    # get se, gff and fasta
    se <- getSummarizedExperiment(.Object,
                                  number,
                                  modifications)
    gff <- .Object@.dataGFF
    fasta <- .Object@.dataFasta
    # create folder
    folder <- paste0(getOutputFolder(.Object),
                     "ModPlots/")
    return(saveModPlotFromSummarizedExperiment(se = se,
                                               gff = gff,
                                               fasta = fasta,
                                               modifications = modifications,
                                               genes = genes,
                                               folder = folder,
                                               focus = focus,
                                               filetype = filetype))
  }
)
#' @rdname getModPlot
#' @export
setMethod(
  f = "saveModPlot", 
  signature = signature(.Object = "RNAmodR",
                        number = "numeric",
                        modifications = "character",
                        genes = "missing"),
  definition = function(.Object,
                        number,
                        modifications,
                        focus,
                        filetype){
    # get se, gff and fasta
    se <- getSummarizedExperiment(.Object,
                                  number,
                                  modifications)
    gff <- .Object@.dataGFF
    fasta <- .Object@.dataFasta
    # get transcript names with modifications
    assays <- SummarizedExperiment::assays(se)[modifications]
    l <- lapply(seq_along(assays), function(i){
      assay <- assays[[i]]
      names(assay[assay > 0,])
    })
    l <- unlist(l)
    if(length(l) == 0){
      stop("No modifications for any genes found.",
           call. = FALSE)
    }
    transcripts <- unique(l)
    # create folder
    folder <- paste0(getOutputFolder(.Object),
                     "ModPlots/")
    # call plotting function
    return(saveModPlotFromSummarizedExperiment(se = se,
                                               gff = gff,
                                               fasta = fasta,
                                               modifications = modifications,
                                               genes = transcripts,
                                               folder = folder,
                                               focus = focus,
                                               filetype = filetype))
  }
)
#' @rdname getModPlot
#' @export
setMethod(
  f = "saveModPlotFromSummarizedExperiment", 
  signature = signature(se = "SummarizedExperiment",
                        gff = "GRanges",
                        fasta = "FaFile",
                        modifications = "character",
                        genes = "character",
                        folder = "character"),
  definition = function(se,
                        gff,
                        fasta,
                        modifications,
                        genes,
                        folder,
                        focus,
                        filetype){
    # Check input
    .check_plot_input(se,
                      modifications,
                      genes,
                      focus,
                      filetype)
    
    # create folder
    folder <- paste0(folder,
                     # unique(SummarizedExperiment::colData(se)$SampleName),
                     unique(SummarizedExperiment::colData(se)$name),
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
    
    # browser()
    positions <- seGetPositions(se,genesAvail)
    mods <- seGetModifications(se,genesAvail,modifications)
    
    modClasses <- .load_mod_classes(modifications)
    
    for(i in seq_along(genesAvail)){
      if(length(positions[[i]][!vapply(positions[[i]],is.null,logical(1))]) == 0){
        message("Skipping plot for '",
                genesAvail[[i]],
                "'. No position data available.")
        next
      }
      if(!focus){
        plot <- .plot_gene_with_modifications(genesAvail[[i]],
                                              positions[[i]],
                                              mods[[i]],
                                              gff,
                                              fasta,
                                              modClasses)
        names <- paste0(paste(modifications, collapse = "-"),
                        "_",
                        genesAvail[[i]])
      } else {
        if(nrow(mods[[i]]) == 0){
          message("Skipping plot for '",
                  genesAvail[[i]],
                  "'. No modification data available.")
          next
        }
        plot <- .plot_per_modifications(genesAvail[[i]],
                                        positions[[i]],
                                        mods[[i]],
                                        gff,
                                        fasta,
                                        modClasses)
        names <- paste0(genesAvail[[i]],
                        "_",
                        plot$names)
      }
      .save_gene_mod_plot(filetype,
                          folder,
                          names,
                          plot$plot,
                          plot$width,
                          plot$height)
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
  gff_sub <- .subset_gff_for_unique_transcript(gff, 
                                               geneName)
  seq <- Rsamtools::getSeq(fasta,gff_sub)[[1]]
  letters <- .get_transcript_sequence(gff,
                                      geneName,
                                      seq)
  # create description layer plot
  layer <- .create_layer_data(gff,gff_sub)
  
  # get data for different plot types
  posData <- .aggregate_pos_data(gff_sub,
                                 positions,
                                 letters,
                                 modClasses)
  mods <- mods[order(mods$start),]
  modData <- .aggregate_mod_data(mods,modClasses)
  
  # to inhibit plotting by grid.arrange
  grDevices::pdf(file = NULL)
  
  dataPlots <- lapply(names(posData$data), function(type){
    return(.get_mod_plot(pos = posData$data[[type]],
                         mods = modData[[type]],
                         label = posData$label[[type]],
                         format = posData$format[[type]])
    )
  })
  #layerPlot <- .get_gene_plot(geneName, layer, posData$data[[1]])
  layerPlot <- list()
  
  nrow <- length(layerPlot) + length(dataPlots)
  width <- 80 + nrow(posData$data[[1]]) * RNAMODR_PLOT_POS_WIDTH
  height <- length(layerPlot) * RNAMODR_PLOT_LAYER_HEIGHT + 
    length(dataPlots) * RNAMODR_PLOT_DATA_HEIGHT
  
  # grid <- do.call(gridExtra::grid.arrange, append(dataPlots,
  #                                                 layerPlot,
  #                                                 list(nrow = nrow)))
  grid <- do.call(gridExtra::grid.arrange, append(dataPlots,
                                                  list(nrow = nrow)))
  
  # to inhibit plotting by grid.arrange
  grDevices::dev.off()
  
  return(list(plot = list(grid),
              width = list(width),
              height = list(height)))
}

# plot per modification --------------------------------------------------------

#' @importFrom gridExtra grid.arrange
#' 
.plot_per_modifications <- function(geneName,
                                    positions,
                                    mods,
                                    gff,
                                    fasta,
                                    modClasses){
  
  gff_sub <- .subset_gff_for_unique_transcript(gff, 
                                               geneName)
  seq <- Rsamtools::getSeq(fasta,gff_sub)[[1]]
  letters <- .get_transcript_sequence(gff,
                                      geneName,
                                      seq)
  # create description layer plot
  layer <- .create_layer_data(gff,gff_sub)
  
  # get data for different plot types
  posData <- .aggregate_pos_data(gff_sub,
                                 positions,
                                 letters,
                                 modClasses)
  mods <- mods[order(mods$start),]
  modData <- .aggregate_mod_data(mods,modClasses)
  modData <- lapply(names(modData), function(x){
    y <- modData[[x]]
    y$plotType <- x
    y
  })
  modData <- do.call(rbind,modData)
  
  plots <- lapply(1:nrow(modData), function(i){
    # subset tp temporary data
    mod <- modData[i,]
    pos <- posData$data[[mod$plotType]]
    # focus on a modification
    localStart <- as.numeric(pos[as.numeric(pos$pos) %in% mods$start,"pos"]) - 
      RNAMODR_PLOT_FOCUS_WINDOW
    localEnd <- as.numeric(pos[as.numeric(pos$pos) %in% mods$end,"pos"]) + 
      RNAMODR_PLOT_FOCUS_WINDOW
    pos <- pos[as.numeric(pos$pos) > localStart &
                 as.numeric(pos$pos) < localEnd,]
    # to inhibit plotting by grid.arrange
    grDevices::pdf(file = NULL)
    # get plots
    # layerPlot <- .get_gene_plot(geneName, layer, pos)
    layerPlot <- list()
    plot <- .get_mod_plot(pos = pos,
                          mods = mod,
                          label = posData$label[[mod$plotType]],
                          format = posData$format[[mod$plotType]])
    # get dimensions of the plot
    nrow <- length(layerPlot) + 1
    width <- 80 + nrow(pos) * RNAMODR_PLOT_POS_WIDTH
    height <- RNAMODR_PLOT_LAYER_HEIGHT + 
      RNAMODR_PLOT_DATA_HEIGHT
    # create grid
    # grid <- gridExtra::grid.arrange(plot, layerPlot, nrow = nrow)
    grid <- gridExtra::grid.arrange(plot, nrow = nrow)
    # to inhibit plotting by grid.arrange
    grDevices::dev.off()
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
.aggregate_pos_data <- function(gff,
                                data, 
                                letters,
                                modClasses){
  plotTypes <- vapply(modClasses, getAnalysisType, character(1))
  plotTypes <- unique(plotTypes)
  data <- lapply(plotTypes, function(type, data){
    if(is.null(data[[type]])){
      stop("No positions data not available for '",
           S4Vectors::mcols(gff)$ID,
           "' in type '",
           type,
           "'",
           call. = FALSE)
    }
    res <- data[[type]]
    # short time fix
    res$pos <- as.numeric(res$pos)[order(as.numeric(res$pos))]
    # add letters to position data
    res$letters <- letters
    return(res)
  }, data)
  # load the analysis classes
  analysisClasses <- .load_analysis_classes(plotTypes)
  label <- lapply(analysisClasses, getDataLabel)
  format <- lapply(analysisClasses, getDataFormat)
  names(data) <- plotTypes
  names(label) <- plotTypes
  names(format) <- plotTypes
  return(list(data = data,
              label = label,
              format = format))
}

# aggregates the modification data
.aggregate_mod_data <- function(data, modClasses){
  modTypes <- vapply(modClasses, getModType, character(1))
  modTypes <- unique(modTypes)
  # get modifications which are requested based on type
  data <- lapply(modTypes, function(type){
    if(nrow(data[data$RNAmodR_type == type,]) == 0){
      return(NULL)
    }
    data[data$RNAmodR_type == type,]
  })
  # set name of plot type
  names(data) <- vapply(modClasses, getAnalysisType, character(1))
  # aggregated data based on plot type
  data <- stats::setNames(lapply(unique(names(data)), function(name){
    x <- do.call(rbind, data[names(data) == name])
    if(is.null(x)) return(NULL)
    return(x[order(x$start,x$end),])
  })
  ,unique(names(data)))
  return(data)
}


# modification visualization ---------------------------------------------------

# returns a plot showing all modifications on one type of position data
.get_mod_plot <- function(pos,
                          mods,
                          label,
                          format){
  requireNamespace("ggplot2", quietly = TRUE)
  # browser()
  break_FUN <- function(lim){
    if(length(lim[is.na(lim)])>0) return(lim)
    x <- 0
    y <- floor(abs(max(lim)/100))*100
    diff <- y - x
    n <- floor(diff / 100)
    as.numeric(c(x,unlist(lapply(1:n,function(z){100*z})),round(max(lim))))
  }
  
  
  # initial plot setup
  plot <- ggplot(pos, aes_(x = ~as.numeric(pos), y = ~mean, label = ~letters)) +
    scale_x_continuous(name = "position of transcript [nt]",
                       expand = c(0,10)) +
    scale_y_continuous(name = label,
                       labels = format,
                       limits = c(NA,max(pos$mean)*1.3)) +
    theme_bw()
  
  # plot position data
  plot <- plot + geom_bar(stat = "identity") +
    geom_text(mapping = aes_(y = 0),
              vjust = 1.5,
              size = RNAMODR_PLOT_SEQ_SIZE)
  
  # if no mod data available
  if(is.null(mods)) return(plot)
  
  # prepare mod data
  mods$localStart <- as.numeric(pos[as.numeric(pos$pos) %in% mods$start,"pos"])
  mods$localEnd <- as.numeric(pos[as.numeric(pos$pos) %in% mods$end,"pos"])
  mods$y <- unlist(lapply(mods$localStart, function(x){
    pos[pos$pos == x,"mean"]*1.02
  }))
  mods$yend <- unlist(lapply(mods$localStart, function(x){
    max(pos$mean)
  }))
  modsPositions <- mods[mods$start == mods$end,]
  modsArea<- mods[mods$start != mods$end,]
  
  # plot modifications for an area
  if( nrow(modsArea) > 0){
    
  }
  
  # if no mod data available
  if( nrow(modsPositions) == 0) return(plot)
    
  # plot modifications for single positions each
  # setup p value text
  p_text <- unlist(lapply(as.numeric(modsPositions$RNAmodR_p.value), 
                          function(p){
    if( p < 0.0001){
      return(paste0("< ",0.0001))
    } else {
      return(paste0(": ",round(p, digits = 4)))
    }
  }))
  # setup modification labels
  modsPositions$label <- paste0(modsPositions$RNAmodR_type,
                                "\n",
                                rownames(modsPositions),
                                "\n \u03C3: ",
                                modsPositions$RNAmodR_signal,
                                " (p ",
                                p_text,
                                ")")
  # plot modification marker
  plot <- plot + geom_segment(data = as.data.frame(modsPositions),
                              mapping = aes_(x = ~localStart,
                                             y = ~y,
                                             xend = ~localStart,
                                             yend = ~yend,
                                             colour = ~RNAmodR_type),
                              inherit.aes = FALSE) +
    scale_colour_brewer(name = "Modification\ntype",
                        palette = "Set1") +
    ggrepel::geom_label_repel(data = as.data.frame(modsPositions),
                              mapping = aes_(x = ~localStart,
                                             y = ~yend,
                                             label = modsPositions$label),
                              segment.color = 'grey50',
                              min.segment.length = 0,
                              box.padding = 0.3,
                              point.padding = 0,
                              direction = "x",
                              ylim = c(max(pos$mean)*1.07,
                                       NA),
                              size = 3)
  
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
                           pos){
  requireNamespace("ggplot2", quietly = TRUE)
  
  xlim <- c(min(pos$pos),max(pos$pos))
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
         call. = TRUE)
  width <- unlist(width)
  height <- unlist(height)
  assertive::assert_all_are_non_empty_character(names)
  assertive::assert_all_are_whole_numbers(width)
  assertive::assert_all_are_whole_numbers(height)
  if(!assertive::is_a_number(getOption("RNAmodR_dpi")))
    options("RNAmodR_dpi",as.numeric(getOption("RNAmodR_dpi")[[1]]))
  dpi <- getOption("RNAmodR_dpi")
  
  fileNames <- paste0(folder,
                     "RNAmodR_",
                     names,
                     ".",
                     filetype)
  for(i in 1:seq_along(length(fileNames))){
    if( assertive::r_has_cairo_capability() & getOption("RNAmodR_use_cairo") ){
      if(filetype == "pdf"){
        grDevices::cairo_pdf(fileNames[[i]],
                             width = (width[[i]]*RNAMODR_PLOT_MM_TO_INCH_F),
                             height = (height[[i]]*RNAMODR_PLOT_MM_TO_INCH_F))
        graphics::plot(plot[[i]])
        grDevices::dev.off()
      }
      if(filetype == "png"){
        grDevices::png(fileNames[[i]],
                       units = "in",
                       width = (width[[i]]*
                                  RNAMODR_PLOT_MM_TO_INCH_F),
                       height = (height[[i]]*
                                   RNAMODR_PLOT_MM_TO_INCH_F),
                       type = "cairo-png",
                       res = dpi)
        graphics::plot(plot[[i]])
        grDevices::dev.off()
      }
    } else {
      ggplot2::ggsave(plot[[i]],
                      filename = fileNames[[i]],
                      units = "mm",
                      width = width[[i]],
                      height = height[[i]],
                      dpi = dpi,
                      limitsize = FALSE)
    }
    
    message("Saving mod plot for '",names[[i]],"'")
  }
}


