
context("Visualization")
test_that("Visualization:",{
  data(msi, package = "RNAmodR")
  data(psd, package = "RNAmodR")
  data(e5sd, package = "RNAmodR")
  data(sds, package = "RNAmodR")
  data(sdl, package = "RNAmodR")
  getDefCoord <- function(){
    GRanges(seqnames = "chr1",
            ranges = IRanges::IRanges(start = 1510L,
                                      end = 1545L),
            strand = "+",
            Parent = "2")
  }
  getDefCoord2 <- function(){
    GRanges(seqnames = "chr2",
            ranges = IRanges::IRanges(start = 10L,
                                      end = 45L),
            strand = "+",
            Parent = "2")
  }
  coord <- getDefCoord()
  # internal functions SequenceData
  # .are_colours
  expect_error(RNAmodR:::.are_colours(),
               'argument "x" is missing, with no default')
  expect_false(RNAmodR:::.are_colours("abc"))
  expect_true(RNAmodR:::.are_colours("red"))
  expect_false(RNAmodR:::.are_colours("ABABAB"))
  expect_true(RNAmodR:::.are_colours("#ABABAB"))
  # .norm_viz_colour
  expect_error(RNAmodR:::.norm_viz_colour(),
               'argument "colour" is missing, with no default')
  expect_equal(RNAmodR:::.norm_viz_colour("red"),"red")
  expect_equal(RNAmodR:::.norm_viz_colour(c(x = "red"),"x"),
               c(x = "red"))
  expect_error(RNAmodR:::.norm_viz_colour(1,"x"),
               "'colour' must be a character vector and contain valid colours")
  expect_error(RNAmodR:::.norm_viz_colour("AB","x"),
               "'colour' must be a character vector and contain valid colours")
  expect_error(RNAmodR:::.norm_viz_colour(c("red","green"),"x"),
               "'colour' must be a named character vector parallel to 'type'")
  expect_error(RNAmodR:::.norm_viz_colour(c(x = "red",y = "green"),"x"),
               "'colour' must be a named character vector parallel to 'type'")
  expect_equal(RNAmodR:::.norm_viz_colour(c(x = "red",y = "green"),c("x","y")),
               c(x = "red",y = "green"))
  # .norm_viz_chromosome
  expect_error(RNAmodR:::.norm_viz_chromosome(),
               'argument "name" is missing, with no default')
  ranges <- ranges(psd)
  expect_error(RNAmodR:::.norm_viz_chromosome(ranges, "ab"),
               "Transcript name 'ab' not found in 'x'")
  expect_equal(RNAmodR:::.norm_viz_chromosome(ranges, "2"),
               "chr1")
  # .norm_coord_for_visualization
  coord2 <- coord
  expect_error(RNAmodR:::.norm_coord_for_visualization(ranges, c(coord,coord)),
               "'coord' must be contain a single range")
  coord2$Parent <- NULL
  expect_error(RNAmodR:::.norm_coord_for_visualization(ranges, coord2),
               "'coord' must contain a metadata column named 'Parent'")
  coord2 <- coord
  coord2$Parent <- "ab"
  expect_error(RNAmodR:::.norm_coord_for_visualization(ranges, coord2),
               "Transcript identifier 'ab' not found in data of 'x'")
  actual <- RNAmodR:::.norm_coord_for_visualization(ranges, coord)
  expect_s4_class(actual,"GRanges")
  expect_equal(coord,actual)
  # .norm_viz_windows.size
  expect_error(RNAmodR:::.norm_viz_windows.size(),
               'argument "window.size" is missing, with no default')
  expect_error(RNAmodR:::.norm_viz_windows.size(1),
               "'window.size' must be a single integer value")
  expect_equal(RNAmodR:::.norm_viz_windows.size(1L),1L)
  # .norm_viz_name
  expect_null(RNAmodR:::.norm_viz_name())
  expect_equal(RNAmodR:::.norm_viz_name("abc"),"abc")
  # .get_viz_from_to_coord
  coord <- GRanges(seqnames = "chr1",
                   ranges = IRanges::IRanges(start = 1510L,
                                             end = 1510L),
                   strand = "+",
                   Parent = "2")
  expect_error(RNAmodR:::.get_viz_from_to_coord(ranges, coord),
               'argument "window.size" is missing, with no default')
  actual <- RNAmodR:::.get_viz_from_to_coord(ranges, coord, 15L)
  expect_type(actual,"list")
  expect_named(actual,c("from","to"))
  expect_equal(actual,list(from = 1500L, to = 1525L))
  coord <- GRanges(seqnames = "chr1",
                   ranges = IRanges::IRanges(start = 1580L,
                                             end = 1580L),
                   strand = "+",
                   Parent = "2")
  actual <- RNAmodR:::.get_viz_from_to_coord(ranges, coord, 15L)
  expect_type(actual,"list")
  expect_named(actual,c("from","to"))
  expect_equal(actual,list(from = 1565L, to = 1595L))
  # .get_viz_from_to
  expect_error(RNAmodR:::.get_viz_from_to(ranges, "2"),
               'argument "from" is missing, with no default')
  expect_error(RNAmodR:::.get_viz_from_to(ranges, "", 1L),
               'argument "to" is missing, with no default')
  actual <- RNAmodR:::.get_viz_from_to(ranges, "2", 1501L, 1530L)
  expect_type(actual,"list")
  expect_named(actual,c("from","to"))
  expect_equal(actual,list(from = 1501L, to = 1530L))
  actual <- RNAmodR:::.get_viz_from_to(ranges, "2", 1L, 2000L)
  expect_type(actual,"list")
  expect_named(actual,c("from","to"))
  expect_equal(actual,list(from = 1500L, to = 1600L))
  # .stitch_chromosome
  seq <- RNAStringSet(c("AGCU","AGCU"))
  chr <- "chr2"
  gr <- GRanges(chr,1:50,"+",Parent = "2")
  gr2 <- GRanges(chr,100:150,"+",Parent = "2")
  ranges <- GRangesList("2" = gr, "3" = gr2)
  expect_error(RNAmodR:::.stitch_chromosome(seq, ranges, "chr3"),
               "No ranges with seqnames = ")
  expect_error(RNAmodR:::.stitch_chromosome(seq, ranges, chr),
               "No sequences for seqnames = ")
  seq <- RNAStringSet(c("2" = "AGCU","3" = "AGCU"))
  expect_error(RNAmodR:::.stitch_chromosome(seq, ranges, chr),
               "width\\(\\) or sequences and ranges does not match")
  gr <- GRanges(chr,1:4,"+",Parent = "2")
  gr2 <- GRanges(chr,100:103,"+",Parent = "2")
  ranges <- GRangesList("2" = gr, "3" = gr2)
  actual <- RNAmodR:::.stitch_chromosome(seq, ranges, chr)
  expect_s4_class(actual,"RNAStringSet")
  expect_equal(width(actual),103L)
  expect_equal(names(actual),chr)
  seq <- ModRNAStringSet(c("2" = "AGCU","3" = "AGCU"))
  actual <- RNAmodR:::.stitch_chromosome(seq, ranges, chr)
  expect_s4_class(actual,"ModRNAStringSet")
  expect_equal(width(actual),103L)
  expect_equal(names(actual),chr)
  # internal functions Modifier
  expect_false(RNAmodR:::.norm_show_argument())
  expect_true(RNAmodR:::.norm_show_argument(default = TRUE))
  expect_error(RNAmodR:::.norm_score_type(),"'type' is missing")
  expect_equal(RNAmodR:::.norm_score_type("x"),"x")
  expect_error(RNAmodR:::.norm_score_type("x",c("a","b")),
               "'type' was not found in data")
  expect_equal(RNAmodR:::.norm_score_type("a",c("a","b")),"a")
  expect_error(RNAmodR:::.norm_score_type(1,c("a","b")),
               "'type' must be a character vector")
  # arguments SequenceData
  expect_error(RNAmodR:::.norm_viz_args_SequenceData(),
               'argument "input" is missing, with no default')
  args <- RNAmodR:::.norm_viz_args_SequenceData(list())
  expect_type(args,"list")
  expect_named(args, c("alias","sequence.track.pars","annotation.track.pars",
                         "plot.pars"))
  expect_null(args[["alias"]])
  expect_type(args[["sequence.track.pars"]],"list")
  expect_type(args[["annotation.track.pars"]],"list")
  expect_type(args[["plot.pars"]],"list")
  # .get_viz_annotation_track
  actual <- RNAmodR:::.get_viz_annotation_track(psd, args)
  expect_s4_class(actual,"AnnotationTrack")
  # .get_viz_sequence_track
  expect_error(RNAmodR:::.get_viz_sequence_track(sequences(psd), ranges(psd), "chr2", args),
               "No ranges with seqnames = 'chr2' found")
  actual <- RNAmodR:::.get_viz_sequence_track(sequences(psd), ranges(psd), "chr1", args)
  expect_s4_class(actual,"RNASequenceTrack")
  actual <- RNAmodR:::.get_viz_sequence_track(sequences(psd), ranges(psd), "chr1", args)
  expect_s4_class(actual,"RNASequenceTrack")
  # .get_data_for_visualization
  expect_error(RNAmodR:::.get_data_for_visualization(),
               'argument "x" is missing, with no default')
  expect_error(RNAmodR:::.get_data_for_visualization(psd),
               'argument "name" is missing, with no default')
  actual <- RNAmodR:::.get_data_for_visualization(psd, "2")
  expect_s4_class(actual, "GRangesList")
  # getDataTrack
  actual <- getDataTrack(psd, name = "2")
  expect_type(actual,"list")
  expect_s4_class(actual[[1]],"DataTrack")
  # arguments Modifier
  expect_error(RNAmodR:::.norm_viz_args_Modifier(),
               'argument "input" is missing, with no default')
  expect_error(RNAmodR:::.norm_viz_args_Modifier(list(modified.seq = "1")),
               "'modified.seq' must be a single logical value.")
  expect_error(RNAmodR:::.norm_viz_args_Modifier(list(additional.mod = "1")),
               "'additional.mod' must be a GRanges or GRangesList object")
  actual <- RNAmodR:::.norm_viz_args_Modifier(list())
  expect_type(actual,"list")
  expect_named(actual, c("alias","sequence.track.pars","annotation.track.pars",
                         "plot.pars","modified.seq","additional.mod"))
  expect_null(actual[["alias"]])
  expect_type(actual[["sequence.track.pars"]],"list")
  expect_type(actual[["annotation.track.pars"]],"list")
  expect_type(actual[["plot.pars"]],"list")
  expect_s4_class(actual[["additional.mod"]],"GRanges")
  # .get_viz_sequence
  args <- 
    expect_error(RNAmodR:::.get_viz_sequence(),
                 'argument "x" is missing, with no default')
  expect_error(RNAmodR:::.get_viz_sequence(msi[[1]]),
               'argument "args" is missing, with no default')
  actual <- RNAmodR:::.get_viz_sequence(msi[[1]],
                                        RNAmodR:::.norm_viz_args_Modifier(list()))
  expect_s4_class(actual,"RNAStringSet")
  actual <- RNAmodR:::.get_viz_sequence(
    msi[[1]],
    RNAmodR:::.norm_viz_args_Modifier(list(modified.seq = TRUE)))
  expect_s4_class(actual,"ModRNAStringSet")
  actual <- separate(actual)
  expect_s4_class(actual,"GRanges")
  expect_equal(unique(actual$mod),"I")
  # SequenceData
  coord <- getDefCoord()
  actual <- visualizeData(psd, "2", from = 1500L, to = 1560L)
  expect_type(actual,"list")
  expect_s4_class(actual[[1]],"DataTrack")
  expect_s4_class(actual[[3]],"ImageMap")
  actual2 <- visualizeDataByCoord(psd, coord)
  expect_type(actual2,"list")
  expect_s4_class(actual2[[1]],"DataTrack")
  expect_s4_class(actual2[[3]],"ImageMap")
  expect_equal(actual,actual2)
  ##############################################################################
  # # SequenceDataSet
  # coord <- getDefCoord()
  # actual <- visualizeData(sds, "2", from = 1500L, to = 1560L)
  # expect_type(actual,"list")
  # expect_s4_class(actual[[1]],"DataTrack")
  # expect_s4_class(actual[[4]],"ImageMap")
  # actual2 <- visualizeDataByCoord(sds, coord)
  # expect_type(actual2,"list")
  # expect_s4_class(actual2[[1]],"DataTrack")
  # expect_s4_class(actual2[[4]],"ImageMap")
  # expect_equal(actual,actual2)
  # # SequenceDataList
  # coord <- getDefCoord()
  # actual <- visualizeData(sdl, "2", from = 1500L, to = 1560L)
  # expect_type(actual,"list")
  # expect_s4_class(actual[[1]],"DataTrack")
  # expect_s4_class(actual[[6]],"ImageMap")
  # actual2 <- visualizeDataByCoord(sdl, coord)
  # expect_type(actual2,"list")
  # expect_s4_class(actual2[[1]],"DataTrack")
  # expect_s4_class(actual2[[6]],"ImageMap")
  # expect_equal(actual,actual2)
  # Modifier
  coord <- getDefCoord2()
  actual <- visualizeData(msi[[1]], "2", from = 1L, to = 60L)
  expect_type(actual,"list")
  expect_s4_class(actual[[1]],"DataTrack")
  expect_s4_class(actual[[2]],"RNASequenceTrack")
  expect_s4_class(actual[[3]],"ImageMap")
  actual2 <- visualizeDataByCoord(msi[[1]], coord)
  expect_type(actual2,"list")
  expect_s4_class(actual2[[1]],"DataTrack")
  expect_s4_class(actual2[[2]],"RNASequenceTrack")
  expect_s4_class(actual2[[3]],"ImageMap")
  expect_equal(actual,actual2)
  # ModifierSet
  actual <- visualizeData(msi, "2", from = 1L, to = 60L)
  expect_type(actual,"list")
  expect_s4_class(actual[[1]],"DataTrack")
  expect_s4_class(actual[[4]],"RNASequenceTrack")
  expect_s4_class(actual[[5]],"ImageMap")
  actual2 <- visualizeDataByCoord(msi, coord)
  expect_type(actual2,"list")
  expect_s4_class(actual2[[1]],"DataTrack")
  expect_s4_class(actual2[[4]],"RNASequenceTrack")
  expect_s4_class(actual2[[5]],"ImageMap")
  expect_equal(actual,actual2)
  #############################################################################
  # data(e3sd, package = "RNAmodR")
  # data(esd, package = "RNAmodR")
  # data(ne3sd, package = "RNAmodR")
  # data(ne5sd, package = "RNAmodR")
  # data(csd, package = "RNAmodR")
  # data(pesd, package = "RNAmodR")
  # #
  # coord <- getDefCoord()
  # actual <- visualizeData(e5sd, "2", from = 1500L, to = 1560L)
  # expect_type(actual,"list")
  # expect_s4_class(actual[[1]],"DataTrack")
  # expect_s4_class(actual[[3]],"ImageMap")
  # actual2 <- visualizeDataByCoord(e5sd, coord)
  # expect_type(actual2,"list")
  # expect_s4_class(actual2[[1]],"DataTrack")
  # expect_s4_class(actual2[[3]],"ImageMap")
  # expect_equal(actual,actual2)
  # #
  # actual <- visualizeData(e3sd, "2", from = 1500L, to = 1560L)
  # expect_type(actual,"list")
  # expect_s4_class(actual[[1]],"DataTrack")
  # expect_s4_class(actual[[3]],"ImageMap")
  # actual2 <- visualizeDataByCoord(e3sd, coord)
  # expect_type(actual2,"list")
  # expect_s4_class(actual2[[1]],"DataTrack")
  # expect_s4_class(actual2[[3]],"ImageMap")
  # expect_equal(actual,actual2)
  # #
  # actual <- visualizeData(esd, "2", from = 1500L, to = 1560L)
  # expect_type(actual,"list")
  # expect_s4_class(actual[[1]],"DataTrack")
  # expect_s4_class(actual[[3]],"ImageMap")
  # actual2 <- visualizeDataByCoord(esd, coord)
  # expect_type(actual2,"list")
  # expect_s4_class(actual2[[1]],"DataTrack")
  # expect_s4_class(actual2[[3]],"ImageMap")
  # expect_equal(actual,actual2)
  # #
  # actual <- visualizeData(ne3sd, "2", from = 1500L, to = 1560L)
  # expect_type(actual,"list")
  # expect_s4_class(actual[[1]],"DataTrack")
  # expect_s4_class(actual[[4]],"ImageMap")
  # actual2 <- visualizeDataByCoord(ne3sd, coord)
  # expect_type(actual2,"list")
  # expect_s4_class(actual2[[1]],"DataTrack")
  # expect_s4_class(actual2[[4]],"ImageMap")
  # expect_equal(actual,actual2)
  # #
  # actual <- visualizeData(ne5sd, "2", from = 1500L, to = 1560L)
  # expect_type(actual,"list")
  # expect_s4_class(actual[[1]],"DataTrack")
  # expect_s4_class(actual[[4]],"ImageMap")
  # actual2 <- visualizeDataByCoord(ne5sd, coord)
  # expect_type(actual2,"list")
  # expect_s4_class(actual2[[1]],"DataTrack")
  # expect_s4_class(actual2[[4]],"ImageMap")
  # expect_equal(actual,actual2)
  # #
  # actual <- visualizeData(csd, "2", from = 1500L, to = 1560L)
  # expect_type(actual,"list")
  # expect_s4_class(actual[[1]],"DataTrack")
  # expect_s4_class(actual[[3]],"ImageMap")
  # actual2 <- visualizeDataByCoord(csd, coord)
  # expect_type(actual2,"list")
  # expect_s4_class(actual2[[1]],"DataTrack")
  # expect_s4_class(actual2[[3]],"ImageMap")
  # expect_equal(actual,actual2)
  # #
  # actual <- visualizeData(pesd, "2", from = 1500L, to = 1560L)
  # expect_type(actual,"list")
  # expect_s4_class(actual[[1]],"DataTrack")
  # expect_s4_class(actual[[3]],"ImageMap")
  # actual2 <- visualizeDataByCoord(pesd, coord)
  # expect_type(actual2,"list")
  # expect_s4_class(actual2[[1]],"DataTrack")
  # expect_s4_class(actual2[[3]],"ImageMap")
  # expect_equal(actual,actual2)
})
