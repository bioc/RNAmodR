
context("Subsetting SequenceData")
test_that("Subsetting SequenceData:",{
  data(msi,package = "RNAmodR")
  data(psd,package = "RNAmodR")
  # arguments
  expect_error(RNAmodR:::.norm_subset_args(),'argument "x" is missing')
  expect_error(RNAmodR:::.norm_subset_args(list()),'argument "x" is missing')
  actual <- RNAmodR:::.norm_subset_args(list(),msi)
  expect_error(RNAmodR:::.norm_subset_args(list(name = 1),msi),
               "'name' must be a character with a width > 0L")
  expect_error(RNAmodR:::.norm_subset_args(list(type = 1),msi),
               "'type' must be a character with a width > 0L")
  expect_error(RNAmodR:::.norm_subset_args(list(type = "meep"),msi),
               "'type' must be one or more elements of shortName(ModRNAString")
  expect_error(RNAmodR:::.norm_subset_args(list(flanking = 1),msi),
        "'flanking' must be a single integer value equal or higher than 0L")
  expect_error(RNAmodR:::.norm_subset_args(list(rawData = 1),msi),
               "'rawData' must be a single logical value")
  expect_error(RNAmodR:::.norm_subset_args(list(perTranscript = 1),msi),
               "'perTranscript' must be a single logical value")
  expect_error(RNAmodR:::.norm_subset_args(list(sequenceData = 1),msi),
               "'sequenceData' must be a single logical value")
  expect_type(actual,"list")
  expect_named(actual,c("name","type","flanking","rawData","perTranscript",
                        "sequenceData"))
  actual <- RNAmodR:::.norm_subset_args(list(),psd)
  expect_type(actual,"list")
  expect_named(actual,c("name","type","flanking","rawData","perTranscript",
                        "sequenceData"))
  expect_equal(RNAmodR:::.norm_subset_args(list(name = "abc"),msi)$name,"abc")
  expect_equal(RNAmodR:::.norm_subset_args(list(type = "I"),msi)$type,"I")
  expect_equal(RNAmodR:::.norm_subset_args(list(flanking = 1L),msi)$flanking,1L)
  expect_equal(RNAmodR:::.norm_subset_args(list(rawData = TRUE),msi)$rawData,
               TRUE)
  expect_equal(RNAmodR:::.norm_subset_args(list(rawData = FALSE),msi)$rawData,
               FALSE)
  expect_equal(RNAmodR:::.norm_subset_args(list(perTranscript = TRUE),msi)$perTranscript,
               TRUE)
  expect_equal(RNAmodR:::.norm_subset_args(list(perTranscript = FALSE),msi)$perTranscript,
               FALSE)
  expect_equal(RNAmodR:::.norm_subset_args(list(sequenceData = TRUE),msi)$sequenceData,
               TRUE)
  expect_equal(RNAmodR:::.norm_subset_args(list(sequenceData = FALSE),msi)$sequenceData,
               FALSE)
  # coords
  coord <- modifications(msi)[[1]]
  expect_error(RNAmodR:::.norm_coord(coord),'argument "type" is missing')
  expect_error(actual <- RNAmodR:::.norm_coord(DataFrame(),NA),
               "Something went wrong")
  coord2 <- coord
  width(coord2) <- 2
  expect_error(RNAmodR:::.norm_coord(coord2,NA),
               "Elements of type 'RNAMOD' with a width != 1L are not supported")
  coord2 <- coord
  coord2$Parent <- c("a","b","c","d","e","f")
  actual <- RNAmodR:::.norm_coord(coord,NA)
  expect_s4_class(actual, "GRangesList")
  expect_equal(length(coord),length(actual))
  actual2 <- RNAmodR:::.norm_coord(coord,"I")
  expect_equal(actual,actual2)
  expect_equal(actual,RNAmodR:::.norm_coord(actual,NA))
  expect_error(RNAmodR:::.norm_coord(coord,"Ix"),
               "No modifications of type 'Ix'")
  coord$Parent <- NULL
  expect_error(RNAmodR:::.norm_coord(coord,NA),"Parent column must be present")
  # subsetting
  actual <- subsetByCoord(sequenceData(msi[[1]]),coord)
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_equal(actual,
               subsetByCoord(sequenceData(msi[[1]]),unlist(coord)))
  expect_equal(colnames(actual@unlistData),
               colnames(aggregate(sequenceData(msi[[1]]))@unlistData))
  expect_error(subsetByCoord(sequenceData(msi[[1]]),coord2))
  # labelling
  actual <- RNAmodR:::.label_SequenceData_by_GRangesList(sequenceData(msi[[1]]),
                                                         coord)
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_equal(colnames(actual@unlistData),
               c(colnames(aggregate(sequenceData(msi[[1]]))@unlistData),
                 "labels"))
  expect_type(actual@unlistData$labels,"logical")
  # SequenceDataList specific
  data(e5sd,package = "RNAmodR")
  sdl <- SequenceDataList(e5sd,psd)
  actual <- subsetByCoord(sdl,coord)
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_equal(actual,subsetByCoord(sdl,unlist(coord)))
  actual <- RNAmodR:::.label_SequenceDataList_by_GRangesList(sdl,coord)
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_type(actual@unlistData$labels,"logical")
  # raw data
  expect_error(subsetByCoord(sequenceData(msi[[1]]),coord, rawData = 1),
               "'rawData' must be a single logical value")
  expect_error(subsetByCoord(sequenceData(msi[[1]]),coord, rawData = c(TRUE,TRUE)),
               "'rawData' must be a single logical value")
  actual <- subsetByCoord(sequenceData(msi[[1]]),coord, rawData = TRUE)
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_equal(unique(ncol(actual)),15)
  actual <- subsetByCoord(sdl,coord, rawData = TRUE)
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_equal(unique(ncol(actual)),18)
  expect_equivalent(cbind(subsetByCoord(sdl[[1]],coord, rawData = TRUE),
                     subsetByCoord(sdl[[2]],coord, rawData = TRUE)),
                    subsetByCoord(sdl,coord, rawData = TRUE))
  actual <- RNAmodR:::.label_SequenceData_by_GRangesList(sdl[[1]],coord)
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_type(actual@unlistData$labels,"logical")
  actual <- RNAmodR:::.label_SequenceDataList_by_GRangesList(sdl,coord)
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_type(actual@unlistData$labels,"logical")
})

context("Subsetting Modifier/ModifierSet")
test_that("Subsetting Modifier/ModifierSet:",{
  data(msi,package = "RNAmodR")
  # element names
  coord <- actual
  expect_error(RNAmodR:::.get_element_names(aggregateData(msi[[1]]),coord),
               'argument "name" is missing')
  expect_error(RNAmodR:::.get_element_names(aggregateData(msi[[1]]),coord,NA),
               'argument "type" is missing')
  actual <- RNAmodR:::.get_element_names(aggregateData(msi[[1]]),coord,NA,NA)
  expect_type(actual,"character")
  actual <- RNAmodR:::.get_element_names(aggregateData(msi[[1]]),coord,"2",NA)
  expect_equal(actual,"2")
  # subsetting
  actual <- subsetByCoord(msi[[1]],coord)
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_equal(colnames(actual@unlistData),mainScore(msi[[1]]))
  expect_equal(actual,subsetByCoord(msi[[1]],unlist(coord)))
  expect_equal(subsetByCoord(sequenceData(msi[[1]]),coord),
               subsetByCoord(msi[[1]],coord,sequenceData = TRUE))
  actual2 <- subsetByCoord(msi,coord)
  expect_named(actual2,names(msi))
  expect_equal(actual2[[1]],actual)
  expect_equal(actual2,subsetByCoord(msi,unlist(coord)))
  actual3 <- subsetByCoord(msi,coord,sequenceData = TRUE)
  expect_type(actual3,"list")
  actual4 <- lapply(msi,subsetByCoord,coord,sequenceData = TRUE)
  expect_equal(actual3,actual4)
  # labelling
  actual <- RNAmodR:::.label_Modifier_by_GRangesList(msi[[1]],coord)
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_equal(colnames(actual@unlistData),c(mainScore(msi[[1]]),"labels"))
  expect_type(actual@unlistData$labels,"logical")
  
})
