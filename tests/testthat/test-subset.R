
context("Subsetting")
test_that("Subsetting:",{
  data(msi,package = "RNAmodR")
  data(psd,package = "RNAmodR")
  # arguments
  expect_error(RNAmodR:::.norm_subset_args(),'argument "x" is missing')
  expect_error(RNAmodR:::.norm_subset_args(list()),'argument "x" is missing')
  actual <- RNAmodR:::.norm_subset_args(list(),msi)
  expect_type(actual,"list")
  expect_named(actual,c("name","type","flanking","rawData","perTranscript",
                        "sequenceData"))
  actual <- RNAmodR:::.norm_subset_args(list(),psd)
  expect_type(actual,"list")
  expect_named(actual,c("name","type","flanking","rawData","perTranscript",
                        "sequenceData"))
  # coords
  coord <- modifications(msi)[[1]]
  expect_error(RNAmodR:::.norm_coord(coord),'argument "type" is missing')
  expect_error(actual <- RNAmodR:::.norm_coord(DataFrame(),NA),"Something went wrong")
  actual <- RNAmodR:::.norm_coord(coord,NA)
  expect_s4_class(actual, "GRangesList")
  expect_equal(length(coord),length(actual))
  actual2 <- RNAmodR:::.norm_coord(coord,"I")
  expect_equal(actual,actual2)
  expect_equal(actual,RNAmodR:::.norm_coord(actual,NA))
  expect_error(RNAmodR:::.norm_coord(coord,"Ix"),"No modifications of type 'Ix'")
  coord$Parent <- NULL
  expect_error(RNAmodR:::.norm_coord(coord,NA),"Parent column must be present")
  # element names
  coord <- actual
  expect_error(RNAmodR:::.get_element_names(aggregateData(msi[[1]]),coord),'argument "name" is missing')
  expect_error(RNAmodR:::.get_element_names(aggregateData(msi[[1]]),coord,NA),'argument "type" is missing')
  actual <- RNAmodR:::.get_element_names(aggregateData(msi[[1]]),coord,NA,NA)
  expect_type(actual,"character")
  actual <- RNAmodR:::.get_element_names(aggregateData(msi[[1]]),coord,"2",NA)
  expect_equal(actual,"2")
  # subsetting
  actual <- subsetByCoord(sequenceData(msi[[1]]),coord)
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_equal(colnames(actual@unlistData),colnames(aggregate(sequenceData(msi[[1]]))@unlistData))
  actual <- subsetByCoord(msi[[1]],coord)
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_equal(colnames(actual@unlistData),mainScore(msi[[1]]))
  # labelling
  actual <- RNAmodR:::.label_SequenceData_by_GRangesList(sequenceData(msi[[1]]),coord)
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_equal(colnames(actual@unlistData),c(colnames(aggregate(sequenceData(msi[[1]]))@unlistData),"labels"))
  expect_type(actual@unlistData$labels,"logical")
  actual <- RNAmodR:::.label_Modifier_by_GRangesList(msi[[1]],coord)
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_equal(colnames(actual@unlistData),c(mainScore(msi[[1]]),"labels"))
  expect_type(actual@unlistData$labels,"logical")
  # SequenceDataList specific
  data(e5sd,package = "RNAmodR")
  sdl <- SequenceDataList(e5sd,psd)
  actual <- subsetByCoord(sdl,coord)
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  actual <- RNAmodR:::.label_SequenceDataList_by_GRangesList(sdl,coord)
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_type(actual@unlistData$labels,"logical")
})