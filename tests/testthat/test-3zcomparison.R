
context("Comparing data")
test_that("Comparing data:",{
  data(msi, package = "RNAmodR")
  getDefCoord <- function(){
    GRanges(seqnames = "chr2",
            ranges = IRanges::IRanges(start = c(10,20,30,40,50), width = 1L),
            strand = "+",
            Parent = c("2","2","2","2","2"))
  }
  coord <- getDefCoord()
  # arguments norm
  expect_error(RNAmodR:::.norm_alias(),
               'argument "input" is missing, with no default')
  actual <- RNAmodR:::.norm_alias(list())
  expect_type(actual,"list")
  expect_named(actual,"alias")
  library(RNAmodR.Data)
  library(Rsamtools)
  library(GenomicFeatures)
  annotation <- GFF3File(RNAmodR.Data.example.AAS.gff3())
  txdb <- makeTxDbFromGFF(annotation)
  alias <- data.frame(tx_id = names(id2name(txdb)),
                      name = id2name(txdb))
  expect_error(RNAmodR:::.norm_alias(list(alias = alias)),
               'argument "x" is missing, with no default')
  actual <- RNAmodR:::.norm_alias(list(alias = alias), msi)
  alias <- data.frame(tx_id1 = names(id2name(txdb)),
                      name1 = id2name(txdb))
  expect_error(RNAmodR:::.norm_alias(list(alias = alias), msi),
               "'alias' has to be a data.frame with 'tx_id' and 'name' columns")
  alias <- data.frame(tx_id = names(id2name(txdb)),
                      name1 = id2name(txdb))
  expect_error(RNAmodR:::.norm_alias(list(alias = alias), msi),
               "'alias' has to be a data.frame with 'tx_id' and 'name' columns")
  alias <- data.frame(tx_id = c(1,1),
                      name = c(1,2))
  expect_error(RNAmodR:::.norm_alias(list(alias = alias), msi),
               "Values in 'tx_id' have to be unique")
  alias <- data.frame(tx_id = c(12),
                      name = c(1))
  expect_error(RNAmodR:::.norm_alias(list(alias = alias), msi),
               "All values in 'tx_id' have to be valid transcript ids")
  alias <- data.frame(tx_id = names(id2name(txdb)),
                      name = id2name(txdb))
  alias <- RNAmodR:::.norm_alias(list(alias = alias), msi)
  expect_error(RNAmodR:::.norm_compare_args(),
               'argument "x" is missing, with no default')
  data <- subsetByCoord(msi, coord)
  actual <- RNAmodR:::.norm_compare_args(list(),data,msi)
  expect_type(actual,"list")
  expect_named(actual,c("alias","compareType","perTranscript","sequenceData"))
  expect_error(RNAmodR:::.norm_compare_args(list(perTranscript = 1),data,msi),
               "'perTranscript' must be a single logical value")
  expect_error(RNAmodR:::.norm_compare_args(list(sequenceData = 1),data,msi),
               "'sequenceData' must be a single logical value")
  expect_error(RNAmodR:::.norm_compare_args(list(compareType = 1),data,msi),
               "'compareType' must be a character and a valid colname")
  expect_error(RNAmodR:::.norm_compare_args(list(allTypes = 1),data,msi),
               "'allTypes' must be a single logical value")
  expect_error(RNAmodR:::.norm_compare_args(list(sequenceData = TRUE, 
                                                 allTypes = FALSE),data,
                                            sequenceData(msi[[1]])),
               "'compareType' must be set if 'sequenceData = TRUE' and")
  actual <- compareByCoord(msi, coord)
  expect_s4_class(actual,"DataFrame")
  actual2 <- compareByCoord(msi, split(coord,seq_along(coord)))
  expect_s4_class(actual2,"DataFrame")
  expect_equal(actual,actual2)
  # perTranscript
  actual3 <- compareByCoord(msi, coord, perTranscript = TRUE)
  expect_s4_class(actual3,"DataFrame")
  expect_equal(actual,actual3)
  # plotting
  coord <- modifications(msi[[1]])
  actual <- plotCompareByCoord(msi, coord)
  expect_s3_class(actual,"ggplot")
  actual2 <- plotCompareByCoord(msi, split(coord,seq_along(coord)))
  expect_s3_class(actual2,"ggplot")
})
