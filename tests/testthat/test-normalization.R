context("argument normalization")
test_that("argument normalization:",{
  gff <- system.file("extdata","example1.gff3",package = "RNAmodR.Data")
  fasta <- system.file("extdata","example1.fasta",package = "RNAmodR.Data")
  bam <- system.file("extdata","example_wt_1.bam",package = "RNAmodR.Data")
  # .norm_gff
  expect_error(RNAmodR:::.norm_gff(),'argument "x" is missing')
  expect_error(RNAmodR:::.norm_gff(""),"The gff3 file does not exist")
  # .norm_annotation
  expect_error(RNAmodR:::.norm_annotation(),'argument "annotation" is missing')
  expect_error(RNAmodR:::.norm_annotation(""),'The gff3 file does not exist')
  actual <- RNAmodR:::.norm_annotation(gff)
  expect_s4_class(actual,"TxDb")
  grl <- GenomicFeatures::exonsBy(actual)
  expect_s4_class(grl,"GRangesList")
  expect_equal(length(grl),2L)
  expect_named(grl,c("1","2"))
  expect_equal(actual,RNAmodR:::.norm_annotation(actual))
  expect_equal(grl,RNAmodR:::.norm_annotation(grl))
  # .norm_annotation_GRangesList
  expect_error(RNAmodR:::.norm_annotation_GRangesList(),
               'argument "annotation" is missing')
  expect_error(RNAmodR:::.norm_annotation_GRangesList(""),
               "Elements of 'annotation' GRangesList")
  expect_error(RNAmodR:::.norm_annotation_GRangesList(grl[c(1,1)]),
               "Names of elements in 'annotation' GRangesList must be unique")
  grl2 <- grl
  GenomicRanges::strand(grl2@unlistData) <- c("*","+")
  expect_error(RNAmodR:::.norm_annotation_GRangesList(grl2),
               "Invalid strand information. Strand must either be")
  expect_equal(grl,RNAmodR:::.norm_annotation_GRangesList(grl))
  # .norm_sequences
  expect_error(RNAmodR:::.norm_sequences(),'argument "seq" is missing')
  expect_error(RNAmodR:::.norm_sequences(""),
               'Some or all of the files specified by')
  actual <- RNAmodR:::.norm_sequences(fasta)
  expect_s4_class(actual,"FaFile")
  fafile <- actual
  # .norm_bamfiles
  expect_error(RNAmodR:::.norm_bamfiles(),'argument "x" is missing')
  expect_error(RNAmodR:::.norm_bamfiles(""),'Bam files do not exists at')
  expect_error(RNAmodR:::.norm_bamfiles(bam),
               "Names of BamFileList must either be 'treated' or 'control'")
  actual <- RNAmodR:::.norm_bamfiles(c(treated = bam))
  expect_s4_class(actual,"BamFileList")
  expect_named(actual,c("treated"))
  bf <- Rsamtools::BamFile(c(treated = bam))
  expect_error(RNAmodR:::.norm_bamfiles(bf),
               "Names of BamFileList must either be 'treated' or 'control'")
  actual <- RNAmodR:::.norm_bamfiles(c(treated = bf))
  expect_s4_class(actual,"BamFileList")
  expect_equal(actual, RNAmodR:::.norm_bamfiles(c(Treated = bf)))
  # .bam_header_to_seqinfo
  expect_error(RNAmodR:::.bam_header_to_seqinfo(),'argument "bfl" is missing')
  expect_error(RNAmodR:::.bam_header_to_seqinfo(""),'BamFileList required')
  actual <- RNAmodR:::.bam_header_to_seqinfo(bf)
  expect_s4_class(actual,"Seqinfo")
  seqinfo <- actual
  # .norm_seqinfo
  expect_error(RNAmodR:::.norm_seqinfo(),'argument "seqinfo" is missing')
  expect_error(RNAmodR:::.norm_seqinfo(""),
               'Input is not a Seqinfo object and could not be coerced to one')
  actual <- RNAmodR:::.norm_seqinfo(seqinfo)
  expect_equal(actual,seqinfo)
  # 
  expect_error(RNAmodR:::.norm_seqnames(),'argument "bamfiles" is missing')
  expect_error(RNAmodR:::.norm_seqnames(""),'BamFileList required')
  expect_error(RNAmodR:::.norm_seqnames(bf,""),'Something went wrong.')
  expect_error(RNAmodR:::.norm_seqnames(bf,grl),'argument "sequences" is missing')
  expect_error(RNAmodR:::.norm_seqnames(bf,grl,fasta),'Something went wrong.')
  actual <- RNAmodR:::.norm_seqnames(bf,grl,fafile)
  expect_equal(actual,seqinfo)
  seqinfo2 <- seqinfo[c("chr1","chr2","chr3"),]
  expect_equal(seqinfo2,RNAmodR:::.norm_seqnames(bf,grl,fafile,seqinfo2))
  # .norm_mod
  expect_error(RNAmodR:::.norm_mod(),'argument "mod" is missing')
  expect_error(RNAmodR:::.norm_mod(""),'argument "className" is missing')
  expect_error(RNAmodR:::.norm_mod("",""),"Modification '' as defined for  does")
  expect_error(RNAmodR:::.norm_mod("II",""),"Modification 'II' as defined for  does")
  expect_equal("I",RNAmodR:::.norm_mod("I",""))
  # .norm_modifiertype
  expect_error(RNAmodR:::.norm_modifiertype(),'argument "x" is missing')
  expect_error(RNAmodR:::.norm_modifiertype(""),"Empty string")
  setClass("Mo2dInosine2",contains = "ModInosine")
  expect_error(RNAmodR:::.norm_modifiertype("Mo2dInosine2"),
               "Invalid class name of Modifier class: the string 'Mod' must be present once at the front")
  setClass("InosineMod2",contains = "ModInosine")
  expect_error(RNAmodR:::.norm_modifiertype("InosineMod2"),
               "Invalid class name of Modifier class: the string 'Mod' can only be present once at the front of")
})