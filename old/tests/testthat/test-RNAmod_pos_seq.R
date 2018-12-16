library(RNAmod)

context("ranges utils functions")
test_that("ranges utils functions:",
          {
            load(system.file("data", file = "RDN18.RData", package = "RNAmodR"))
            expect_equal(RNAmodR:::.get_strand(RDN18$gr), "+")
            expect_equal(RNAmodR:::.get_strand(c(RDN18$gr,RDN18$gr)), c("+","+"))
            expect_equal(RNAmodR:::.get_unique_strand(c(RDN18$gr,RDN18$gr)), "+")
            expect_true(RNAmodR:::.is_on_correct_strand(RDN18$gr, "+"))
            expect_false(RNAmodR:::.is_on_correct_strand(RDN18$gr, "-"))
            expect_true(RNAmodR:::.is_on_correct_strand2(RDN18$gr, RDN18$gr))
            expect_false(RNAmodR:::.is_minus_strand(RDN18$gr))
            expect_false(RNAmodR:::.is_on_minus_strand(RDN18$gr))
            
            load(system.file("data", file = "gff.RData", package = "RNAmodR"))
            expect_equal(RNAmodR:::.get_gr_ids(RDN18$gr), "RDN18-1")
            expect_length(RNAmodR:::.get_gr_ids(gff),31698)
            expect_equal(RNAmodR:::.get_unique_identifiers(RDN18$gr), "RDN18-1")
            expect_length(RNAmodR:::.get_unique_identifiers(gff),32082)
            expect_equal(length(RNAmodR:::.subset_rnamod_transcript_features(gff)),
                         30422)
            expect_equal(length(RNAmodR:::.subset_rnamod_containing_features(gff)),
                         13385)
            expect_equal(length(RNAmodR:::.subset_rnamod_transcript_features(RDN18$gr)),
                         0)
            expect_equal(length(RNAmodR:::.subset_rnamod_containing_features(RDN18$gr)),
                         1)
            
            expect_equal(length(RNAmodR:::.subset_gff_for_unique_transcript(gff,
                                                                     "YAL030W")),
                          1)
            gr2 <- RNAmodR:::.subset_gff_for_unique_transcript(gff,
                                                              "YAL030W",
                                                              wo.childs = FALSE)
            expect_equal(length(gr2),
                         7)
            expect_equal(length(RNAmodR:::.get_parent_annotations(gr2)),
                         1)
            expect_equal(length(RNAmodR:::.get_parent_annotations(gr2,
                                                                 "RDN18-1")),
                         1)
            expect_equal(length(RNAmodR:::.get_parent_annotations(gr2)),
                         1)
            expect_equal(length(RNAmodR:::.get_parent_annotations(gr2,
                                                                 doRecursiveSearch = TRUE,
                                                                 IDs = c("RDN18-1"))),
                         0)
            expect_equal(length(RNAmodR:::.get_parent_annotations(gr2,
                                                                 doRecursiveSearch = TRUE,
                                                                 forceSingle = TRUE)),
                         1)
            expect_equal(length(RNAmodR:::.get_parent_annotations(gr2,
                                                                 forceSingle = TRUE,
                                                                 doRecursiveSearch = TRUE,
                                                                 IDs = c("YAL030W"))),
                         1)
            
            expect_equal(RNAmodR:::.order_GRanges(gr2[order(BiocGenerics::end(gr2)),]),
                         gr2)
          })


context("Transcript position functions")
test_that("Transcript position functions:",
          {
            load(system.file("data", file = "gff.RData", package = "RNAmodR"))
            load(system.file("data", file = "RDN18.RData", package = "RNAmodR"))
            # .convert_global_to_local_position
            expect_error(
              RNAmodR:::.convert_global_to_local_position(gff,"","")
            )
            expect_error(
              RNAmodR:::.convert_global_to_local_position(gff,RDN18$gr,"")
            )
            
            n <- 10
            actual <- RNAmodR:::.convert_global_to_local_position(gff,
                                                                 RDN18$gr,
                                                                 RDN18$data[1:n])
            expect_length(actual, n)
            expect_equal(actual, rep(1,n))
            actual <- RNAmodR:::.convert_global_to_local_position(gff,
                                                                 RDN18$gr,
                                                                 RDN18$data[(length(RDN18$data)-n):length(RDN18$data)])
            expect_length(actual, n+1)
            expect_equal(actual, c(rep(1785,10),1787))
            
            # .move_postion
            n <- 20
            actual <- RNAmodR:::.move_postion(1:n,c(12,13,14),"+")
            expect_length(actual, n-3)
            expect_equal(max(actual), 17)
            expect_length(actual[actual == 12],1)
            actual <- RNAmodR:::.move_postion(1:n,c(12,13,14),"-")
            expect_length(actual, n-3)
            expect_equal(min(actual), 4)
            expect_length(actual[actual == 14],1)
            
            # no introns
            actual <- RNAmodR:::.get_intron_positions(gff,"RDN18-1")
            expect_equal(actual, list())
            actual <- RNAmodR:::.get_intron_positions(gff,"YBR056W")
            expect_equal(actual, list())
            # one intron
            actual <- RNAmodR:::.get_intron_positions(gff,"YAL030W")
            expect_length(actual, 1)
            
            # .get_transcript_sequence
            actual <- RNAmodR:::.get_transcript_sequence(gff,"RDN18-1",RDN18$seq)
            expected <- as.character(split(RDN18$seq, 1:length(RDN18$seq)))
            names(expected) <- as.character(1:1800)
            expect_length(actual, 1800)
            expect_equal(actual, expected)
            expect_named(actual, as.character(1:1800))
          })
