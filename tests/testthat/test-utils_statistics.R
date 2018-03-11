library(RNAmod)

context("General functionality")

test_that("ranges utils functions:",
          {
            # RiboMethScore
            n <- c("100" = 100)
            dataL <- c("94" = 1000,
                       "95" = 1000,
                       "96" = 1000,
                       "97" = 1000,
                       "98" = 1000,
                       "99" = 1000) 
            dataR <- c("101" = 1000,
                       "102" = 1000,
                       "103" = 1000,
                       "104" = 1000,
                       "105" = 1000,
                       "106" = 1000) 
            weights <- c("-6" = 0.5,
                         "-5" = 0.6,
                         "-4" = 0.7,
                         "-3" = 0.8,
                         "-2" = 0.9,
                         "-1" = 1,
                         "0" = 0,
                         "1" = 1,
                         "2" = 0.9,
                         "3" = 0.8,
                         "4" = 0.7,
                         "5" = 0.6,
                         "6" = 0.5)
            actual <- .calculate_ribometh_score_A(n,
                                                  mean(dataL),
                                                  sd(dataL),
                                                  mean(dataR),
                                                  sd(dataR))
            expect_equal(actual, 0.8174387)
            actual <- .calculate_ribometh_score_A(0,
                                                  mean(dataL),
                                                  sd(dataL),
                                                  mean(dataR),
                                                  sd(dataR))
            expect_equal(actual, 0.999001)
            actual <- .calculate_ribometh_score_B(n,
                                                  dataL,
                                                  dataR,
                                                  weights)
            expect_equal(actual, 1.210121)
            actual <- .calculate_ribometh_score_C(n,
                                                  dataL,
                                                  dataR,
                                                  weights)
            expect_equal(actual, 0.55)
          }
)