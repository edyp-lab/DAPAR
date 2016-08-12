context("Missing values filter")

test_that("Missing values filter, whole Matrix, no NA lines", {
data(test)

Mout <- matrix(c(23.70,23.49,23.43,24.53,24.96,24.89,
                NA,23.84,23.81,24.67,24.88,24.76,
                NA,NA,23.58,24.86,25.15,25.07,
                NA,NA,NA,25.72,25.73,25.71,
                NA,NA,NA,NA,25.39,25.37,
                NA,NA,NA,NA,NA,23.73,
                21.10,21.24,21.24,21.87,22.40,22.45,
                20.44,20.53,20.62,21.34,21.70,21.88,
                NA,24.10,23.66,24.72,25.26,25.00),
                byrow=TRUE, ncol=6,
                dimnames = list(c(1,2,3,4,5,6,8,9,10), 
                                c("25fmolR1", "25fmolR2", "25fmolR3", 
                                "50fmolR1", "50fmolR2", "50fmolR3"))
)

expect_equal(Biobase::exprs(mvFilter(test, "wholeMatrix", th=as.integer(1))), Mout)
})


test_that("Missing values filter, 
        at least 2 intensities for each condition", {
data(test)


Mout <- matrix(c(23.70,23.49,23.43,24.53,24.96,24.89,
                    NA,23.84,23.81,24.67,24.88,24.76,
                    21.10,21.24,21.24,21.87,22.40,22.45,
                    20.44,20.53,20.62,21.34,21.70,21.88,
                    NA,24.10,23.66,24.72,25.26,25.00),
                byrow=TRUE, ncol=6,
                dimnames = list(c(1,2,8,9,10), c("25fmolR1", "25fmolR2", 
                                                "25fmolR3", "50fmolR1", 
                                                "50fmolR2", "50fmolR3"))
)

expect_equal(Biobase::exprs(mvFilter(test, "allCond", th=as.integer(2))), Mout)
})


test_that("Missing values filter, 
        at least 2 intensities ine one condition", {
data(test)

Mout <- matrix(c(23.70,23.49,23.43,24.53,24.96,24.89,
                    NA,23.84,23.81,24.67,24.88,24.76,
                    NA,NA,23.58,24.86,25.15,25.07,
                    NA,NA,NA,25.72,25.73,25.71,
                    NA,NA,NA,NA,25.39,25.37,
                    21.10,21.24,21.24,21.87,22.40,22.45,
                    20.44,20.53,20.62,21.34,21.70,21.88,
                    NA,24.10,23.66,24.72,25.26,25.00),
                byrow=TRUE, ncol=6,
                dimnames = list(c(1,2,3,4,5,8,9,10), 
                                c("25fmolR1", "25fmolR2", "25fmolR3", 
                                    "50fmolR1", "50fmolR2", "50fmolR3"))
)

expect_equal(Biobase::exprs(mvFilter(test, "atLeastOneCond", th=as.integer(2))), Mout)
})
