context("Agregation peptides to proteins")

#########################################################
test_that("Build Adjacency Matrix", {

data(testWithoutNA)

matShared <- matrix(rep(0,40), 10,4, dimnames=list(1:10, 1:4))
matShared[1:2,1] <- matShared[c(5:8,10),2] <-1 
matShared[c(3:5, 9),3] <- 1
matShared[8:10,4] <- 1
matUnique <- matrix(rep(0,18), 6,3, dimnames=list(c(1:4, 6,7), 1:3))
matUnique[1:2,1] <- matUnique[3:4,3] <- matUnique[5:6,2] <-  1

computedMatUnique  <- BuildAdjacencyMatrix(testWithoutNA, 
                                        "Protein.group.IDs", 
                                        unique=TRUE)
computedMatShared  <- BuildAdjacencyMatrix(testWithoutNA, 
                                        "Protein.group.IDs", 
                                        unique=FALSE)

expect_equal(matUnique, computedMatUnique[,sort(colnames(computedMatUnique))])
expect_equal(matShared, computedMatShared[,sort(colnames(computedMatShared))])
})


#########################################################
test_that("Sum of shared peptides", {

data(testWithoutNA)

protID <- "Protein.group.IDs"

computedMShared <- BuildAdjacencyMatrix(testWithoutNA, protID, unique=FALSE)
sumShared <- matrix(c( 41.13,47.33,47.24,49.20,49.84,49.65,
                        107.67,113.81,113.45,116.73,120.52,118.71,
                        93.59,94.15,93.53,97.31,97.97,98.03,
                        61.90,65.87,65.52,67.93,69.36,69.33),byrow=TRUE, 4,6, 
dimnames=list(1:4, c("25fmolR1", "25fmolR2", "25fmolR3", 
                    "50fmolR1", "50fmolR2", "50fmolR3")))

peptSharedUsed <- matrix(c( 2,2,2,2,2,2,
                            5,5,5,5,5,5,
                            4,4,4,4,4,4,
                            3,3,3,3,3,3),byrow=TRUE, 4,6, 
                    dimnames=list(1:4, c("nb.pep.used.25fmolR1", 
                                        "nb.pep.used.25fmolR2", 
                                        "nb.pep.used.25fmolR3",
                                        "nb.pep.used.50fmolR1", 
                                        "nb.pep.used.50fmolR2", 
                                        "nb.pep.used.50fmolR3")))


sumOfMatShared <- SumPeptides(computedMShared, exprs(testWithoutNA))
expect_equal(sumShared, 
            sumOfMatShared$matfin[sort(rownames(sumOfMatShared$matfin)),])
expect_equal(peptSharedUsed,
            sumOfMatShared$nbpep[sort(rownames(sumOfMatShared$nbpep)),])
})


#########################################################
test_that("Sum of unique peptides", {

data(testWithoutNA)
protID <- "Protein.group.IDs"

computedMUnique <- BuildAdjacencyMatrix(testWithoutNA, protID, unique=TRUE)
sumUnique <- matrix(c(41.13,47.33,47.24,49.20,49.84,49.65,
                        42.00,43.86,43.92,44.75,47.47,45.89,
                        48.94,49.01,48.28,50.58,50.88,50.78),byrow=TRUE, 3,6, 
                    dimnames=list(1:3, c("25fmolR1", "25fmolR2", 
                                        "25fmolR3", "50fmolR1", "50fmolR2", 
                                        "50fmolR3")))

peptUniqueUsed <- matrix(c( 2,2,2,2,2,2,
                            2,2,2,2,2,2,
                            2,2,2,2,2,2),byrow=TRUE, 3,6, 
                            dimnames=list(1:3, c("nb.pep.used.25fmolR1", 
                                                "nb.pep.used.25fmolR2", 
                                                "nb.pep.used.25fmolR3",
                                                "nb.pep.used.50fmolR1", 
                                                "nb.pep.used.50fmolR2", 
                                                "nb.pep.used.50fmolR3")))


sumOfMatUnique <- SumPeptides(computedMUnique, exprs(testWithoutNA))
expect_equal(sumUnique, 
                sumOfMatUnique$matfin[sort(rownames(sumOfMatUnique$matfin)),])
expect_equal(peptUniqueUsed,
                sumOfMatUnique$nbpep[sort(rownames(sumOfMatUnique$nbpep)),])
})



#########################################################
test_that("Mean of unique peptides", {

data(testWithoutNA)
protID <- "Protein.group.IDs"

computedMUnique <- BuildAdjacencyMatrix(testWithoutNA, protID, unique=TRUE)
meanUnique <- matrix(c(20.565,23.665,23.62,24.600,24.920,24.825,
                        21.000,21.930,21.96,22.375,23.735,22.945,
                        24.470,24.505,24.14,25.290,25.440,25.390
                        ),byrow=TRUE, 3,6, 
                    dimnames=list(1:3, c("25fmolR1", "25fmolR2", 
                                        "25fmolR3", 
                                        "50fmolR1", "50fmolR2", "50fmolR3")))

peptUniqueUsed <- matrix(c( 2,2,2,2,2,2,
                            2,2,2,2,2,2,
                            2,2,2,2,2,2),byrow=TRUE, 3,6, 
                            dimnames=list(1:3, c("nb.pep.used.25fmolR1", 
                                                "nb.pep.used.25fmolR2", 
                                                "nb.pep.used.25fmolR3",
                                                "nb.pep.used.50fmolR1", 
                                                "nb.pep.used.50fmolR2", 
                                                "nb.pep.used.50fmolR3")))

meanOfMatUnique <- MeanPeptides(computedMUnique, exprs(testWithoutNA))
expect_equal(meanUnique, 
            meanOfMatUnique$matfin[sort(rownames(meanOfMatUnique$matfin)),])
expect_equal(peptUniqueUsed,
                meanOfMatUnique$nbpep[sort(rownames(meanOfMatUnique$nbpep)),])
})



#########################################################
test_that("Mean of SHARED peptides", {

data(testWithoutNA)
protID <- "Protein.group.IDs"

computedMShared <- BuildAdjacencyMatrix(testWithoutNA, protID, unique=FALSE)
meanShared <- 
    matrix(c(20.56500, 23.66500,  23.6200, 24.60000,  24.9200,  24.8250,
            21.53400, 22.76200,  22.6900, 23.34600,  24.1040,  23.7420,
            23.39750, 23.53750,  23.3825, 24.32750,  24.4925,  24.5075,
            20.63333, 21.95667,  21.8400, 22.64333,  23.1200 , 23.1100),
            byrow=TRUE, 4,6, 
dimnames=list(1:4, c("25fmolR1", "25fmolR2", "25fmolR3", 
                    "50fmolR1", "50fmolR2", "50fmolR3")))

peptSharedUsed <- matrix(c( 2,2,2,2,2,2,
                            5,5,5,5,5,5,
                            4,4,4,4,4,4,
                            3,3,3,3,3,3),byrow=TRUE, 4,6, 
                            dimnames=list(1:4, c("nb.pep.used.25fmolR1", 
                                                "nb.pep.used.25fmolR2", 
                                                "nb.pep.used.25fmolR3",
                                                "nb.pep.used.50fmolR1", 
                                                "nb.pep.used.50fmolR2", 
                                                "nb.pep.used.50fmolR3")))
meanOfMatShared <- MeanPeptides(computedMShared, exprs(testWithoutNA))
expect_equal(meanShared, 
            meanOfMatShared$matfin[sort(rownames(meanOfMatShared$matfin)),]
            , tolerance=1e-5)
expect_equal(peptSharedUsed,
            meanOfMatShared$nbpep[sort(rownames(meanOfMatShared$nbpep)),])
})





#########################################################
test_that("Top 3 of SHARED peptides", {

data(testWithoutNA)
n <- 3
protID <- "Protein.group.IDs"

computedMShared <- BuildAdjacencyMatrix(testWithoutNA, protID, unique=FALSE)
topnShared <- matrix(c( 41.13,47.33,47.24,49.20,49.84,49.65,
                        65.44,71.47,70.93,73.63,74.56,74.10,
                        73.15,73.62,72.91,75.97,76.27,76.15,
                        61.90,65.87,65.52,67.93,69.36,69.33),
                        byrow=TRUE, 4,6, 
                        dimnames=list(1:4, 
                                    c("25fmolR1", "25fmolR2", "25fmolR3"
                                    , "50fmolR1", "50fmolR2", "50fmolR3")))

peptSharedUsed <- matrix(c( 2,2,2,2,2,2,
                            3,3,3,3,3,3,
                            3,3,3,3,3,3,
                            3,3,3,3,3,3),byrow=TRUE, 4,6, 
                            dimnames=list(1:4, c("nb.pep.used.25fmolR1", 
                                                "nb.pep.used.25fmolR2", 
                                                "nb.pep.used.25fmolR3",
                                                "nb.pep.used.50fmolR1", 
                                                "nb.pep.used.50fmolR2", 
                                                "nb.pep.used.50fmolR3")))


topnOfMatShared <- TopnPeptides(computedMShared, exprs(testWithoutNA), n)
expect_equal(topnShared, 
            topnOfMatShared$matfin[sort(rownames(topnOfMatShared$matfin)),]
            , tolerance=1e-5)
expect_equal(peptSharedUsed,
                topnOfMatShared$nbpep[sort(rownames(topnOfMatShared$nbpep)),])
})


#########################################################
test_that("Top 3 of UNIQUE peptides", {

data(testWithoutNA)
n <- 3
protID <- "Protein.group.IDs"

computedMUnique <- BuildAdjacencyMatrix(testWithoutNA, protID, unique=TRUE)
topnUnique <- matrix(c( 41.13,47.33,47.24,49.20,49.84,49.65,
                        42.00,43.86,43.92,44.75,47.47,45.89,
                        48.94,49.01,48.28,50.58,50.88,50.78),
                        byrow=TRUE, 3,6, 
                        dimnames=list(1:3, 
                                    c("25fmolR1", "25fmolR2", "25fmolR3"
                                    , "50fmolR1", "50fmolR2", "50fmolR3")))

peptUniqueUsed <- matrix(rep(2,18),byrow=TRUE, 3,6, 
                            dimnames=list(1:3, c("nb.pep.used.25fmolR1", 
                                                "nb.pep.used.25fmolR2", 
                                                "nb.pep.used.25fmolR3",
                                                "nb.pep.used.50fmolR1", 
                                                "nb.pep.used.50fmolR2", 
                                                "nb.pep.used.50fmolR3")))


topnOfMatUnique <- TopnPeptides(computedMUnique, exprs(testWithoutNA), n)
expect_equal(topnUnique, 
            topnOfMatUnique$matfin[sort(rownames(topnOfMatUnique$matfin)),]
                , tolerance=1e-5)
expect_equal(peptUniqueUsed,
                topnOfMatUnique$nbpep[sort(rownames(topnOfMatUnique$nbpep)),])
})
