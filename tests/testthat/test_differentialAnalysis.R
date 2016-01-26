context("Differential analysis FDR")

test_that("compare FDR", {
  datalimma <- matrix(c(1.144844e-07, 1.055972e-07, 8.192425e-05, 2.450354e-06, 7.641383e-06, 5.534496e-08,
                        -1.559135, -1.554906, -1.683824, -1.841546, -1.786713, -1.468805), 6,2)
  colnames(datalimma) <- c("P.Value","logFC")
  rownames(datalimma) <- c("2","10","11","12","13","14")
  datalimma <- as.data.frame(datalimma)
  
  funcFdr <- diffAnaComputeFDR(datalimma)
  FDR <- 8.19242509813458e-05
  
  expect_equal(funcFdr, FDR)
})



test_that("Compute limma", {
  data(testWithoutNA)
  limmaRes <- data.frame(P.Value=c(0.0006,0.1669,0.0038,0.0002,0.0018,0.0334,0.1197,0.0027,0.0013,0.0844),
                         logFC=c(1.2533,3.0767,1.0633,0.9400,0.9000, 1.6300,1.1467,1.0467,1.1100,2.2867),
                         row.names = as.character(seq(1:10)))
  
  expect_equal(round(wrapper.diffAnaLimma(testWithoutNA, "25fmol", "50fmol"),4), limmaRes)
})


test_that("Compute Welch", {
  data(testWithoutNA)
  welchRes <- data.frame(P.Value=c(0.0028,0.2857,0.0194,0.0021,0.0221,0.1120,0.2319,0.0244,0.0131,0.1904),
                         logFC=c(1.2533,3.0767,1.0633,0.9400,0.9000, 1.6300,1.1467,1.0467,1.1100,2.2867),
                         row.names =seq(1:10))
  
  expect_equal(round(wrapper.diffAnaWelch(testWithoutNA, "25fmol", "50fmol"),4), welchRes)
})
