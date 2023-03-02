context("Imputation")

require(DAPARdata)

library(testthat)



test_that("MLE", {
  data(Exp1_R25_pept, package="DAPARdata")
  obj <- Exp1_R25_pept[seq_len(10)]
  qdata <- exprs(obj)
  metadata <- GetMetacell(obj)
  
  # Imputation of 'Missing POV'
  obj.imputed <- wrapper.impute.mle(obj)
  qdata.imputed <- exprs(obj.imputed)
  metadata.imputed <- GetMetacell(obj.imputed)
  
  
  ind.missing.pov <- which(metadata=='Missing POV', arr.ind = TRUE)
  res <- metadata.imputed[ind.missing.pov]
  expect_equal(res, rep("Imputed POV", length(res)))
  
  res <- qdata.imputed[ind.missing.pov]
  expect_equal(is.na(res), rep(FALSE, length(res)))
  
  # Check imputation of 'Missing MEC'
  ind.missing.mec <- which(metadata=='Missing MEC', arr.ind = TRUE)
  res <- qdata.imputed[ind.missing.mec]
  expect_equal(is.na(res), rep(TRUE, length(res)))
  
  # Check in the qdata if there are no Missing POV now
  
  
})



test_that("impute.mi", {
  data(Exp1_R25_pept, package="DAPARdata")
  obj <- Exp1_R25_pept[seq_len(100)]
  qdata <- exprs(obj)
  metadata <- GetMetacell(obj)

  #
  # Imputation of 'Missing POV' only
  #
  obj.imputed.pov <- wrapper.dapar.impute.mi(obj, nb.iter = 1, lapala = FALSE)
  qdata.imputed.pov <- exprs(obj.imputed.pov)
  metadata.imputed.pov <- GetMetacell(obj.imputed.pov)


  ind.missing.pov <- which(metadata=='Missing POV', arr.ind = TRUE)
  res <- metadata.imputed.pov[ind.missing.pov]
  expect_equal(res, rep("Imputed POV", length(res)))

  res <- qdata.imputed.pov[ind.missing.pov]
  expect_equal(is.na(res), rep(FALSE, length(res)))
  
  ind.missing.mec <- which(metadata=='Missing MEC', arr.ind = TRUE)
  res <- metadata.imputed.pov[ind.missing.mec]
  expect_equal(res, rep("Missing MEC", length(res)))
  
  res <- qdata.imputed.pov[ind.missing.mec]
  expect_equal(is.na(res), rep(TRUE, length(res)))

  
  
  
  #
  # Check imputation of both 'Missing POV' and 'Missing MEC'
  #
  obj.imputed.pov.mec <- wrapper.dapar.impute.mi(obj, nb.iter = 1, lapala = TRUE)
  qdata.imputed.pov.mec <- exprs(obj.imputed.pov.mec)
  metadata.imputed.pov.mec <- GetMetacell(obj.imputed.pov.mec)
  
  
  ind.missing.pov <- which(metadata=='Missing POV', arr.ind = TRUE)
  res <- metadata.imputed.pov.mec[ind.missing.pov]
  expect_equal(res, rep("Imputed POV", length(res)))
  
  res <- qdata.imputed.pov.mec[ind.missing.pov]
  expect_equal(is.na(res), rep(FALSE, length(res)))
  
  ind.missing.mec <- which(metadata=='Missing MEC', arr.ind = TRUE)
  res <- metadata.imputed.pov.mec[ind.missing.mec]
  expect_equal(res, rep("Imputed MEC", length(res)))
  
  res <- qdata.imputed.pov.mec[ind.missing.mec]
  expect_equal(is.na(res), rep(FALSE, length(res)))

})



test_that("impute.pa2", {
  data(Exp1_R25_pept, package="DAPARdata")
  obj <- Exp1_R25_pept[seq_len(100)]
  qdata <- exprs(obj)
  metadata <- GetMetacell(obj)
  
  # Imputation of 'Missing POV'
  obj.imputed <- wrapper.impute.pa2(obj)
  qdata.imputed <- exprs(obj.imputed)
  metadata.imputed <- GetMetacell(obj.imputed)
  
  
  ind.missing.pov <- which(metadata=='Missing POV', arr.ind = TRUE)
  res <- metadata.imputed[ind.missing.pov]
  expect_equal(res, rep("Imputed POV", length(res)))
  
  res <- qdata.imputed[ind.missing.pov]
  expect_equal(is.na(res), rep(FALSE, length(res)))
  
  # Check imputation of 'Missing MEC'
  ind.missing.mec <- which(metadata=='Missing MEC', arr.ind = TRUE)
  res <- metadata.imputed[ind.missing.mec]
  expect_equal(res, rep("Imputed MEC", length(res)))
  
  
  res <- qdata.imputed[ind.missing.mec]
  expect_equal(is.na(res), rep(FALSE, length(res)))
  
  
})
