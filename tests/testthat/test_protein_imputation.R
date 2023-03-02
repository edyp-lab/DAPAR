context("Imputation")

require(DAPARdata)

library(testthat)

test_that("KNN", {
  data(Exp1_R25_pept, package="DAPARdata")
  obj <- Exp1_R25_pept[seq_len(10)]
  qdata <- exprs(obj)
  metadata <- GetMetacell(obj)
  
  # Imputation of 'Missing POV' only
  obj.imputed.pov <- wrapper.impute.KNN(obj, K = 3)
  qdata.imputed.pov <- exprs(obj.imputed.pov)
  metadata.imputed.pov <- GetMetacell(obj.imputed.pov)
  
  
  ind.missing.pov <- which(metadata=='Missing POV', arr.ind = TRUE)
  res <- metadata.imputed.pov[ind.missing.pov]
  expect_equal(res, rep("Imputed POV", length(res)))
  res <- qdata.imputed.pov[ind.missing.pov]
  expect_equal(is.na(res), rep(FALSE, length(res)))
  
  #Check if 'Missing MEC' has not been imputed and nor the corresponding metadata
  ind.missing.mec <- which(metadata=='Missing MEC', arr.ind = TRUE)
  res <- metadata.imputed.pov[ind.missing.mec]
  expect_equal(res, rep("Missing MEC", length(res)))
  res <- qdata.imputed.pov[ind.missing.mec]
  expect_equal(is.na(res), rep(TRUE, length(res)))
  
})


test_that("slsa", {
  data(Exp1_R25_pept, package="DAPARdata")
  obj <- Exp1_R25_pept[seq_len(10)]
  qdata <- exprs(obj)
  metadata <- GetMetacell(obj)
  
  obj.imputed <- wrapper.impute.slsa(obj)
  qdata.imputed <- exprs(obj.imputed)
  metadata.imputed <- GetMetacell(obj.imputed)
  
  
  # Check if 'Missing POV' are really imputed and corresponding metadata updated
  ind.missing.pov <- which(metadata=='Missing POV', arr.ind = TRUE)
  res <- metadata.imputed[ind.missing.pov]
  expect_equal(res, rep("Imputed POV", length(res)))
  res <- qdata.imputed[ind.missing.pov]
  expect_equal(is.na(res), rep(FALSE, length(res)))
  
  #Check if 'Missing MEC' has not been imputed and nor the corresponding metadata
  ind.missing.mec <- which(metadata=='Missing MEC', arr.ind = TRUE)
  res <- metadata.imputed[ind.missing.mec]
  expect_equal(res, rep("Missing MEC", length(res)))
  res <- qdata.imputed[ind.missing.mec]
  expect_equal(is.na(res), rep(TRUE, length(res)))
   
  
})





test_that("fixedValue", {
  data(Exp1_R25_pept, package="DAPARdata")
  obj <- Exp1_R25_pept[seq_len(10)]
  qdata <- exprs(obj)
  metadata <- GetMetacell(obj)

  #
  #
  # Imputation of 'Missing POV' only
  #
  #
  obj.imputed.pov <- wrapper.impute.fixedValue(obj, 0.001, na.type = "Missing POV")
  qdata.imputed.pov <- exprs(obj.imputed.pov)
  metadata.imputed.pov <- GetMetacell(obj.imputed.pov)


  # Check if 'Missing POV' are really imputed and corresponding metadata updated
  ind.missing.pov <- which(metadata=='Missing POV', arr.ind = TRUE)
  res <- metadata.imputed.pov[ind.missing.pov]
  expect_equal(res, rep("Imputed POV", length(res)))
  res <- qdata.imputed.pov[ind.missing.pov]
  expect_equal(is.na(res), rep(FALSE, length(res)))
  
  
  # Check if 'Missing MEC' has not been imputed and nor the corresponding metadata
  ind.missing.pov <- which(metadata=='Missing MEC', arr.ind = TRUE)
  res <- metadata.imputed.pov[ind.missing.pov]
  expect_equal(res, rep("Missing MEC", length(res)))
  res <- qdata.imputed.pov[ind.missing.pov]
  expect_equal(is.na(res), rep(TRUE, length(res)))
  
  

  #
  #
  # Check imputation of 'Missing MEC' only
  #
  #
  obj.imputed.mec <- wrapper.impute.fixedValue(obj, 0.001, na.type = "Missing MEC")
  qdata.imputed.mec <- exprs(obj.imputed.mec)
  metadata.imputed.mec <- GetMetacell(obj.imputed.mec)
  
  
  # Check if 'Missing MEC' are really imputed and corresponding metadata updated
  ind.missing.mec <- which(metadata=='Missing MEC', arr.ind = TRUE)
  res <- metadata.imputed.mec[ind.missing.mec]
  expect_equal(res, rep("Imputed MEC", length(res)))
  res <- qdata.imputed.mec[ind.missing.mec]
  expect_equal(is.na(res), rep(FALSE, length(res)))
  
  # Check if 'Missing POV' has not been imputed and nor the corresponding metadata
  ind.missing.pov <- which(metadata=='Missing POV', arr.ind = TRUE)
  res <- metadata.imputed.mec[ind.missing.pov]
  expect_equal(res, rep("Missing POV", length(res)))
  res <- qdata.imputed.mec[ind.missing.pov]
  expect_equal(is.na(res), rep(TRUE, length(res)))
  
  
  
  #
  #
  # Check imputation of both 'Missing POV', and 'Missing MEC'
  #
  #
  obj.imputed.pov.mec <- wrapper.impute.fixedValue(obj, 0.001, na.type = c("Missing POV", "Missing MEC"))
  qdata.imputed.pov.mec <- exprs(obj.imputed.pov.mec)
  metadata.imputed.pov.mec <- GetMetacell(obj.imputed.pov.mec)
  
  
  # Check if 'Missing MEC' are really imputed and corresponding metadata updated
  ind.missing.mec <- which(metadata=='Missing MEC', arr.ind = TRUE)
  res <- metadata.imputed.pov.mec[ind.missing.mec]
  expect_equal(res, rep("Imputed MEC", length(res)))
  res <- qdata.imputed.pov.mec[ind.missing.mec]
  expect_equal(is.na(res), rep(FALSE, length(res)))
  
  # Check if 'Missing POV' are really imputed and corresponding metadata updated
  ind.missing.pov <- which(metadata=='Missing POV', arr.ind = TRUE)
  res <- metadata.imputed.pov.mec[ind.missing.pov]
  expect_equal(res, rep("Imputed POV", length(res)))
  res <- qdata.imputed.pov.mec[ind.missing.pov]
  expect_equal(is.na(res), rep(FALSE, length(res)))

})




test_that("impute.pa", {
  data(Exp1_R25_pept, package="DAPARdata")
  obj <- Exp1_R25_pept[seq_len(10)]
  qdata <- exprs(obj)
  metadata <- GetMetacell(obj)
  
  # Imputation of 'Missing POV'
  obj.imputed <- wrapper.impute.pa(obj)
  qdata.imputed <- exprs(obj.imputed)
  metadata.imputed <- GetMetacell(obj.imputed)
  
  
  # Check if 'Missing POV' are really imputed and corresponding metadata updated
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




test_that("detQuant", {
  data(Exp1_R25_pept, package="DAPARdata")
  obj <- Exp1_R25_pept[seq_len(10)]
  qdata <- exprs(obj)
  metadata <- GetMetacell(obj)
  
  #
  #
  # Imputation of 'Missing POV' only
  #
  #
  obj.imputed.pov <- wrapper.impute.detQuant(obj,na.type = "Missing POV")
  qdata.imputed.pov <- exprs(obj.imputed.pov)
  metadata.imputed.pov <- GetMetacell(obj.imputed.pov)
  
  
  # Check if 'Missing POV' are really imputed and corresponding metadata updated
  ind.missing.pov <- which(metadata=='Missing POV', arr.ind = TRUE)
  res <- metadata.imputed.pov[ind.missing.pov]
  expect_equal(res, rep("Imputed POV", length(res)))
  res <- qdata.imputed.pov[ind.missing.pov]
  expect_equal(is.na(res), rep(FALSE, length(res)))
  
  
  # Check if 'Missing MEC' has not been imputed and nor the corresponding metadata
  ind.missing.pov <- which(metadata=='Missing MEC', arr.ind = TRUE)
  res <- metadata.imputed.pov[ind.missing.pov]
  expect_equal(res, rep("Missing MEC", length(res)))
  res <- qdata.imputed.pov[ind.missing.pov]
  expect_equal(is.na(res), rep(TRUE, length(res)))
  
  
  
  #
  #
  # Check imputation of 'Missing MEC' only
  #
  #
  obj.imputed.mec <- wrapper.impute.detQuant(obj,na.type = "Missing MEC")
  qdata.imputed.mec <- exprs(obj.imputed.mec)
  metadata.imputed.mec <- GetMetacell(obj.imputed.mec)
  
  
  # Check if 'Missing MEC' are really imputed and corresponding metadata updated
  ind.missing.mec <- which(metadata=='Missing MEC', arr.ind = TRUE)
  res <- metadata.imputed.mec[ind.missing.mec]
  expect_equal(res, rep("Imputed MEC", length(res)))
  res <- qdata.imputed.mec[ind.missing.mec]
  expect_equal(is.na(res), rep(FALSE, length(res)))
  
  # Check if 'Missing POV' has not been imputed and nor the corresponding metadata
  ind.missing.pov <- which(metadata=='Missing POV', arr.ind = TRUE)
  res <- metadata.imputed.mec[ind.missing.pov]
  expect_equal(res, rep("Missing POV", length(res)))
  res <- qdata.imputed.mec[ind.missing.pov]
  expect_equal(is.na(res), rep(TRUE, length(res)))
  
  
  
  #
  #
  # Check imputation of both 'Missing POV', and 'Missing MEC'
  #
  #
  obj.imputed.pov.mec <- wrapper.impute.detQuant(obj, na.type = c("Missing POV", "Missing MEC"))
  qdata.imputed.pov.mec <- exprs(obj.imputed.pov.mec)
  metadata.imputed.pov.mec <- GetMetacell(obj.imputed.pov.mec)
  
  
  # Check if 'Missing MEC' are really imputed and corresponding metadata updated
  ind.missing.mec <- which(metadata=='Missing MEC', arr.ind = TRUE)
  res <- metadata.imputed.pov.mec[ind.missing.mec]
  expect_equal(res, rep("Imputed MEC", length(res)))
  res <- qdata.imputed.pov.mec[ind.missing.mec]
  expect_equal(is.na(res), rep(FALSE, length(res)))
  
  # Check if 'Missing POV' are really imputed and corresponding metadata updated
  ind.missing.pov <- which(metadata=='Missing POV', arr.ind = TRUE)
  res <- metadata.imputed.pov.mec[ind.missing.pov]
  expect_equal(res, rep("Imputed POV", length(res)))
  res <- qdata.imputed.pov.mec[ind.missing.pov]
  expect_equal(is.na(res), rep(FALSE, length(res)))
  
})


