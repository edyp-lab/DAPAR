context("Imputation")

require(DAPARdata)

library(testthat)

test_that("KNN", {
  data(Exp1_R25_pept, package="DAPARdata")
  obj <- Exp1_R25_pept[seq_len(10)]
  qdata <- exprs(obj)
  metadata <- GetMetacell(obj)
  
  # Imputation of 'Missing POV'
  obj.imputed.pov <- wrapper.impute.KNN(obj, K = 3)
  qdata.imputed.pov <- exprs(obj.imputed.pov)
  metadata.imputed.pov <- GetMetacell(obj.imputed.pov)
  
  
  ind.missing.pov <- which(metadata=='Missing POV', arr.ind = TRUE)
  res <- metadata.imputed.pov[ind.missing.pov]
  expect_equal(res, rep("Imputed POV", length(res)))
  
  res <- qdata.imputed.pov[ind.missing.pov]
  expect_equal(is.na(res), rep(FALSE, length(res)))
  
  # Check in the qdata if there are no Missing POV now
  
  
})


test_that("slsa", {
  data(Exp1_R25_pept, package="DAPARdata")
  obj <- Exp1_R25_pept[seq_len(10)]
  qdata <- exprs(obj)
  metadata <- GetMetacell(obj)
  
  # Imputation of 'Missing POV'
  obj.imputed <- wrapper.impute.slsa(obj)
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




obj.imp.pov <- wrapper.impute.fixedValue(obj, 0.001, na.type = "Missing POV")
 obj.imp.mec <- wrapper.impute.fixedValue(obj, 0.001, na.type = "Missing MEC")

test_that("fixedValue", {
  data(Exp1_R25_pept, package="DAPARdata")
  obj <- Exp1_R25_pept[seq_len(10)]
  qdata <- exprs(obj)
  metadata <- GetMetacell(obj)

  #
  # Imputation of 'Missing POV'
  #
  obj.imputed.pov <- wrapper.impute.fixedValue(obj, 0.001, na.type = "Missing POV")
  qdata.imputed.pov <- exprs(obj.imputed.pov)
  metadata.imputed.pov <- GetMetacell(obj.imputed.pov)


  ind.missing.pov <- which(metadata=='Missing POV', arr.ind = TRUE)
  res <- metadata.imputed.pov[ind.missing.pov]
  expect_equal(res, rep("Imputed POV", length(res)))

  res <- qdata.imputed.pov[ind.missing.pov]
  expect_equal(is.na(res), rep(FALSE, length(res)))

  #
  # Check imputation of 'Missing MEC'
  #
  obj.imputed.mec <- wrapper.impute.fixedValue(obj, 0.001, na.type = "Missing MEC")
  qdata.imputed.mec <- exprs(obj.imputed.mec)
  metadata.imputed.mec <- GetMetacell(obj.imputed.mec)
  
  
  ind.missing.mec <- which(metadata=='Missing MEC', arr.ind = TRUE)
  res <- metadata.imputed.mec[ind.missing.mec]
  expect_equal(res, rep("Imputed MEC", length(res)))
  
  res <- qdata.imputed.mec[ind.missing.mec]
  expect_equal(is.na(res), rep(FALSE, length(res)))

})