context("Utils")

require(DAPARdata)
data(Exp1_R2_pept)
obj <- Exp1_R2_pept
test <- Exp1_R2_pept[130:140]
test2 <- Exp1_R2_pept[20:28]




test_that("getProcessingInfo", {
    expect_equal(getProcessingInfo(Exp1_R2_pept), "Log2 tranformed data")
})


test_that("getNumberOfEmptyLines", {
    expect_equal(getNumberOfEmptyLines(exprs(Exp1_R2_pept)), 715)
})


test_that("getIndicesConditions", {
    labels <- Biobase::pData(Exp1_R2_pept)[,"Label"]
    l <- list(iCond1=c(1,2,3), iCond2=c(4,5,6))
    expect_equal(getIndicesConditions(labels, "10fmol", "5fmol"), l)
})


test_that("getPaletteForLabels", {
    labels <- Biobase::pData(Exp1_R2_pept)[,"Label"]
    expect_equal(getPaletteForLabels(labels), c("#1B9E77", "#1B9E77", "#1B9E77", "#D95F02", "#D95F02", "#D95F02"))
})



test_that("getPaletteForReplicates", {
    n <- nrow(Biobase::pData(Exp1_R2_pept))
expect_equal(getPaletteForReplicates(n), c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"))
})

