context("normalization method")
data(test)
test <- test[-7,]

test_that("Global Alignement, sum by columns", {

norm <- matrix(c(-0.3441,NA,NA,NA,NA,NA,-2.9441,-3.6041,NA,
-2.0495,-1.6995,NA,NA,NA,NA,-4.2995,-5.0095,-1.4395,
-2.3061,-1.9261,-2.1561,NA,NA,NA,-4.4961,-5.1161,-2.0761,
-2.8174,-2.6774,-2.4874,-1.6274,NA,NA,-5.4774,-6.0074,-2.6274,
-2.9344,-3.0144,-2.7444,-2.1644,-2.5044,NA,-5.4944,-6.1944,-2.6344,
-3.0076,-3.1376,-2.8276,-2.1876,-2.5276,-4.1676,-5.4476,-6.0176,-2.8976),9,6)

norm <- round(norm, 4)
colnames(norm) <- c("25fmolR1","25fmolR2","25fmolR3",
                    "50fmolR1","50fmolR2","50fmolR3")
rownames(norm) <- c("1","2","3","4","5","6","8","9","10")

funcNorm <- wrapper.normalizeD2(test, "Global Alignment", "sum by columns")
expect_equal(round(Biobase::exprs(funcNorm),4), norm)

labels <- Biobase::pData(test)[,"Label"]
funcNorm <- normalizeD(Biobase::exprs(test), labels, 
                        "Global Alignment", 
                        "sum by columns")
expect_equal(round(funcNorm,4), norm)
})


test_that("Global Alignment, Alignment on all quantiles", {

norm <- matrix(c(24.795,NA,NA,NA,NA,NA,23.785,21.085,NA,
23.7850,24.1721,NA,NA,NA,NA,22.4979,21.0850,24.7950,
23.5185,24.7950,23.9442,NA,NA,NA,22.0986,21.0850,24.2775,
23.1339,23.7850,24.3478,24.7950,NA,NA,21.8324,21.0850,24.0467,
23.5946,22.7705,23.8987,24.7950,24.3980,NA,21.6422,21.0850,24.1183,
23.7850,23.4519,24.1721,24.7950,24.4356,22.4979,21.4996,21.0850,23.9840),9,6)
colnames(norm) <- c("25fmolR1","25fmolR2","25fmolR3",
                    "50fmolR1","50fmolR2","50fmolR3")
rownames(norm) <- c("1","2","3","4","5","6","8","9","10")

funcNorm <- wrapper.normalizeD2(test, "Global Alignment", "Alignment on all quantiles")
expect_equal(round(Biobase::exprs(funcNorm),4), norm)

labels <- Biobase::pData(test)[,"Label"]
funcNorm <- normalizeD2(Biobase::exprs(test), labels, "Global Alignment", "Alignment on all quantiles")
expect_equal(round(funcNorm,4), norm)
})


# test_that("Global Rescaling, sum by columns", {
# 
# norm <- matrix(c(-0.3441,NA,NA,NA,NA,NA,-2.9441,-3.6041,NA,
# -2.0495,-1.6995,NA,NA,NA,NA,-4.2995,-5.0095,-1.4395,
# -2.3061,-1.9261,-2.1561,NA,NA,NA,-4.4961,-5.1161,-2.0761,
# -2.8174,-2.6774,-2.4874,-1.6274,NA,NA,-5.4774,-6.0074,-2.6274,
# -2.9344,-3.0144,-2.7444,-2.1644,-2.5044,NA,-5.4944,-6.1944,-2.6344,
# -3.0076,-3.1376,-2.8276,-2.1876,-2.5276,-4.1676,-5.4476,-6.0176,-2.8976),9,6)
# norm <- round(norm, 4)
# colnames(norm) <- c("25fmolR1","25fmolR2","25fmolR3",
#                     "50fmolR1","50fmolR2","50fmolR3")
# rownames(norm) <- c("1","2","3","4","5","6","8","9","10")
# 
# funcNorm <- wrapper.normalizeD(test, "Global Rescaling", "sum by columns")
# expect_equal(round(Biobase::exprs(funcNorm),4), norm)
# 
# labels <- Biobase::pData(test)[,"Label"]
# funcNorm <- normalizeD(Biobase::exprs(test), labels, 
#                         "Global Rescaling", 
#                         "sum by columns")
# expect_equal(round(funcNorm,4), norm)
# })


test_that("Quantile Centering, overall", {

norm <- matrix(c(23.9385,NA,NA,NA,NA,NA,21.3385,20.6785,NA,
                 23.4105,23.7605,NA,NA,NA,NA,21.1605,20.4505,24.0205,
                 23.2215,23.6015,23.3715,NA,NA,NA,21.0315,20.4115,23.4515,
                 23.5895,23.7295,23.9195,24.7795,NA,NA,20.9295,20.3995,23.7795,
                 23.3125,23.2325,23.5025,24.0825,23.7425,NA,20.7525,20.0525,23.6125,
                 23.0605,22.9305,23.2405,23.8805,23.5405,21.9005,20.6205,20.0505,23.1705),9,6)
colnames(norm) <- c("25fmolR1","25fmolR2","25fmolR3",
                    "50fmolR1","50fmolR2","50fmolR3")
rownames(norm) <- c("1","2","3","4","5","6","8","9","10")

funcNorm <- wrapper.normalizeD2(test, "Quantile Centering", "overall")
expect_equal(round(Biobase::exprs(funcNorm),4), norm)

labels <- Biobase::pData(test)[,"Label"]
funcNorm <- normalizeD2(Biobase::exprs(test), labels, "Quantile Centering", "overall")
expect_equal(round(funcNorm,4), norm)
})


test_that("Quantile Centering, within conditions", {

norm <- matrix(c(23.7954,NA,NA,NA,NA,NA,21.1954,20.5354,NA,
                 23.2674,23.6174,NA,NA,NA,NA,21.0174,20.3074,23.8774,
                 23.0784,23.4584,23.2284,NA,NA,NA,20.8884,20.2684,23.3084,
                 24.7421,24.8821,25.0721,25.9321,NA,NA,22.0821,21.5521,24.9321,
                 24.4651,24.3851,24.6551,25.2351,24.8951,NA,21.9051,21.2051,24.7651,
                 24.2131,24.0831,24.3931,25.0331,24.6931,23.0531,21.7731,21.2031,24.3231),9,6)
norm <- round(norm, 4)
colnames(norm) <- c("25fmolR1","25fmolR2","25fmolR3",
                    "50fmolR1","50fmolR2","50fmolR3")
rownames(norm) <- c("1","2","3","4","5","6","8","9","10")

funcNorm <- wrapper.normalizeD2(test, "Quantile Centering", "within conditions")
expect_equal(round(Biobase::exprs(funcNorm),4), norm)

labels <- Biobase::pData(test)[,"Label"]
funcNorm <- normalizeD2(Biobase::exprs(test), labels,
                        "Quantile Centering", 
                        "within conditions")
expect_equal(round(funcNorm,4), norm)
})


test_that("Mean Centering, overall", {

norm <- matrix(c(  25.2567,NA,NA,NA,NA,NA,22.6567,21.9967,NA,
                   24.1533,24.5033,NA,NA,NA,NA,21.9033,21.1933,24.7633,
                   24.01,24.39,24.16,NA,NA,NA,21.82,21.20,24.24,
                   23.8748,24.0148,24.2048,25.0648,NA,NA,21.2148,20.6848,24.0648,
                   23.8296,23.7496,24.0196,24.5996,24.2596,NA,21.2696,20.5696,24.1296,
                   23.8756,23.7456,24.0556,24.6956,24.3556,22.7156,21.4356,20.8656,23.9856),9,6)
colnames(norm) <- c("25fmolR1","25fmolR2","25fmolR3",
                    "50fmolR1","50fmolR2","50fmolR3")
rownames(norm) <- c("1","2","3","4","5","6","8","9","10")

funcNorm <- wrapper.normalizeD2(test, "Mean Centering", "overall")
expect_equal(round(Biobase::exprs(funcNorm),4), norm)

labels <- Biobase::pData(test)[,"Label"]
funcNorm <- normalizeD2(Biobase::exprs(test), labels, "Mean Centering", "overall")
expect_equal(round(funcNorm,4), norm)
})


test_that("Mean Centering, within conditions", {

norm <- matrix(c(24.3233,NA,NA,NA,NA,NA,21.7233,21.0633,NA,
23.22,23.57,NA,NA,NA,NA,20.97,20.26,23.83,
23.0767,23.4567,23.2267,NA,NA,NA,20.8867,20.2667,23.3067,
24.8081,24.9481,25.1381,25.9981,NA,NA,22.1481,21.6181,24.9981,
24.7629,24.6829,24.9529,25.5329,25.1929,NA,22.2029,21.5029,25.0629,
24.8089,24.6789,24.9889,25.6289,25.2889,23.6489,22.3689,21.7989,24.9189),9,6)
colnames(norm) <- c("25fmolR1","25fmolR2","25fmolR3",
                    "50fmolR1","50fmolR2","50fmolR3")
rownames(norm) <- c("1","2","3","4","5","6","8","9","10")

funcNorm <- wrapper.normalizeD2(test,  "Mean Centering", "within conditions")
expect_equal(round(Biobase::exprs(funcNorm),4), norm)

labels <- Biobase::pData(test)[,"Label"]
funcNorm <- normalizeD2(Biobase::exprs(test), labels, 
                        "Mean Centering",
                        "within conditions")
expect_equal(round(funcNorm,4), norm)
})


test_that("Mean Centering , overall, scaling", {

norm <- matrix(c(24.4367,NA,NA,NA,NA,NA,22.9281,22.5452,NA,
23.8229,24.0368,NA,NA,NA,NA,22.4476,22.0136,24.1958,
23.8052,24.0750,23.9117,NA,NA,NA,22.2500,21.8098,23.9685,
23.6474,23.7318,23.8462,24.3640,NA,NA,22.0457,21.7265,23.7619,
23.6528,23.5997,23.7790,24.1641,23.9383,NA,21.9529,21.4881,23.8520,
23.7302,23.6332,23.8645,24.3419,24.0883,22.8649,21.9101,21.4849,23.8123),9,6)
colnames(norm) <- c("25fmolR1","25fmolR2","25fmolR3",
                    "50fmolR1","50fmolR2","50fmolR3")
rownames(norm) <- c("1","2","3","4","5","6","8","9","10")

funcNorm <- wrapper.normalizeD2(test, "Mean Centering", "overall", scaling=TRUE)
expect_equal(round(Biobase::exprs(funcNorm),4), norm)

labels <- Biobase::pData(test)[,"Label"]
funcNorm <- normalizeD2(Biobase::exprs(test), labels, 
                        "Mean Centering", 
                        "overall", scaling=TRUE)
expect_equal(round(funcNorm,4), norm)
})


test_that("Mean Centering, within conditions, scaling", {

norm <- matrix(c(23.5033,NA,NA,NA,NA,NA,21.9948,21.6119,NA,
22.8896,23.1035,NA,NA,NA,NA,21.5143,21.0803,23.2624,
22.8718,23.1416,22.9783,NA,NA,NA,21.3167,20.8764,23.0351,
24.5808,24.6651,24.7795,25.2974,NA,NA,22.9790,22.6599,24.6952,
24.5861,24.5330,24.7123,25.0974,24.8717,NA,22.8863,22.4214,24.7853,
24.6636,24.5666,24.7978,25.2752,25.0216,23.7982,22.8434,22.4182,24.7456),9,6)
colnames(norm) <- c("25fmolR1","25fmolR2","25fmolR3",
                    "50fmolR1","50fmolR2","50fmolR3")
rownames(norm) <- c("1","2","3","4","5","6","8","9","10")

funcNorm <- wrapper.normalizeD2(test, 
                                "Mean Centering", 
                                "within conditions", scaling=TRUE)
expect_equal(round(Biobase::exprs(funcNorm),4), norm)

labels <- Biobase::pData(test)[,"Label"]
funcNorm <- normalizeD2(Biobase::exprs(test), labels, 
                        "Mean Centering", 
                        "within conditions", scaling=TRUE)
expect_equal(round(funcNorm,4), norm)

})
