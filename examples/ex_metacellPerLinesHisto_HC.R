data(Exp1_R25_pept, package="DAPARdata")
obj <- Exp1_R25_pept

obj <- obj[1:10]

metacellPerLinesHisto_HC(obj, pattern = "Missing POV")

metacellPerLinesHisto_HC(obj)
metacellPerLinesHisto_HC(obj, pattern = "Quantified")
metacellPerLinesHisto_HC(obj, pattern = "Quant. by direct id")
metacellPerLinesHisto_HC(obj, pattern = "Quant. by recovery")
metacellPerLinesHisto_HC(obj, pattern = c("Quantified", "Quant. by direct id", "Quant. by recovery"))