data(Exp1_R25_pept, package="DAPARdata")
obj <- Exp1_R25_pept

obj <- obj[1:10]

metacellPerLinesHisto_HC(obj, pattern = "Missing")

metacellPerLinesHisto_HC(obj)
metacellPerLinesHisto_HC(obj, pattern = "Quantified")
metacellPerLinesHisto_HC(obj, pattern = "Quant. by direct id")
metacellPerLinesHisto_HC(obj, pattern = "Quant. by recovery")