##' This function is a wrappper to the function adjust.p from the
##' cp4p package. It returns the FDR corresponding to the p-values of the 
##' differential analysis.
##' The FDR is computed with the function \code{p.adjust}\{stats\}..
##' 
##' @title Computes the FDR corresponding to the p-values of the 
##' differential analysis using 
##' @param data The result of the differential analysis processed 
##' by \code{\link{diffAna}} 
##' @param threshold_PVal The threshold on p-pvalue to
##' distinguish between differential and non-differential data 
##' @param threshold_LogFC The threshold on log(Fold Change) to
##' distinguish between differential and non-differential data 
##' @param pi0Method The parameter pi0.method of the method adjust.p 
##' in the package \code{cp4p}
##' @return The computed FDR value (floating number)
##' @author Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' obj <- wrapper.mvImputation(UPSpep25, "QRILC")
##' condition1 <- '25fmol'
##' condition2 <- '10fmol'
##' qData <- Biobase::exprs(obj)
##' samplesData <- Biobase::pData(obj)
##' labels <- Biobase::pData(obj)[,"Label"]
##' limma <- diffAnaLimma(qData,samplesData, labels, condition1, condition2)
##' diffAnaComputeFDR(limma)
diffAnaComputeFDR <- function(data,threshold_PVal=0, threshold_LogFC = 0, 
                            pi0Method=1){
    upItems <- which(abs(data$logFC) >= threshold_LogFC)
    selectedItems <- data[upItems,]$P.Value

    padj <- adjust.p(selectedItems,  pi0Method)
    
    items <- which(-log10(padj$adjp[,1]) >= threshold_PVal)
    
    BH.fdr <- max(padj$adjp[items,2])

    return(BH.fdr)
}


##' This method returns a \code{\link{MSnSet}} object with the results
##' of differential analysis.
##' 
##' @title Returns a \code{\link{MSnSet}} object with the results of
##' the differential analysis performed with \code{\link{limma}} package. 
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param data The result of the differential analysis processed 
##' by \code{\link{diffAna}} 
##' @param method The method used for differential analysis. 
##' Available choices are : "limma", "Welch"
##' @param condition1 A vector containing the names (some values of the slot 
##' "Label" of \code{pData()} of the first condition.
##' @param condition2 A vector containing the names (some values of the slot 
##' "Label" of \code{pData()} of the second condition.
##' @param threshold_pVal A float that indicates the threshold on p-value 
##' choosen to discriminate differential proteins.
##' @param threshold_logFC  A float that indicates the threshold on 
##' log(Fold Change) to discriminatedifferential proteins.
##' @param fdr The FDR based on the values of threshold_pVal and 
##' threshold_logFC
##' @param calibrationMethod The calibration method used to compute the 
##' calibration plot
##' @return A MSnSet
##' @author Alexia Dorffer, Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' condition1 <- '25fmol'
##' condition2 <- '10fmol'
##' limma <- wrapper.diffAnaLimma(UPSpep25, condition1, condition2)
##' obj <- diffAnaSave(UPSpep25, limma, "limma", condition1, condition2)
diffAnaSave <- function (obj, data, method="limma", condition1, condition2, 
                        threshold_pVal=1e-60, threshold_logFC=0, fdr=0, 
                        calibrationMethod = "pounds"){
    if (is.null(data)){
        warning("The differential analysis has not been completed. Maybe there 
            are some missing values in the dataset. If so, please impute before
            running differential analysis")
        return(NULL)}
    
    # temp <- obj
    #####################################################################
    Biobase::fData(obj)$P.Value <- data$P.Value
    Biobase::fData(obj)$logFC <- data$logFC
    Biobase::fData(obj)$Significant <- FALSE

    text <- paste("Differential analysis with",method)
    obj@processingData@processing <- c(obj@processingData@processing, text)
    
    
    ##setSignificant
    x <- Biobase::fData(obj)$logFC
    y <- -log10(Biobase::fData(obj)$P.Value)
    
    ipval <- which(y >= threshold_pVal)
    ilogfc <- which(abs(x) >= threshold_logFC)
    Biobase::fData(obj)[intersect(ipval, ilogfc),]$Significant <- TRUE
    
    # obj@experimentData@other <- list(obj@experimentData@other,
    #                                 method = method,
    #                                     condition1 = condition1,
    #                                     condition2 = condition2,
    #                                     threshold.p.value = threshold_pVal,
    #                                     threshold.logFC = threshold_logFC,
    #                                     fdr = fdr,
    #                                 calibrationMethod = calibrationMethod)
    # 
    
    obj@experimentData@other$method = method
    obj@experimentData@other$condition1 = condition1
    obj@experimentData@other$condition2 = condition2
    obj@experimentData@other$threshold.p.value = threshold_pVal
    obj@experimentData@other$threshold.logFC = threshold_logFC
    obj@experimentData@other$fdr = fdr
    obj@experimentData@other$calibrationMethod = calibrationMethod
    
    
    text <- paste("Differential analysis : Selection with the following 
                    threshold values :logFC =",threshold_logFC,
                    ", -log10(p-value) = ", threshold_pVal,
                    ", FDR = ", fdr, sep=" ")
    
    obj@processingData@processing <- c(obj@processingData@processing, text)
    return(obj)
}


##' Returns a MSnSet object with only proteins significant after 
##' differential analysis.
##' 
##' @title Returns a MSnSet object with only proteins significant after 
##' differential analysis.
##' @param obj An object of class \code{\link{MSnSet}}.
##' @return A MSnSet
##' @author Alexia Dorffer
##' @examples
##' data(UPSpep25)
##' condition1 <- "25fmol"
##' condition2 <- "10fmol"
##' resLimma <- wrapper.diffAnaLimma(UPSpep25, condition1, condition2)
##' obj <-diffAnaSave(UPSpep25, resLimma, "limma", condition1, condition2)
##' signif <- diffAnaGetSignificant(obj)

diffAnaGetSignificant <- function (obj){
    if (is.null(obj)){
        warning("The dataset contains no data")
        return(NULL)
    }
    if (!("Significant" %in% colnames(Biobase::fData(obj)))) {
        warning ("Please Set Significant data before")
        return(NULL)
    }
    temp <- obj
    signif <- which(Biobase::fData(temp)$Significant == TRUE)
    return (temp[signif,])
}



##' Performs a differential analysis on an \code{\link{MSnSet}} object,
##' based on \code{\link{limma}} functions.
##'
##' @title This function performs a differential analysis on an MSnSet
##' object (adapted from \code{\link{limma}})
##' @param qData A dataframe that contains quantitative data.
##' @param design The design matrix as described in the limma package 
##' documentation
##' @return A dataframe with the p-value and log(Fold Change) associated 
##' to each element (peptide/protein)
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' qData <- Biobase::exprs(UPSpep25)
##' design <- cbind(cond1=1, cond2 = rep(0,nrow(Biobase::pData(UPSpep25))))
##' rownames(design) <- rownames(Biobase::pData(UPSpep25))
##' labels <- Biobase::pData(UPSpep25)[,"Label"]
##' indices <- getIndicesConditions(labels, "25fmol", "10fmol")
##' design[indices$iCond2,2] <- 1
##' diffAna(qData, design)
diffAna <- function(qData, design){

    fit <- lmFit(qData, design)
    fit <- eBayes(fit)
    #  fit$df.prior
    diffAna.res <- topTable(fit,
                            coef = 2, 
                            sort.by = "none",
                            number=nrow(qData))
    return (diffAna.res)
}




##' Method to perform differential analysis on
##' a \code{\link{MSnSet}} object (calls the \code{limma} package function).  
##' 
##' @title Performs differential analysis on
##' an MSnSet object, calling the \code{limma} package functions 
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param condition1 A vector that contains the names of the conditions 
##' considered as condition 1.
##' @param condition2 A vector that contains the names of the conditions 
##' considered as condition 2.
##' @return A dataframe as returned by the \code{limma} package
##' @author Alexia Dorffer
##' @examples
##' data(UPSpep25)
##' condition1 <- '25fmol'
##' condition2 <- '10fmol'
##' wrapper.diffAnaLimma(UPSpep25, condition1, condition2)
wrapper.diffAnaLimma <- function(obj, condition1, condition2){

qData <- Biobase::exprs(obj)
samplesData <- Biobase::pData(obj)
labels <- Biobase::pData(obj)[,"Label"]
p <- diffAnaLimma(qData, samplesData, labels, condition1, condition2)
return(p)
}




##' Method to perform differential analysis on
##' an \code{\link{MSnSet}} object (calls the \code{limma} package function).  
##' 
##' @title Performs differential analysis on
##' an MSnSet object, calling the \code{limma} package functions 
##' @param qData A dataframe that contains quantitative data.
##' @param samplesData A dataframe where lines correspond to samples and 
##' columns to the meta-data for those samples.
##' @param labels A vector of the conditions (labels) (one label per sample).
##' @param condition1 A vector that contains the names of the conditions 
##' considered as condition 1
##' @param condition2 A vector that contains the names of the conditions 
##' considered as condition 2
##' @return A dataframe as returned by the \code{limma} package
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' condition1 <- '25fmol'
##' condition2 <- '10fmol'
##' qData <- Biobase::exprs(UPSpep25)
##' samplesData <- Biobase::pData(UPSpep25)
##' labels <- Biobase::pData(UPSpep25)[,"Label"]
##' diffAnaLimma(qData, samplesData, labels, condition1, condition2)
diffAnaLimma <- function(qData, samplesData, labels, condition1, condition2){
if( sum(is.na(qData == TRUE))>0) {
    warning("There are some missing values. Please impute before.")
    return (NULL)
}
# if (condition1 == condition2){
#   warning("The two conditions are identical.")
#   return (NULL)
# }

indices <- getIndicesConditions(labels, condition1, condition2)
flatIndices <- unlist(indices)

tempexprs <- qData[,flatIndices]
design <- cbind(cond1=1, cond2 = rep(0,length(flatIndices)))

rownames(design) <- rownames(samplesData[flatIndices,])
design[flatIndices == indices$iCond2,2] <- 1
res <- diffAna(tempexprs, design)
#print(labels)
#print(c(1:length(labels)))
#res <- limmaCompleteTest(tempexprs, labels, c(1:length(labels)), c(1:length(labels)))
#print(length(res$logFC))
#print(length(rownames(qData)))


p <- data.frame(P.Value = res$P.Value, 
                logFC = res$logFC,
                row.names = rownames(qData))

return(p)
}


##' Computes differential analysis on
##' a \code{\link{MSnSet}} object, using the Welch t-test
##' (\code{\link{t.test}{stats}}). 
##' 
##' @title Performs a differential analysis on a \code{\link{MSnSet}} object
##' using the Welch t-test
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param condition1 A vector containing the names of the conditions 
##' considered as condition 1.
##' @param condition2 A vector containing the names of the conditions 
##' considered as condition 2.
##' @return A dataframe with two slots : P.Value (for the p-value) and logFC
##' (the log of the Fold Change).
##' @author Alexia Dorffer
##' @examples
##' data(UPSpep25)
##' condition1 <- '25fmol'
##' condition2 <- '10fmol'
##' wrapper.diffAnaWelch(UPSpep25, condition1, condition2)
wrapper.diffAnaWelch <- function(obj, condition1, condition2){

qData <- Biobase::exprs(obj)
labels <- Biobase::pData(obj)[,"Label"]
p <- diffAnaWelch(qData, labels, condition1, condition2)
return(p)
}





##' Computes differential analysis on
##' an \code{\link{MSnSet}} object, using the Welch t-test
##' (\code{\link{t.test}{stats}}). 
##' 
##' @title Performs a differential analysis on a \code{\link{MSnSet}} object
##' using the Welch t-test
##' @param qData A dataframe that contains quantitative data.
##' @param labels A vector of the conditions (labels) (one label per sample).
##' @param condition1 A vector containing the names of the conditions 
##' qData as condition 1
##' @param condition2 A vector containing the names of the conditions 
##' considered as condition 2
##' @return A dataframe with two slots : P.Value (for the p-value) and logFC
##' (the log of the Fold Change).
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' condition1 <- '25fmol'
##' condition2 <- '10fmol'
##' qData <- Biobase::exprs(UPSpep25)
##' labels <- Biobase::pData(UPSpep25)[,"Label"]
##' diffAnaWelch(qData, labels, condition1, condition2)
diffAnaWelch <- function(qData, labels, condition1, condition2){

if( sum(is.na(qData == TRUE))>0) {
    warning("There are some missing values. Please impute before.")
    return (NULL)
}
# if (condition1 == condition2){
#   warning("The two conditions are identical.")
#   return (NULL)
# }


indices <- getIndicesConditions(labels, condition1, condition2)
t <- NULL
logRatio <- NULL
for (i in 1:nrow(qData)){
    res <- t.test(x=qData[i,indices$iCond1],
                y=qData[i,indices$iCond2])
    t <- c(t,res$p.value)
    logRatio <- c(logRatio, (res$estimate[2] - res$estimate[1]))
}
p <- data.frame(P.Value=t, logFC = logRatio)
return(p)
}


##' This function is a wrapper to the calibration.plot method of the 
##' \code{cp4p} package for use with \code{\link{MSnSet}} objects.
##'
##' @title Performs a calibration plot on
##' an \code{\link{MSnSet}} object, calling the \code{cp4p} package functions. 
##' @param vPVal A dataframe that contains quantitative data.
##' @param pi0Method A vector of the conditions (labels) (one label 
##' per sample).
##' @return A plot
##' @author Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' condition1 <- '25fmol'
##' condition2 <- '10fmol'
##' qData <- Biobase::exprs(UPSpep25)
##' labels <- Biobase::pData(UPSpep25)[,"Label"]
##' diffAnaWelch(qData, labels, condition1, condition2)
wrapperCalibrationPlot <- function(vPVal, pi0Method="pounds"){

if (is.null(vPVal)){return(NULL)}

p <- calibration.plot(vPVal, pi0.method=pi0Method)

return(p)
}
