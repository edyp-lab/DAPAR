##' This function show the density plots of Fold Change (the same as calculated by limma) for a list 
##' of the comparisons of conditions in a differnetial analysis.
##' 
##' @title Density plots of FC values
##' @param df_FC A dataframe that contains the FC values
##' @param threshold_LogFC The threshold on log(Fold Change) to
##' distinguish between differential and non-differential data 
##' @return A highcharts density plot
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept[1:1000]
##' lapala <- findLapalaBlock(obj)
##' obj <- wrapper.impute.detQuant(obj)
##' obj <- reIntroduceLapala(obj, lapala)
##' obj <- wrapper.impute.detQuant(obj)
##' limma <- wrapper.limmaCompleteTest(obj, 1)
##' hc_FC_DensityPlot(limma$FC)
hc_FC_DensityPlot <-function(df_FC, threshold_LogFC = 0){
    
    if (is.null(df_FC)){return()}

     hc <-  highchart() %>% 
         hc_title(text = "log(FC) repartition") %>% 
         my_hc_chart(chartType = "spline", zoomType="x") %>%
         hc_legend(enabled = TRUE) %>%
         hc_xAxis(title = list(text = "log(FC)"),
                  plotBands = list(list(from= -threshold_LogFC, to = threshold_LogFC, color = "lightgrey")),
                  plotLines=list(list(color= "grey" , width = 2, value = 0, zIndex = 5)))%>%
        hc_yAxis(title = list(text="Density")) %>%
         hc_tooltip(headerFormat= '',
                    pointFormat = "<b> {series.name} </b>: {point.y} ",
                    valueDecimals = 2) %>%
         my_hc_ExportMenu(filename = "densityplot") %>%
         hc_plotOptions(
             series=list(
                 animation=list(
                     duration = 100
                 ),
                 connectNulls= TRUE,
                 marker=list(
                     enabled = FALSE)
             )
         )
     
     myColors <- getPaletteForLabels_HC(ncol(df_FC))
     
    for (i in 1:ncol(df_FC)){
        tmp <- density(df_FC[,i])
        hc <- hc_add_series(hc,
                            data.frame(x = tmp$x,  y = tmp$y), 
                            name=colnames(df_FC)[i], 
                            color=myColors[i])
    }
     
     
 return(hc)
    
}






##' This function is a wrappper to the function adjust.p from the
##' cp4p package. It returns the FDR corresponding to the p-values of the 
##' differential analysis.
##' The FDR is computed with the function \code{p.adjust}\{stats\}..
##' 
##' @title Computes the FDR corresponding to the p-values of the 
##' differential analysis using 
##' @param FC The result (FC values) of the differential analysis processed 
##' by \code{\link{wrapper.limmaCompleteTest}} 
##' @param pval The result (p-values) of the differential analysis processed 
##' by \code{\link{wrapper.limmaCompleteTest}} 
##' @param threshold_PVal The threshold on p-pvalue to
##' distinguish between differential and non-differential data 
##' @param threshold_LogFC The threshold on log(Fold Change) to
##' distinguish between differential and non-differential data 
##' @param pi0Method The parameter pi0.method of the method adjust.p 
##' in the package \code{cp4p}
##' @return The computed FDR value (floating number)
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept[1:1000]
##' lapala <- findLapalaBlock(obj)
##' obj <- wrapper.impute.detQuant(obj)
##' obj <- reIntroduceLapala(obj, lapala)
##' obj <- wrapper.impute.detQuant(obj)
##' limma <- wrapper.limmaCompleteTest(obj, 1)
##' diffAnaComputeFDR(limma$FC,limma$P_Value)
diffAnaComputeFDR <- function(FC, pval,threshold_PVal=0, threshold_LogFC = 0, 
                            pi0Method=1){
    if (is.null(FC) || is.null(pval)){return()}
    
    upItems <- which(abs(FC) >= threshold_LogFC)
    
    selectedItems <- pval[upItems,1]

    padj <- adjust.p(selectedItems,  pi0Method)
    
    items <- which(-log10(padj$adjp[,1]) >= threshold_PVal)
    
    BH.fdr <- max(padj$adjp[items,2])

    return(BH.fdr)
}




##' This method returns a class \code{MSnSet} object with the RAW results (logFC and pValue for all
##' comparisons) of a differential analysis processed by \code{\link{wrapper.limmaCompleteTest}}.
##' 
##' @title Returns a \code{MSnSet} object with the RAW results of
##' the differential analysis 
##' @param obj An object of class \code{MSnSet}.
##' @param data The result of the differential analysis processed 
##' by \code{\link{wrapper.limmaCompleteTest}}
##' @param method The method used for differential analysis. 
##' @return A MSnSet
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' limma <- wrapper.limmaCompleteTest(Exp1_R25_pept[1:1000], 1)
##' obj <- diffAnaSaveRAW_Data(Exp1_R25_pept[1:1000], limma)
diffAnaSaveRAW_Data <- function(obj, data,method="limma"){
    if (is.null(data)){
        warning("The differential analysis has not been completed. Maybe there 
                are some missing values in the dataset. If so, please impute before
                running differential analysis")
        return(NULL)}
    
    .fc <- as.data.frame(data$FC)
    .pval <- as.data.frame(data$P_Value)
    for (i in 1:ncol(.fc)){
        Biobase::fData(obj) <- cbind(Biobase::fData(obj), .fc[,i], .pval[,i])
        coln <- colnames(Biobase::fData(obj))
        colnames(Biobase::fData(obj))[(length(coln)-1):length(coln)] <- c(colnames(data$FC)[i],colnames(data$P_Value)[i])
    }
   
    text <- paste("Differential analysis with",method)
    obj@processingData@processing <- c(obj@processingData@processing, text)
    obj@experimentData@other$method = method
    obj@experimentData@other$RawPValues <- TRUE
    return(obj)
}

##' This method returns a class \code{MSnSet} object with the results
##' of differential analysis.
##' 
##' @title Returns a \code{MSnSet} object with the results of
##' the differential analysis performed with \code{\link{limma}} package. 
##' @param obj An object of class \code{MSnSet}.
##' @param data The result of the differential analysis processed 
##' by \code{\link{wrapper.limmaCompleteTest}} 
##' @param method The method used for differential analysis. 
##' Available choices are : "limma", "Welch"
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
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept
##' lapala <- findLapalaBlock(obj)
##' obj <- wrapper.impute.detQuant(obj)
##' obj <- reIntroduceLapala(obj, lapala)
##' obj <- wrapper.impute.detQuant(obj)
##' limma <- wrapper.limmaCompleteTest(obj, 1)
##' fc <- limma$FC[1]
##' pval <- limma$P_Value[1]
##' diffAnaSave(obj, list(FC=fc, P_Value = pval))
diffAnaSave <- function (obj, data, method="limma", 
                        threshold_pVal=1e-60, threshold_logFC=0, fdr=0, 
                        calibrationMethod = "pounds"){
    if (is.null(data)){
        warning("The differential analysis has not been completed. Maybe there 
            are some missing values in the dataset. If so, please impute before
            running differential analysis")
        return(NULL)}
    
    # temp <- obj
    #####################################################################
    
    Biobase::fData(obj)$P_Value <- data$P_Value
    Biobase::fData(obj)$FC <- data$FC
    Biobase::fData(obj)$Significant <- 0

    text <- paste("Differential analysis with",method)
    obj@processingData@processing <- c(obj@processingData@processing, text)
    
    
    ##setSignificant
    x <- Biobase::fData(obj)$FC
    y <- -log10(Biobase::fData(obj)$P_Value)
    
    ipval <- which(y >= threshold_pVal)
    ilogfc <- which(abs(x) >= threshold_logFC)
    Biobase::fData(obj)[intersect(ipval, ilogfc),]$Significant <- 1
    
    # obj@experimentData@other <- list(obj@experimentData@other,
    #                                 method = method,
    #                                     condition1 = condition1,
    #                                     condition2 = condition2,
    #                                     threshold_p_value = threshold_pVal,
    #                                     threshold.logFC = threshold_logFC,
    #                                     fdr = fdr,
    #                                 calibrationMethod = calibrationMethod)
    # 
    
    obj@experimentData@other$method = method
    obj@experimentData@other$condition1 = data$condition1
    obj@experimentData@other$condition2 = data$condition2
    obj@experimentData@other$threshold_p_value = threshold_pVal
    obj@experimentData@other$threshold_logFC = threshold_logFC
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
##' @param obj An object of class \code{MSnSet}.
##' @return A MSnSet
##' @author Alexia Dorffer
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept
##' lapala <- findLapalaBlock(obj)
##' obj <- wrapper.impute.detQuant(obj)
##' obj <- reIntroduceLapala(obj, lapala)
##' obj <- wrapper.impute.detQuant(obj)
##' limma <- wrapper.limmaCompleteTest(obj, 1)
##' fc <- limma$FC[1]
##' pval <- limma$P_Value[1]
##' obj <- diffAnaSave(obj, list(FC=fc, P_Value = pval))
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



##' This function is a wrapper to the calibration.plot method of the 
##' \code{cp4p} package for use with \code{MSnSet} objects.
##'
##' @title Performs a calibration plot on an \code{MSnSet} object, 
##' calling the \code{cp4p} package functions. 
##' @param vPVal A dataframe that contains quantitative data.
##' @param pi0Method A vector of the conditions (labels) (one label 
##' per sample).
##' @return A plot
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept[1:1000]
##' lapala <- findLapalaBlock(obj)
##' obj <- wrapper.impute.detQuant(obj)
##' obj <- reIntroduceLapala(obj, lapala)
##' obj <- wrapper.impute.detQuant(obj)
##' limma <- wrapper.limmaCompleteTest(obj, 1)
##' wrapperCalibrationPlot(limma$P_Value[,1])
wrapperCalibrationPlot <- function(vPVal, pi0Method="pounds"){

if (is.null(vPVal)){return(NULL)}

p <- calibration.plot(vPVal, pi0.method=pi0Method)

return(p)
}
