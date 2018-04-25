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
##' lapala <- findMECBlock(obj)
##' obj <- wrapper.impute.detQuant(obj)
##' obj <- reIntroduceMEC(obj, lapala)
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
##' lapala <- findMECBlock(obj)
##' obj <- wrapper.impute.detQuant(obj)
##' obj <- reIntroduceMEC(obj, lapala)
##' obj <- wrapper.impute.detQuant(obj)
##' limma <- wrapper.limmaCompleteTest(obj, 1)
##' diffAnaComputeFDR(limma$FC[,1],limma$P_Value[,1])
diffAnaComputeFDR <- function(FC, pval,threshold_PVal=0, threshold_LogFC = 0, 
                            pi0Method=1){
    if (is.null(FC) || is.null(pval)){return()}
    
    upItems <- which(abs(FC) >= threshold_LogFC)
    
    selectedItems <- pval[upItems]

    padj <- adjust.p(selectedItems,  pi0Method)
    
    items <- which(-log10(padj$adjp[,1]) >= threshold_PVal)
    
    BH.fdr <- max(padj$adjp[items,2])

    return(BH.fdr)
}




##' This method returns a class \code{MSnSet} object with the results
##' of differential analysis.
##' 
##' @title Returns a \code{MSnSet} object with the results of
##' the differential analysis performed with \code{\link{limma}} package. 
##' @param obj An object of class \code{MSnSet}.
##' @param allComp A list of two items which is the result of the function wrapper.limmaCompleteTest or xxxx 
##' @param data The result of the differential analysis processed 
##' by \code{\link{wrapper.limmaCompleteTest}} 
##' @param l.params A list of parameters:
##' comp The name of the comparison 
##' th_pVal A float that indicates the threshold on p-value choosen to discriminate differential proteins.
##' fdr The FDR based on the values of threshold_pVal and threshold_logFC
##' calibMethod The calibration method used to compute the calibration plot
##' design xxxxxx
##' th_logFC xxxx
##' @return A MSnSet
##' @author Alexia Dorffer, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept
##' lapala <- findMECBlock(obj)
##' obj <- wrapper.impute.detQuant(obj)
##' obj <- reIntroduceMEC(obj, lapala)
##' obj <- wrapper.impute.detQuant(obj)
##' allComp <- wrapper.limmaCompleteTest(obj, 1)
##' data <- list(FC=allComp$FC[1], P_Value = allComp$P_Value[1])
##' params <- list(design="OnevsOne", method="limma", th_logFC=0)
##' diffAnaSave(obj, allComp, data, params)
diffAnaSave <- function (obj, allComp, data=NULL, l.params){
    if (is.null(allComp)){
        warning("The differential analysis has not been completed. Maybe there 
            are some missing values in the dataset. If so, please impute before
            running differential analysis")
        return(NULL)}
  
  
  #### Check parameters in l.params
  checkRaw <- is.null(l.params$design) &&
            is.null(l.params$method) &&
          is.null(l.params$th_logFC)

  if (is.null(l.params$th_pval)) l.params$th_pval<-0
  if (is.null(l.params$th_logFC)) l.params$th_logFC<-0
  
    
  ####### SAVE ALL THEPAIRWISE COMPARISON RESULTS
  
  .fc <- as.data.frame(allComp$FC)
  .pval <- as.data.frame(allComp$P_Value)
  cnames <- c(colnames(allComp$FC), colnames(allComp$P_Value))
  ind <- which(colnames(fData(obj)) %in% cnames)
  if (length(ind) > 0) {
      Biobase::fData(obj) <- Biobase::fData(obj)[,-ind]
  }
  for (i in 1:ncol(.fc)){
    Biobase::fData(obj) <- cbind(Biobase::fData(obj), .fc[,i], .pval[,i])
    coln <- colnames(Biobase::fData(obj))
    colnames(Biobase::fData(obj))[(length(coln)-1):length(coln)] <- c(colnames(allComp$FC)[i],colnames(allComp$P_Value)[i])
  }
  
  text <- paste("Differential analysis with",l.params$method)
  obj@processingData@processing <- c(obj@processingData@processing, text)
  #Save parameters
  
  obj@experimentData@other$RawPValues <- TRUE
  
  #### SAVE A COMPARISON ANALYSIS IF EXISTS
  if (!(is.null(data$FC) && is.null(data$P_Value))){
  
        Biobase::fData(obj)$P_Value <- data$P_Value
        Biobase::fData(obj)$FC <- data$FC
        Biobase::fData(obj)$Significant <- 0

        ##setSignificant info
        x <- data$FC
        y <- -log10(data$P_Value)
    
        ipval <- which(y >= l.params$th_pval)
        ilogfc <- which(abs(x) >= l.params$th_logFC)
        Biobase::fData(obj)[intersect(ipval, ilogfc),]$Significant <- 1
   

        l.params[["condition1"]] <-  data$condition1
        l.params[["condition2"]] <-  data$condition2
    
        # text <- paste("Differential analysis : Selection with the following 
        #                 threshold values :logFC =",threshold_logFC,
        #                 ", -log10(p-value) = ", threshold_pVal,
        #                 ", FDR = ", fdr, sep=" ")
        # 
        # obj@processingData@processing <- c(obj@processingData@processing, text)
        # 
    }
  
  obj <- saveParameters(obj, "anaDiff", l.params)
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
##' lapala <- findMECBlock(obj)
##' obj <- wrapper.impute.detQuant(obj)
##' obj <- reIntroduceMEC(obj, lapala)
##' obj <- wrapper.impute.detQuant(obj)
##' limma <- wrapper.limmaCompleteTest(obj, 1)
##' fc <- limma$FC[1]
##' pval <- limma$P_Value[1]
##' params <- list(design="OnevsOne", method="limma", th_logFC=0)
##' obj <- diffAnaSave(obj, limma,list(FC=fc, P_Value = pval), params)
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
##' lapala <- findMECBlock(obj)
##' obj <- wrapper.impute.detQuant(obj)
##' obj <- reIntroduceMEC(obj, lapala)
##' obj <- wrapper.impute.detQuant(obj)
##' limma <- wrapper.limmaCompleteTest(obj, 1)
##' wrapperCalibrationPlot(limma$P_Value[,1])
wrapperCalibrationPlot <- function(vPVal, pi0Method="pounds"){

if (is.null(vPVal)){return(NULL)}

p <- calibration.plot(vPVal, pi0.method=pi0Method)

return(p)
}
