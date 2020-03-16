
#' #' This method plots from a \code{MSnSet} object a histogram of 
#' #' missing values.
#' #' 
#' #' @title Histogram of missing values from a \code{MSnSet} object
#' #' @param obj An object of class \code{MSnSet}.
#' #' @param indLegend The indices of the column name's in \code{pData()} tab.
#' #' @param showValues A logical that indicates wether numeric values should be
#' #' drawn above the bars.
#' #' @return A histogram
#' #' @author Alexia Dorffer
#' #' @examples
#' #' utils::data(Exp1_R25_pept, package='DAPARdata')
#' #' wrapper.mvHisto(Exp1_R25_pept, showValues=TRUE)
#' #' @export
#' wrapper.mvHisto <- function(obj, indLegend="auto", showValues=FALSE){
#'   qData <- Biobase::exprs(obj)
#'   samplesData <- Biobase::pData(obj)
#'   conds <- samplesData[,"Condition"]
#'   mvHisto(qData, samplesData, conds, indLegend, showValues)
#' }


#' This method plots from a \code{MSnSet} object a histogram of 
#' missing values.
#' 
#' @title Histogram of missing values from a \code{MSnSet} object
#' @param obj An object of class \code{MSnSet}.
#' @param indLegend The indices of the column name's in \code{pData()} tab.
#' @param showValues A logical that indicates wether numeric values should be
#' drawn above the bars.
#' @param ... xxx
#' @return A histogram
#' @author Alexia Dorffer
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' wrapper.mvHisto_HC(Exp1_R25_pept, showValues=TRUE)
#' @export
wrapper.mvHisto_HC <- function(obj, indLegend="auto", showValues=FALSE, ...){
  if (is.null(obj)){
    warning("The dataset in NULL. Cannot continue.")
    return(NULL)
  }
  
  qData <- Biobase::exprs(obj)
  samplesData <- Biobase::pData(obj)
  conds <- samplesData[,"Condition"]
  mvHisto_HC(qData, samplesData, conds, indLegend, showValues, ...)
}



#' #' This method plots a histogram of missing values.
#' #' 
#' #' @title Histogram of missing values
#' #' @param qData A dataframe that contains quantitative data.
#' #' @param samplesData A dataframe where lines correspond to samples and 
#' #' columns to the meta-data for those samples.
#' #' @param conds A vector of the conditions (one condition per sample).
#' #' @param indLegend The indices of the column name's in \code{pData()} tab
#' #' @param showValues A logical that indicates wether numeric values should be
#' #' drawn above the bars.
#' #' @param palette xxx
#' #' @return A histogram
#' #' @author Florence Combes, Samuel Wieczorek
#' #' @examples
#' #' utils::data(Exp1_R25_pept, package='DAPARdata')
#' #' qData <- Biobase::exprs(Exp1_R25_pept)
#' #' samplesData <- Biobase::pData(Exp1_R25_pept)
#' #' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
#' #' mvHisto(qData, samplesData, conds, indLegend="auto", showValues=TRUE)
#' #' @export
#' mvHisto <- function(qData, samplesData, conds, indLegend="auto", showValues=FALSE, palette=NULL){
#'   
#'   if (is.null(palette)){
#'     palette <- RColorBrewer::brewer.pal(length(unique(conds)),"Dark2")
#'   }else{
#'     if (length(palette) != ncol(qData)){
#'       warning("The color palette has not the same dimension as the number of samples")
#'       return(NULL)
#'     }
#'   }
#'   
#'   if (identical(indLegend,"auto")) { 
#'     indLegend <- c(2:length(colnames(samplesData)))
#'   }
#'   
#'   
#'   colnames(qData) <- samplesData[,"Condition"]
#'   
#'   coeffMax <- .1
#'   
#'   NbNAPerCol <- colSums(is.na(qData))
#'   NbNAPerRow <- rowSums(is.na(qData))
#'   
#'   if (sum(NbNAPerCol) == 0) {if (sum(NbNAPerCol) == 0){
#'     NbNAPerCol <- rep(0,1+ncol(qData))
#'   }} 
#'   x <- barplot(NbNAPerCol, 
#'                main = "# NA per columns",
#'                col=palette,
#'                las=1,
#'                ylim = c(0, 1.2*max(1,NbNAPerCol)),
#'                names.arg = c(1:18), 
#'                cex.names=1.5,
#'                cex.axis=1.5,
#'                axisnames = FALSE
#'   )
#'   
#'   par(xpd = TRUE)
#'   graphics::text(x, -3,
#'                  label = colnames(qData),
#'                  srt = 45,
#'                  adj=1,
#'                  cex=1.4)
#'   
#' }



#' This method plots a histogram of missing values. Same as the function \code{mvHisto}
#' but uses the package \code{highcharter}
#' 
#' @title Histogram of missing values
#' @param qData A dataframe that contains quantitative data.
#' @param samplesData A dataframe where lines correspond to samples and 
#' columns to the meta-data for those samples.
#' @param conds A vector of the conditions (one condition per sample).
#' @param indLegend The indices of the column name's in \code{pData()} tab
#' @param showValues A logical that indicates wether numeric values should be
#' drawn above the bars.
#' @param palette xxx
#' @return A histogram
#' @author Florence Combes, Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' qData <- Biobase::exprs(Exp1_R25_pept)
#' samplesData <- Biobase::pData(Exp1_R25_pept)
#' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
#' mvHisto_HC(qData, samplesData, conds, indLegend="auto", showValues=TRUE)
#' @export
#' @import highcharter
mvHisto_HC <- function(qData, samplesData, conds, indLegend="auto", 
                       showValues=FALSE, palette = NULL){
  
  palette <- BuildPalette(conds, palette)
  
  if (identical(indLegend,"auto")) { 
    indLegend <- c(2:length(colnames(samplesData)))
  }
  
  
  NbNAPerCol <- colSums(is.na(qData))
  NbNAPerRow <- rowSums(is.na(qData))
  
  df <- data.frame(NbNAPerCol)
  names(df) <- 'y'
  
  
  h1 <-  highchart() %>%
    dapar_hc_chart(chartType = "column") %>%
    hc_title(text = "#NA by replicate") %>%
    hc_add_series(df,type="column", colorByPoint = TRUE) %>%
    hc_colors(palette) %>%
    hc_plotOptions( column = list(stacking = "normal"),
                    animation=list(duration = 100)) %>%
    hc_legend(enabled = FALSE) %>%
    hc_xAxis(categories = conds, title = list(text = "Replicates")) %>%
    dapar_hc_ExportMenu(filename = "missingValuesPlot_3") %>%
    hc_tooltip(headerFormat= '',
               pointFormat = "{point.y}")
  
  
  return(h1)
  
  
  
  
  
}

