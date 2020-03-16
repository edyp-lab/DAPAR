

#' This method is a wrapper to plots from a \code{MSnSet} object a 
#' bar plot which represents the distribution of the 
#' number of missing values (NA) per lines (ie proteins) and per conditions.
#' 
#' @title Bar plot of missing values per lines and per conditions from an 
#' object \code{MSnSet}
#' @param obj An object of class \code{MSnSet}.
#' @param indLegend The indice of the column name's in \code{pData()} tab .
#' @param showValues A logical that indicates wether numeric values should be
#' drawn above the bars.
#' @return A bar plot
#' @author Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' wrapper.mvPerLinesHistoPerCondition(Exp1_R25_pept)
#' @export
wrapper.mvPerLinesHistoPerCondition <- function(obj, indLegend="auto", 
                                                showValues=FALSE){
  qData <- Biobase::exprs(obj)
  samplesData <- Biobase::pData(obj)
  mvPerLinesHistoPerCondition(qData, samplesData, indLegend, showValues)
}


#' This method is a wrapper to plots (using highcharts) from a \code{MSnSet} object a 
#' bar plot which represents the distribution of the 
#' number of missing values (NA) per lines (ie proteins) and per conditions.
#' 
#' @title Bar plot of missing values per lines and per conditions from an 
#' object \code{MSnSet}
#' @param obj An object of class \code{MSnSet}.
#' @param indLegend The indice of the column name's in \code{pData()} tab .
#' @param showValues A logical that indicates wether numeric values should be
#' drawn above the bars.
#' @param ... xxx
#' @return A bar plot
#' @author Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' wrapper.mvPerLinesHistoPerCondition_HC(Exp1_R25_pept)
#' @export
wrapper.mvPerLinesHistoPerCondition_HC <- function(obj, indLegend="auto", 
                                                   showValues=FALSE, ...){
  if (is.null(obj)){
    warning("The dataset in NULL. Cannot continue.")
    return(NULL)
  }
  qData <- Biobase::exprs(obj)
  samplesData <- Biobase::pData(obj)
  mvPerLinesHistoPerCondition_HC(qData, samplesData, indLegend, showValues, ...)
}


#' This method plots a bar plot which represents the distribution of the 
#' number of missing values (NA) per lines (ie proteins) and per conditions.
#' 
#' @title Bar plot of missing values per lines and per condition
#' @param qData A dataframe that contains quantitative data.
#' @param samplesData A dataframe where lines correspond to samples and 
#' columns to the meta-data for those samples.
#' @param indLegend The indice of the column name's in \code{pData()} tab 
#' @param showValues A logical that indicates wether numeric values should be
#' drawn above the bars.
#' @param palette xxx
#' @return A bar plot
#' @author Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' qData <- Biobase::exprs(Exp1_R25_pept)
#' samplesData <- Biobase::pData(Exp1_R25_pept)
#' mvPerLinesHistoPerCondition(qData, samplesData)
#' @export
mvPerLinesHistoPerCondition <- function(qData, samplesData, indLegend="auto", 
                                        showValues=FALSE, palette=NULL){
  
  
  if (is.null(palette)){
    palette <-RColorBrewer::brewer.pal(ncol(qData),"Dark2")[1:ncol(qData)]
  }else{
    if (length(palette) != ncol(qData)){
      warning("The color palette has not the same dimension as the number of samples")
      return(NULL)
    }
  }
  if (identical(indLegend,"auto")) { indLegend <- c(2:length(colnames(samplesData)))}
  
  nbConditions <- length(unique(samplesData[,"Condition"]))
  
  ncolMatrix <- max(unlist(lapply(unique(samplesData[,"Condition"]), function(x){length(which(samplesData[,"Condition"]==x))})))
  m <- matrix(rep(0, nbConditions*(1+ncolMatrix)), 
              ncol = nbConditions, 
              dimnames=list(seq(0:(ncolMatrix)),unique(samplesData[,"Condition"])))
  
  for (i in unique(samplesData[,"Condition"]))
  {
    nSample <- length(which(samplesData[,"Condition"] == i))
    t <- NULL
    if (nSample == 1) {
      t <- table(as.integer(is.na(qData[,which(samplesData[,"Condition"] == i)])))
    } else {t <- table(rowSums(is.na(qData[,which(samplesData[,"Condition"] == i)])))}
    
    m[as.integer(names(t))+1,i] <- t
  }
  
  m <- t(m)
  
  x <- barplot(m, 
               main = "# lines by # of NA",
               xlab = "# NA per lines",
               names.arg = as.character(0:ncolMatrix), 
               col = palette,
               ylim = c(0, 1.2*max(m)), 
               xpd = FALSE,
               las=1,
               cex.names=1.5,
               cex.axis=1.5,
               beside=TRUE
  )
  
}




#' This method plots a bar plot which represents the distribution of the 
#' number of missing values (NA) per lines (ie proteins) and per conditions.
#' Same as the function \link{mvPerLinesHistoPerCondition} but uses the package
#' \code{highcharter}.
#' 
#' @title Bar plot of missing values per lines and per condition
#' @param qData A dataframe that contains quantitative data.
#' @param samplesData A dataframe where lines correspond to samples and 
#' columns to the meta-data for those samples.
#' @param indLegend The indice of the column name's in \code{pData()} tab 
#' @param showValues A logical that indicates wether numeric values should be
#' drawn above the bars.
#' @param palette xxx
#' @return A bar plot
#' @author Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' qData <- Biobase::exprs(Exp1_R25_pept)
#' samplesData <- Biobase::pData(Exp1_R25_pept)
#' mvPerLinesHistoPerCondition_HC(qData, samplesData)
#' @export
#' @import highcharter
mvPerLinesHistoPerCondition_HC <- function(qData, samplesData, indLegend="auto", 
                                           showValues=FALSE, palette=NULL){
  
  conds <- samplesData[,"Condition"]
  palette <- BuildPalette(conds, palette)
  
  
  if (identical(indLegend,"auto")) { indLegend <- c(2:length(colnames(samplesData)))}
  
  nbConditions <- length(unique(samplesData[,"Condition"]))
  
  ncolMatrix <- max(unlist(lapply(unique(samplesData[,"Condition"]), function(x){length(which(samplesData[,"Condition"]==x))})))
  m <- matrix(rep(0, nbConditions*(1+ncolMatrix)), 
              ncol = nbConditions, 
              dimnames=list(seq(0:(ncolMatrix)),unique(samplesData[,"Condition"])))
  
  for (i in unique(samplesData[,"Condition"]))
  {
    nSample <- length(which(samplesData[,"Condition"] == i))
    t <- NULL
    if (nSample == 1) {
      t <- table(as.integer(is.na(qData[,which(samplesData[,"Condition"] == i)])))
    } else {t <- table(rowSums(is.na(qData[,which(samplesData[,"Condition"] == i)])))}
    
    m[as.integer(names(t))+1,i] <- t
  }
  m <- as.data.frame(m)
  
  rownames(m) <- 0:(nrow(m)-1)
  
  h1 <-  highchart() %>% 
    hc_title(text = "#[lines] with X NA values (condition-wise)") %>% 
    dapar_hc_chart(chartType = "column") %>%
    hc_plotOptions( column = list(stacking = ""),
                    dataLabels = list(enabled = FALSE),
                    animation=list(duration = 100)) %>%
    hc_colors(unique(palette)) %>%
    hc_legend(enabled = FALSE) %>%
    hc_xAxis(categories = row.names(m), title = list(text = "#[NA values] per line (condition-wise)")) %>%
    dapar_hc_ExportMenu(filename = "missingValuesPlot_2") %>%
    hc_tooltip(headerFormat= '',
               pointFormat = "{point.y} ")
  
  for (i in 1:nbConditions){
    h1 <- h1 %>% hc_add_series(data=m[,unique(samplesData[,"Condition"])[i]]) }
  
  
  return(h1)
  
  
  
}

