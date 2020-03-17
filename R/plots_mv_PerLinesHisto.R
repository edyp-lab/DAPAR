#' This method is a wrapper to plots from a \code{MSnSet} object a 
#' histogram which represents the distribution of the 
#' number of missing values (NA) per lines (ie proteins).
#' 
#' @title Histogram of missing values per lines from an object using highcharter
#' \code{MSnSet}
#' @param obj An object of class \code{MSnSet}.
#' @param indLegend The indice of the column name's in \code{pData()} tab .
#' @param showValues A logical that indicates wether numeric values should be
#' drawn above the bars.
#' @return A histogram
#' @author Alexia Dorffer
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' wrapper.mvPerLinesHisto(Exp1_R25_pept)
#' @export
wrapper.mvPerLinesHisto_HC <- function(obj, indLegend="auto", showValues=FALSE){
  if (is.null(obj)){
    warning("The dataset in NULL. Cannot continue.")
    return(NULL)
  }
  
  qData <- Biobase::exprs(obj)
  samplesData <- Biobase::pData(obj)
  hc <- mvPerLinesHisto_HC(qData, samplesData, indLegend, showValues)
  return(hc)
}




#' This method plots a bar plot which represents the distribution of the 
#' number of missing values (NA) per lines (ie proteins).
#' 
#' @title Bar plot of missing values per lines using highcharter
#' @param qData A dataframe that contains the data to plot.
#' @param samplesData A dataframe which contains informations about 
#' the replicates.
#' @param indLegend The indice of the column name's in \code{pData()} tab 
#' @param showValues A logical that indicates wether numeric values should be
#' drawn above the bars.
#' @return A bar plot
#' @author Florence Combes, Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' qData <- Biobase::exprs(Exp1_R25_pept)
#' samplesData <- Biobase::pData(Exp1_R25_pept)
#' mvPerLinesHisto_HC(qData, samplesData)
#' @export
#' @import highcharter
mvPerLinesHisto_HC <- function(qData, samplesData, indLegend="auto", showValues=FALSE){
  
  if (identical(indLegend,"auto")) { indLegend <- c(2:length(colnames(samplesData)))}
  
  for (j in 1:length(colnames(qData))){
    noms <- NULL
    for (i in 1:length(indLegend)){
      noms <- paste(noms, samplesData[j,indLegend[i]], sep=" ")
    }
    colnames(qData)[j] <- noms
  }
  
  coeffMax <- .1
  
  NbNAPerCol <- colSums(is.na(qData))
  NbNAPerRow <- rowSums(is.na(qData))
  #par(mar = c(10,3, 3, 3))
  
  nb.col <- dim(qData)[2] 
  nb.na <- NbNAPerRow
  temp <- table(NbNAPerRow)
  nb.na2barplot <- c(temp, rep(0,1+ncol(qData)-length(temp)))
  
  if (sum(NbNAPerRow) == 0){
    nb.na2barplot <- rep(0,1+ncol(qData))
  }
  
  df <- data.frame(y=nb.na2barplot[-1])
  myColors = rep("lightgrey",nrow(df))
  myColors[nrow(df)] <- "red"
  
  #df1 <- df2 <- df
  #df2[1:(nrow(df)-1),] <- 0
  #df1 [nrow(df),] <- 0
  
  
  #, series = list( pointWidth = 50)
  
  h1 <-  highchart() %>% 
    hc_title(text = "#[lines] with X NA values") %>% 
    hc_add_series(data = df, type="column", colorByPoint = TRUE) %>%
    hc_colors(myColors) %>%
    hc_plotOptions( column = list(stacking = "normal"),
                    animation=list(duration = 100)) %>%
    hc_legend(enabled = FALSE) %>%
    hc_xAxis(categories = row.names(df), title = list(text = "#[NA values] per line")) %>%
    dapar_hc_ExportMenu(filename = "missingValuesPlot1") %>%
    hc_tooltip(enabled = TRUE,
               headerFormat= '',
               pointFormat = "{point.y} ")
  
  return(h1)
  
}
