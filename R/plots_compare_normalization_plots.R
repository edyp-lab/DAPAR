
#' Wrapper to the function that plot to compare the quantitative proteomics 
#' data before and after normalization
#' 
#' @title Builds a plot from a dataframe
#' @param objBefore A dataframe that contains quantitative data before 
#' normalization.
#' @param objAfter A dataframe that contains quantitative data after 
#' normalization.
#' @param condsForLegend A vector of the conditions (one condition per sample).
#' @param indData2Show A vector of the indices of the columns to show in the 
#' plot. The indices are those of indices of 
#' the columns int the data.frame qDataBefore.
#' @param ... arguments for palette
#' @return A plot
#' @author Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
#' objAfter <- wrapper.normalizeD(Exp1_R25_pept, "QuantileCentering","within conditions")
#' wrapper.compareNormalizationD(Exp1_R25_pept, objAfter, conds)
#' @importFrom Biobase exprs pData
#' @export
wrapper.compareNormalizationD <- function(objBefore, objAfter, 
                                          condsForLegend=NULL,
                                          indData2Show=NULL,
                                          ...){
  
  qDataBefore <- Biobase::exprs(objBefore)
  qDataAfter <- Biobase::exprs(objAfter)
  if (is.null(condsForLegend)){
    condsForLegend <- Biobase::pData(objBefore)[,"Condition"]}
  
  compareNormalizationD(qDataBefore, qDataAfter, condsForLegend, indData2Show, ...)
}

#' Wrapper to the function that plot to compare the quantitative proteomics 
#' data before and after normalization. Same as the function \link{wrapper.compareNormalizationD}
#' but uses the package \code{highcharter}
#' 
#' @title Builds a plot from a dataframe
#' @param objBefore A dataframe that contains quantitative data before 
#' normalization.
#' @param objAfter A dataframe that contains quantitative data after 
#' normalization.
#' @param condsForLegend A vector of the conditions (one condition 
#' per sample).
#' @param indData2Show A vector of the indices of the columns to show in the 
#' plot. The indices are those of indices of 
#' the columns int the data.frame qDataBefore.
#' @param ... arguments for palette
#' @return A plot
#' @author Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
#' objAfter <- wrapper.normalizeD(Exp1_R25_pept, "QuantileCentering", 
#' "within conditions")
#' wrapper.compareNormalizationD_HC(Exp1_R25_pept, objAfter, conds)
#' @importFrom Biobase exprs
#' @export
wrapper.compareNormalizationD_HC <- function(objBefore, objAfter, 
                                             condsForLegend=NULL,
                                             indData2Show=NULL,
                                             ...){
  
  qDataBefore <- Biobase::exprs(objBefore)
  qDataAfter <- Biobase::exprs(objAfter)
  
  DAPAR::compareNormalizationD_HC(qDataBefore, qDataAfter, condsForLegend, indData2Show,...)
}

#' Plot to compare the quantitative proteomics data before and after 
#' normalization
#' 
#' @title Builds a plot from a dataframe
#' @param qDataBefore A dataframe that contains quantitative data before 
#' normalization.
#' @param qDataAfter A dataframe that contains quantitative data after 
#' normalization.
#' @param condsForLegend A vector of the conditions (one condition 
#' per sample).
#' @param indData2Show A vector of the indices of the columns to show in 
#' the plot. The indices are those of indices of 
#' the columns int the data.frame qDataBefore.
#' @param palette xxx
#' @return A plot
#' @author Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' qDataBefore <- Biobase::exprs(Exp1_R25_pept)
#' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
#' objAfter <- wrapper.normalizeD(Exp1_R25_pept,"QuantileCentering","within conditions")
#' compareNormalizationD(qDataBefore, Biobase::exprs(objAfter), conds)
#' @importFrom RColorBrewer brewer.pal
#' @export
compareNormalizationD <- function(qDataBefore,
                                  qDataAfter,
                                  condsForLegend=NULL,
                                  indData2Show=NULL,
                                  palette = NULL){
  
  if (is.null(condsForLegend)) return(NULL)
  if (is.null(indData2Show)) {indData2Show <- c(1:ncol(qDataAfter)) }
  
  if (is.null(palette)){
    tmp <- RColorBrewer::brewer.pal(length(unique(condsForLegend)),"Dark2")[1:length(unique(condsForLegend))]
    
    for (i in 1:ncol(qDataBefore)){
      palette[i] <- tmp[ which(condsForLegend[i] == unique(condsForLegend))]
    }
    
  }else{
    if (length(palette) != ncol(qDataBefore)){
      warning("The color palette has not the same dimension as the number of samples")
      return(NULL)
    }
  }
  
  x <- qDataBefore
  y <- qDataAfter/qDataBefore
  
  lim.x <- range(min(x, na.rm=TRUE), max(x, na.rm=TRUE))
  lim.y <- range(min(y, na.rm=TRUE), max(y, na.rm=TRUE))
  
  
  ##Colors definition
  legendColor <- unique(palette)
  txtLegend <- unique(condsForLegend)
  
  
  
  plot(x=NULL
       ,xlim = lim.x
       ,ylim = lim.y
       , cex = 1
       , axes=TRUE
       , xlab = "Intensities before normalization"
       , ylab = "Intensities after normalization / Intensities before 
    normalization"
       ,cex.lab = 1
       ,cex.axis = 1
       ,cex.main = 3)
  
  
  for (i in indData2Show){
    points(x[,i], y[,i], col = palette[i], cex = 1,pch=16)
  }
  
  legend("topleft"
         , legend = txtLegend
         , col = legendColor
         , pch = 15 
         , bty = "n"
         , pt.cex = 2
         , cex = 1
         , horiz = FALSE
         , inset=c(0,0)
  )
  
  
}





#' Wrapper to the function that plot to compare the quantitative proteomics 
#' data before and after normalization
#' 
#' @title Builds a plot from a dataframe
#' @param objBefore A dataframe that contains quantitative data before 
#' normalization.
#' @param objAfter A dataframe that contains quantitative data after 
#' normalization.
#' @param condsForLegend A vector of the conditions (one condition per sample).
#' @param indData2Show A vector of the indices of the columns to show in the 
#' plot. The indices are those of indices of 
#' the columns int the data.frame qDataBefore.
#' @param ... arguments for palette
#' @return A plot
#' @author Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
#' objAfter <- wrapper.normalizeD(Exp1_R25_pept, "QuantileCentering","within conditions")
#' ids <- Biobase::fData(Exp1_R25_pept)[,Exp1_R25_pept@experimentData@other$proteinId]
#' wrapper.compareNormalizationDSubset(Exp1_R25_pept, objAfter, conds, idsForLegend=ids, subset.view=1:10)
#' @importFrom Biobase exprs pData
#' @export
wrapper.compareNormalizationDSubset <- function(objBefore, objAfter, 
                                                condsForLegend=NULL,
                                                indData2Show=NULL,
                                                ...){
  
  qDataBefore <- Biobase::exprs(objBefore)
  qDataAfter <- Biobase::exprs(objAfter)
  #if (is.null(condsForLegend)){
  #condsForLegend <- Biobase::pData(objBefore)[,"Condition"]}
  
  compareNormalizationDSubset(qDataBefore, qDataAfter, indData2Show, ...)
}





#' Plot to compare the quantitative proteomics data before and after 
#' normalization for a subset of protein
#' 
#' @title Builds a plot from a dataframe
#' @param qDataBefore A dataframe that contains quantitative data before 
#' normalization.
#' @param qDataAfter A dataframe that contains quantitative data after 
#' normalization.
#' @param indData2Show A vector of the indices of the columns to show in 
#' the plot. The indices are those of indices of 
#' the columns int the data.frame qDataBefore.
#' @param idsForLegend A vector of the ids of the row in the data.frame qDataBefore.
#' @param palette xxx
#' @param subset.view A vector of index indicating rows to highlight
#' @return A plot
#' @author Samuel Wieczorek, Anais Courtier
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept
#' qDataBefore <- Biobase::exprs(Exp1_R25_pept)
#' ids <-obj@featureData@data[,obj@experimentData@other$proteinId]
#' objAfter <- wrapper.normalizeD(obj,"QuantileCentering","within conditions")
#' compareNormalizationDSubset(qDataBefore, Biobase::exprs(objAfter), idsForLegend=ids,subset.view=1:10)
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export
compareNormalizationDSubset <- function(qDataBefore,
                                        qDataAfter,
                                        indData2Show=NULL, idsForLegend=NULL,
                                        palette = NULL,subset.view=NULL){
  
  
  if (is.null(idsForLegend)) return(NULL)
  if (is.null(indData2Show)) {indData2Show <- c(1:ncol(qDataAfter)) }
  
  if (is.null(palette)){
    palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(subset.view))
    
  }else{
    if (length(palette) != length(subset.view)){
      warning("The color palette has not the same dimension as the number of proteins")
      return(NULL)
    }
  }
  
  x <- qDataBefore
  y <- qDataAfter/qDataBefore
  
  lim.x <- range(min(x, na.rm=TRUE), max(x, na.rm=TRUE))
  lim.y <- range(min(y, na.rm=TRUE), max(y, na.rm=TRUE))
  
  
  ##Colors definition
  legendColor <- palette
  txtLegend <- idsForLegend[subset.view]
  
  
  plot(x=NULL
       ,xlim = lim.x
       ,ylim = lim.y
       , cex = 1
       , axes=TRUE
       , xlab = "Intensities before normalization"
       , ylab = "Intensities after normalization / Intensities before 
       normalization"
       ,cex.lab = 1
       ,cex.axis = 1
       ,cex.main = 3)
  
  
  for (i in indData2Show){
    points(x[subset.view,i], y[subset.view,i], col = palette, cex = 1,pch=16)
  }
  
  legend("topleft"
         , legend = txtLegend
         , col = legendColor
         , pch = 15 
         , bty = "n"
         , pt.cex = 2
         , cex = 1
         , horiz = FALSE
         , inset=c(0,0)
  )
  
}









#' Plot to compare the quantitative proteomics data before and after 
#' normalization using the library \code{highcharter}
#' 
#' @title Builds a plot from a dataframe. Same as compareNormalizationD but 
#' uses the library \code{highcharter}
#' @param qDataBefore A dataframe that contains quantitative data before 
#' normalization.
#' @param qDataAfter A dataframe that contains quantitative data after 
#' normalization.
#' @param condsForLegend A vector of the conditions (one condition 
#' per sample).
#' @param indData2Show A vector of the indices of the columns to show in 
#' the plot. The indices are those of indices of 
#' the columns int the data.frame qDataBefore.
#' @param palette xxx
#' @return A plot
#' @author Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept[1:1000]
#' qDataBefore <- Biobase::exprs(obj)
#' conds <- Biobase::pData(obj)[,"Condition"]
#' objAfter <- wrapper.normalizeD(obj,"QuantileCentering","within conditions")
#' compareNormalizationD_HC(qDataBefore, Biobase::exprs(objAfter), conds)
#' @importFrom RColorBrewer brewer.pal
#' @import highcharter
#' @export
compareNormalizationD_HC <- function(qDataBefore,
                                     qDataAfter,
                                     condsForLegend=NULL,
                                     indData2Show=NULL,
                                     palette = NULL){
  
  if (is.null(condsForLegend)) return(NULL)
  if (is.null(indData2Show)) {indData2Show <- c(1:ncol(qDataAfter)) }
  
  
  if (is.null(palette)){
    tmp <- RColorBrewer::brewer.pal(length(unique(condsForLegend)),"Dark2")[1:length(unique(condsForLegend))]
    
    for (i in 1:ncol(qDataBefore)){
      palette[i] <- tmp[ which(condsForLegend[i] == unique(condsForLegend))]
    }
    
  }else{
    if (length(palette) != ncol(qDataBefore)){
      warning("The color palette has not the same dimension as the number of samples")
      return(NULL)
    }
  }
  
  x <- qDataBefore
  y <- qDataAfter/qDataBefore
  
  ##Colors definition
  legendColor <- unique(palette)
  txtLegend <- unique(condsForLegend)
  
  
  series <- list()
  for (i in 1:length(indData2Show)){
    tmp <- list(name=condsForLegend[i], data =list_parse(data.frame(x=x[,indData2Show[i]],
                                                                    y=y[,indData2Show[i]])))
    series[[i]] <- tmp
  }
  
  h1 <-  highchart2() %>% 
    dapar_hc_chart( chartType = "scatter") %>%
    hc_add_series_list(series) %>%
    hc_tooltip(enabled= "false" ) %>%
    dapar_hc_ExportMenu(filename = "compareNormalization")
  h1
  return(h1)
  
}
