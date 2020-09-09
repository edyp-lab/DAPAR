#' 
#' 
#' #' Wrapper to the function that plot to compare the quantitative proteomics 
#' #' data before and after normalization
#' #' 
#' #' @title Builds a plot from a dataframe
#' #' 
#' #' @param objBefore A dataframe that contains quantitative data before 
#' #' normalization.
#' #' 
#' #' @param objAfter A dataframe that contains quantitative data after 
#' #' normalization.
#' #' 
#' #' @param conds A vector of the conditions (one condition per sample).
#' #' 
#' #' @param indData2Show A vector of the indices of the columns to show in the 
#' #' plot. The indices are those of indices of 
#' #' the columns int the data.frame qDataBefore.
#' #' 
#' #' @param ... arguments for palette
#' #' 
#' #' @return A plot
#' #' 
#' #' @author Samuel Wieczorek
#' #' 
#' #' @examples
#' #' utils::data(Exp1_R25_pept, package='DAPARdata')
#' #' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
#' #' objAfter <- wrapper.normalizeD(obj=Exp1_R25_pept, method="QuantileCentering",conds = conds, type="within conditions")
#' #' wrapper.compareNormalizationD(Exp1_R25_pept, objAfter, conds)
#' #' 
#' #' @importFrom Biobase exprs pData
#' #' 
#' #' @export
#' #' 
#' #' @importFrom Biobase exprs fData pData
#' #' 
#' wrapper.compareNormalizationD <- function(objBefore, objAfter, 
#'                                           conds=NULL,
#'                                           indData2Show=NULL,
#'                                           ...){
#'   
#'   qDataBefore <- Biobase::exprs(objBefore)
#'   qDataAfter <- Biobase::exprs(objAfter)
#'   if (is.null(conds)){
#'     conds <- Biobase::pData(objBefore)[,"Condition"]}
#'   
#'   compareNormalizationD(qDataBefore, qDataAfter, conds, indData2Show, ...)
#' }

#' Wrapper to the function that plot to compare the quantitative proteomics 
#' data before and after normalization.
#' 
#' @title Builds a plot from a dataframe
#' 
#' @param objBefore A dataframe that contains quantitative data before 
#' normalization.
#' 
#' @param objAfter A dataframe that contains quantitative data after 
#' normalization.
#' 
#' @param condsForLegend A vector of the conditions (one condition 
#' per sample).
#' 
#' @param ... arguments for palette
#' 
#' @return A plot
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
#' objAfter <- wrapper.normalizeD(obj = Exp1_R25_pept, method = "QuantileCentering",  conds=conds,type = "within conditions")
#' wrapper.compareNormalizationD_HC(Exp1_R25_pept, objAfter, conds)
#' 
#' @importFrom Biobase exprs
#' 
#' @export
#' 
#' @importFrom Biobase exprs fData pData
#' 
wrapper.compareNormalizationD_HC <- function(objBefore, objAfter, 
                                             condsForLegend=NULL,
                                             ...){
  
  qDataBefore <- Biobase::exprs(objBefore)
  qDataAfter <- Biobase::exprs(objAfter)
  
  compareNormalizationD_HC(qDataBefore, qDataAfter, condsForLegend, ...)
}

#' 
#' 
#' #' Plot to compare the quantitative proteomics data before and after 
#' #' normalization
#' #' 
#' #' @title Builds a plot from a dataframe
#' #' 
#' #' @param qDataBefore A dataframe that contains quantitative data before 
#' #' normalization.
#' #' 
#' #' @param qDataAfter A dataframe that contains quantitative data after 
#' #' normalization.
#' #' 
#' #' @param condsForLegend A vector of the conditions (one condition 
#' #' per sample).
#' #' 
#' #' @param indData2Show A vector of the indices of the columns to show in 
#' #' the plot. The indices are those of indices of 
#' #' the columns int the data.frame qDataBefore.
#' #' 
#' #' @param palette xxx
#' #' 
#' #' @return A plot
#' #' 
#' #' @author Samuel Wieczorek
#' #' 
#' #' @examples
#' #' utils::data(Exp1_R25_pept, package='DAPARdata')
#' #' qDataBefore <- Biobase::exprs(Exp1_R25_pept)
#' #' conds <- Biobase::pData(Exp1_R25_pept)$Condition
#' #' objAfter <- wrapper.normalizeD(obj=Exp1_R25_pept, method="QuantileCentering", conds = conds, type="within conditions")
#' #' compareNormalizationD(qDataBefore, Biobase::exprs(objAfter), conds)
#' #'  
#' #' @importFrom RColorBrewer brewer.pal
#' #' 
#' #' @export
#' #' 
#' compareNormalizationD <- function(qDataBefore,
#'                                   qDataAfter,
#'                                   condsForLegend=NULL,
#'                                   indData2Show=NULL,
#'                                   palette = NULL){
#'   
#'   if (is.null(condsForLegend)) return(NULL)
#'   if (is.null(indData2Show)) {indData2Show <- c(1:ncol(qDataAfter)) }
#'   
#'   if (is.null(palette)){
#'     tmp <- RColorBrewer::brewer.pal(length(unique(condsForLegend)),"Dark2")[1:length(unique(condsForLegend))]
#'     
#'     for (i in 1:ncol(qDataBefore)){
#'       palette[i] <- tmp[ which(condsForLegend[i] == unique(condsForLegend))]
#'     }
#'     
#'   }else{
#'     if (length(palette) != ncol(qDataBefore)){
#'       warning("The color palette has not the same dimension as the number of samples")
#'       return(NULL)
#'     }
#'   }
#'   
#'   x <- qDataBefore
#'   y <- qDataAfter/qDataBefore
#'   
#'   lim.x <- range(min(x, na.rm=TRUE), max(x, na.rm=TRUE))
#'   lim.y <- range(min(y, na.rm=TRUE), max(y, na.rm=TRUE))
#'   
#'   
#'   ##Colors definition
#'   legendColor <- unique(palette)
#'   txtLegend <- unique(condsForLegend)
#'   
#'   
#'   
#'   plot(x=NULL
#'        ,xlim = lim.x
#'        ,ylim = lim.y
#'        , cex = 1
#'        , axes=TRUE
#'        , xlab = "Intensities before normalization"
#'        , ylab = "Intensities after normalization / Intensities before 
#'     normalization"
#'        ,cex.lab = 1
#'        ,cex.axis = 1
#'        ,cex.main = 3)
#'   
#'   
#'   for (i in indData2Show){
#'     points(x[,i], y[,i], col = palette[i], cex = 1,pch=16)
#'   }
#'   
#'   legend("topleft"
#'          , legend = txtLegend
#'          , col = legendColor
#'          , pch = 15 
#'          , bty = "n"
#'          , pt.cex = 2
#'          , cex = 1
#'          , horiz = FALSE
#'          , inset=c(0,0)
#'   )
#'   
#'   
#' }



#' Plot to compare the quantitative proteomics data before and after 
#' normalization using the library \code{highcharter}
#' 
#' @title Builds a plot from a dataframe. Same as compareNormalizationD but 
#' uses the library \code{highcharter}
#' 
#' @param qDataBefore A dataframe that contains quantitative data before 
#' normalization.
#' 
#' @param qDataAfter A dataframe that contains quantitative data after 
#' normalization.
#' 
#' @param conds A vector of the conditions (one condition 
#' per sample).
#' 
#' @param palette xxx
#' 
#' @param subset.view xxx
#' 
#' @param n An integer that is equal to the maximum number of displayed points. 
#' This number must be less or equal to the size 
#' of the dataset. If it is less than it, it is a random selection
#' 
#' @param type scatter or line
#' 
#' @return A plot
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_prot, package='DAPARdata')
#' obj <- Exp1_R25_prot
#' qDataBefore <- Biobase::exprs(obj)
#' conds <- Biobase::pData(obj)[,"Condition"]
#' objAfter <- wrapper.normalizeD(obj, method = "QuantileCentering", conds =conds, type = "within conditions")
#' compareNormalizationD_HC(qDataBefore=qDataBefore, qDataAfter=Biobase::exprs(objAfter), conds=conds, n=100)
#' 
#' @import highcharter
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tibble as_tibble
#' @importFrom utils str
#' 
#' @export
#' 
compareNormalizationD_HC <- function(qDataBefore,
                                     qDataAfter,
                                     conds =NULL,
                                     palette = NULL,
                                     subset.view = NULL,
                                     n = 100,
                                     type = 'scatter'){
  
  if (is.null(conds)){
    warning("'conds' is null.")
    return(NULL)
  }
  
  print(str(subset.view))
  print(paste0('subset.view:',subset.view))
  
  if (!is.null(subset.view) && length(subset.view) > 0)
  {
    if (nrow(qDataBefore) > 1)
      if (length(subset.view)==1){
        qDataBefore <- as_tibble(cbind(t(qDataBefore[subset.view,])))
        qDataAfter <- as_tibble(cbind(t(qDataAfter[subset.view,])))
      } else {
        qDataBefore <- as_tibble(cbind(qDataBefore[subset.view,]))
        qDataAfter <- as_tibble(cbind(qDataBefore[subset.view,]))
      }
  }
  
  
  if (!match(type, c('scatter', 'line') )){
    warning("'type' must be equal to 'scatter' or 'line'.")
    return(NULL)
  }
  
  
  if (is.null(n)){
    n <- seq_len(nrow(qDataBefore))
  } else {
    if (n > nrow(qDataBefore)){
      warning("'n' is higher than the number of rows of datasets. Set to number 
              of rows.")
      n <- nrow(qDataBefore)
    }
    
    ind <- sample(seq_len(nrow(qDataBefore)),n)
    if (nrow(qDataBefore) > 1)
      if (length(ind) == 1){
        qDataBefore <- as_tibble(cbind(t(qDataBefore[ind,])))
        qDataAfter <- as_tibble(cbind(t(qDataAfter[ind,])))
      } else {
        qDataBefore <- as_tibble(cbind(qDataBefore[ind,]))
        qDataAfter <- as_tibble(cbind(qDataAfter[ind,]))
      }
  }
  
  palette <- BuildPalette(conds, palette)
  # if (is.null(palette)){palette <- rep("#FFFFFF", ncol(qDataBefore))
  # } else {
  #   if (length(palette) != ncol(qDataBefore)){
  #     warning("The color palette has not the same dimension as the number of samples")
  #     return(NULL)
  #   }
  # }
  
  x <- qDataBefore
  y <- qDataAfter/qDataBefore
  
  ##Colors definition
  legendColor <- unique(palette)
  txtLegend <- unique(conds)
  
  
  series <- list()
  for (i in 1:length(conds)){
    tmp <- list(name=colnames(x)[i],
                data =list_parse(data.frame(x = x[,i],
                                            y = y[,i])
    )
    )
    series[[i]] <- tmp
  }
  
  h1 <-  highchart() %>% 
    dapar_hc_chart( chartType = type) %>%
    hc_add_series_list(series) %>%
    hc_colors(palette) %>%
    hc_tooltip(list(enabled=F)) %>%
    dapar_hc_ExportMenu(filename = "compareNormalization")
  h1
  
}