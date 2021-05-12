



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
#' obj <- Exp1_R25_pept
#' conds <- Biobase::pData(obj)[,"Condition"]
#' objAfter <- wrapper.normalizeD(obj = obj, method = "QuantileCentering",  
#' conds=conds, type = "within conditions")
#' wrapper.compareNormalizationD_HC(obj, objAfter, conds)
#' 
#' @importFrom Biobase exprs
#' 
#' @export
#' 
#' @importFrom Biobase exprs fData pData
#' 
wrapper.compareNormalizationD_HC <- function(objBefore, 
                                             objAfter, 
                                             condsForLegend=NULL,
                                             ...){
  
  qDataBefore <- Biobase::exprs(objBefore)
  qDataAfter <- Biobase::exprs(objAfter)
  
  compareNormalizationD_HC(qDataBefore, qDataAfter, condsForLegend, ...)
}



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
#' @param keyId xxx
#' 
#' @param conds A vector of the conditions (one condition 
#' per sample).
#' 
#' @param pal xxx
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
#' id <- Biobase::fData(obj)[,obj@experimentData@other$proteinId]
#' pal <- ExtendPalette(2)
#' objAfter <- wrapper.normalizeD(obj, method = "QuantileCentering", 
#' conds =conds, type = "within conditions")
#' compareNormalizationD_HC(qDataBefore=qDataBefore, 
#' qDataAfter=Biobase::exprs(objAfter), keyId = id, conds=conds, n=1000)
#' 
#' compareNormalizationD_HC(qDataBefore=qDataBefore, 
#' qDataAfter=Biobase::exprs(objAfter), keyId = id, pal=pal, subset.view=1:4, 
#' conds=conds, n=100)
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
                                     keyId = NULL,
                                     conds =NULL,
                                     pal = NULL,
                                     subset.view = NULL,
                                     n = 100,
                                     type = 'scatter'){
  
  if (is.null(conds)){
    warning("'conds' is null.")
    return(NULL)
  }
  if (is.null(keyId))
    keyId <- 1:length(qDataBefore)
  
  if (!is.null(subset.view) && length(subset.view) > 0)
  {
    keyId <- keyId[subset.view]
    if (nrow(qDataBefore) > 1)
      if (length(subset.view)==1){
        qDataBefore <- t(qDataBefore[subset.view,])
        qDataAfter <- t(qDataAfter[subset.view,])
      } else {
        qDataBefore <- qDataBefore[subset.view,]
        qDataAfter <- qDataAfter[subset.view,]
      }
  } 
  # else {
  #       qDataBefore <- as_tibble(cbind(t(qDataBefore)))
  #       qDataAfter <- as_tibble(cbind(t(qDataAfter)))
  # }
  
  
  if (!match(type, c('scatter', 'line') )){
    warning("'type' must be equal to 'scatter' or 'line'.")
    return(NULL)
  }
  
  # browser()
  if (is.null(n)){
    n <- seq_len(nrow(qDataBefore))
  } else {
    if (n > nrow(qDataBefore)){
      warning("'n' is higher than the number of rows of datasets. Set to number 
              of rows.")
      n <- nrow(qDataBefore)
    }
    
    ind <- sample(seq_len(nrow(qDataBefore)),n)
    keyId <- keyId[ind]
    if (nrow(qDataBefore) > 1)
      if (length(ind) == 1){
       # qDataBefore <- as_tibble(cbind(t(qDataBefore[ind,])))
       # qDataAfter <- as_tibble(cbind(t(qDataAfter[ind,])))
        qDataBefore <- t(qDataBefore[ind,])
        qDataAfter <- t(qDataAfter[ind,])
      } else {
        #qDataBefore <- as_tibble(cbind(qDataBefore[ind,]))
       # qDataAfter <- as_tibble(cbind(qDataAfter[ind,]))
        qDataBefore <- qDataBefore[ind,]
         qDataAfter <- qDataAfter[ind,]
      }
  }
  
  myColors <- NULL
  if (is.null(pal)){
    warning("Color palette set to default.")
    myColors <-   GetColorsForConditions(conds, ExtendPalette(length(unique(conds))))
  } else {
    if (length(pal) != length(unique(conds))){
      warning("The color palette has not the same dimension as the number of samples")
      myColors <- GetColorsForConditions(conds, ExtendPalette(length(unique(conds))))
    } else 
      myColors <- GetColorsForConditions(conds, pal)
  }
  
  x <- qDataBefore
  y <- qDataAfter/qDataBefore
  
  ##Colors definition
  legendColor <- unique(myColors)
  txtLegend <- unique(conds)
  
  
  series <- list()
  for (i in 1:length(conds)){
    series[[i]] <- list(name=colnames(x)[i],
                        data =list_parse(data.frame(x = x[,i],
                                                    y = y[,i],
                                                    name=keyId)
                        )
    )
  }
  
  h1 <-  highchart() %>% 
    dapar_hc_chart( chartType = type) %>%
    hc_add_series_list(series) %>%
    hc_colors(myColors) %>%
    hc_tooltip(headerFormat= '',pointFormat = "Id: {point.name}") %>%
    dapar_hc_ExportMenu(filename = "compareNormalization")
  h1
  
}