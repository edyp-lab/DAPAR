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
#' id <- Biobase::fData(obj)[,obj@experimentData@other$proteinId]
#' objAfter <- wrapper.normalizeD(obj, method = "QuantileCentering", conds =conds, type = "within conditions")
#' compareNormalizationD_HC(qDataBefore=qDataBefore, qDataAfter=Biobase::exprs(objAfter), keyId = id, conds=conds, n=100)
#' 
#' compareNormalizationD_HC(qDataBefore=qDataBefore, qDataAfter=Biobase::exprs(objAfter), keyId = id, subset.view=1:4, conds=conds, n=100)
#' 
#' @import highcharter
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tibble as_tibble
#' @importFrom utils str
#' 
#' @export
#' 
# compareNormalizationD_HC <- function(qDataBefore,
#                                      qDataAfter,
#                                      keyId = NULL,
#                                      conds =NULL,
#                                      palette = NULL,
#                                      subset.view = NULL,
#                                      n = 100,
#                                      type = 'scatter'){
#   
#   if (is.null(conds)){
#     warning("'conds' is null.")
#     return(NULL)
#   }
#  if (is.null(keyId))
#    keyId <- 1:length(qDataBefore)
# 
#   if (!is.null(subset.view) && length(subset.view) > 0)
#   {
#     keyId <- keyId[subset.view]
#     if (nrow(qDataBefore) > 1)
#       if (length(subset.view)==1){
#         qDataBefore <- as_tibble(cbind(t(qDataBefore[subset.view,])))
#         qDataAfter <- as_tibble(cbind(t(qDataAfter[subset.view,])))
#       } else {
#         qDataBefore <- as_tibble(cbind(qDataBefore[subset.view,]))
#         qDataAfter <- as_tibble(cbind(qDataAfter[subset.view,]))
#       }
#   } 
#   # else {
#   #       qDataBefore <- as_tibble(cbind(t(qDataBefore)))
#   #       qDataAfter <- as_tibble(cbind(t(qDataAfter)))
#   # }
#   
#   
#   if (!match(type, c('scatter', 'line') )){
#     warning("'type' must be equal to 'scatter' or 'line'.")
#     return(NULL)
#   }
#   
#  # browser()
#   if (is.null(n)){
#     n <- seq_len(nrow(qDataBefore))
#   } else {
#     if (n > nrow(qDataBefore)){
#       warning("'n' is higher than the number of rows of datasets. Set to number 
#               of rows.")
#       n <- nrow(qDataBefore)
#     }
#     
#     ind <- sample(seq_len(nrow(qDataBefore)),n)
#     keyId <- keyId[ind]
#     if (nrow(qDataBefore) > 1)
#       if (length(ind) == 1){
#         qDataBefore <- as_tibble(cbind(t(qDataBefore[ind,])))
#         qDataAfter <- as_tibble(cbind(t(qDataAfter[ind,])))
#       } else {
#         qDataBefore <- as_tibble(cbind(qDataBefore[ind,]))
#         qDataAfter <- as_tibble(cbind(qDataAfter[ind,]))
#       }
#   }
#   
#   palette <- BuildPalette(conds, palette)
#   # if (is.null(palette)){palette <- rep("#FFFFFF", ncol(qDataBefore))
#   # } else {
#   #   if (length(palette) != ncol(qDataBefore)){
#   #     warning("The color palette has not the same dimension as the number of samples")
#   #     return(NULL)
#   #   }
#   # }
#   
#   x <- qDataBefore
#   y <- qDataAfter/qDataBefore
#   
#   ##Colors definition
#   legendColor <- unique(palette)
#   txtLegend <- unique(conds)
#   
#   
#   series <- list()
#   for (i in 1:length(conds)){
#     tmp <- list(name=colnames(x)[i],
#                 data =list_parse(data.frame(x = x[,i],
#                                             y = y[,i],
#                                             name=keyId)
#     )
#     )
#     series[[i]] <- tmp
#   }
#   
#   h1 <-  highchart() %>% 
#     dapar_hc_chart( chartType = type) %>%
#     hc_add_series_list(series) %>%
#     hc_colors(palette) %>%
#     hc_tooltip(headerFormat= '',pointFormat = "Id: {point.name}") %>%
#     dapar_hc_ExportMenu(filename = "compareNormalization")
#   h1
#   
# }

compareNormalizationD_HC <- function(qDataBefore,
                                     qDataAfter,
                                     keyId = NULL,
                                     conds =NULL,
                                     palette = NULL,
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
  
  myColors <- BuildPalette(conds, palette)
  
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