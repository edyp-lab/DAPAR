
#' Densityplot of quantitative proteomics data over samples.
#' 
#' @title Builds a densityplot from a dataframe
#' @param obj xxx
#' @param conds xxx
#' @param legend A vector of the conditions (one condition per sample).
#' @param palette xxx
#' @return A density plot
#' @author Florence Combes, Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
#' densityPlotD(Exp1_R25_pept, conds)
#' @importFrom Biobase exprs pData
#' @importFrom RColorBrewer brewer.pal
#' @import graphics
#' @export
densityPlotD <- function(obj, conds, legend=NULL,palette = NULL){
  
  qData <- Biobase::exprs(obj)
  
  if (is.null(legend) ) { legend <- Biobase::pData(obj)[,"Condition"]}
  
  if (is.null(palette)){
    palette <- RColorBrewer::brewer.pal(length(unique(conds)),"Dark2")
  }else{
    if (length(palette) != ncol(qData)){
      warning("The color palette has not the same dimension as the number of samples")
      return(NULL)
    }
  }
  
  
  ### Range of axis definition
  axis.limits <- matrix(data = 0, nrow = 4, ncol = ncol(qData))
  for (i in 1:ncol(qData)){
    dens <- density(qData[,i], na.rm = TRUE)
    axis.limits[,i] <- c(min(dens$x), max(dens$x), min(dens$y), max(dens$y))
  }
  lim.x <- range(min(axis.limits[1,]), max(axis.limits[2,]))
  lim.y <- range(min(axis.limits[3,]), max(axis.limits[4,]))
  
  
  
  plot(x =NULL
       , ylab ="Density"
       , xlab = "log(intensity)"
       , col = palette
       ,xlim = lim.x
       ,ylim = lim.y
       ,las = 1
       ,cex.lab = 1
       ,cex.axis = 1
       ,cex.main = 3)
  
  for (i in ncol(qData)){
    lines(density(qData[,i], na.rm=TRUE), col = palette[i])
  }
  
  
  legend("topleft"         
         , legend = unique(legend)
         , col = unique(palette)
         , pch = 15 
         , bty = "n"
         , pt.cex = 2
         , cex = 1
         , horiz = FALSE
         , inset=c(0,0)
  )
}




#' Densityplot of quantitative proteomics data over samples. Same as the function \code{\link{densityPlotD}}
#' but uses the package \code{highcharter}
#' 
#' @title Builds a densityplot from a dataframe
#' @param obj xxx
#' @param legend A vector of the conditions (one condition 
#' per sample).
#' @param palette xxx
#' @return A density plot
#' @author Samuel WieczorekistD}}
#' @examples 
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' densityPlotD_HC(Exp1_R25_pept)
#' @importFrom Biobase exprs pData
#' @import highcharter
#' @export
densityPlotD_HC <- function(obj, legend=NULL, palette = NULL){
  
  qData <- Biobase::exprs(obj)
  
  if (is.null(legend) ) { legend<- Biobase::pData(obj)[,"Condition"]}
  
  palette <- BuildPalette(Biobase::pData(obj)[,"Condition"], palette)
  
  
  h1 <-  highcharter::highchart() %>% 
    hc_title(text = "Density plot") %>% 
    dapar_hc_chart(chartType = "spline", zoomType="x") %>%
    hc_legend(enabled = TRUE) %>%
    hc_xAxis(title = list(text = "log(Intensity)")) %>%
    hc_yAxis(title = list(text = "Density")) %>%
    hc_tooltip(headerFormat= '',
               pointFormat = "<b> {series.name} </b>: {point.y} ",
               valueDecimals = 2) %>%
    dapar_hc_ExportMenu(filename = "densityplot") %>%
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
  
  if (!is.null(palette)) {
    if (length(palette) != ncol(qData)){
      warning("The color palette has not the same dimension as the number of samples")
      return(NULL)
    }
    h1 <- h1 %>% hc_colors(palette)
  }
  
  if (is.null(legend)) {
    legend <- paste0("series", 1:ncol(qData))
  }
  
  for (i in 1:ncol(qData)){
    
    tmp <- data.frame(x = density(qData[,i], na.rm = TRUE)$x, 
                      y = density(qData[,i], na.rm = TRUE)$y)
    
    h1 <- h1 %>% hc_add_series(data=list_parse(tmp), name=legend[i]) 
    
  }
  
  
  return(h1)
  
}


