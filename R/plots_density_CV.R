
#' Builds a densityplot of the CV of entities in the exprs() table
#' of an object \code{MSnSet}. The variance is calculated for each 
#' condition present
#' in the dataset (see the slot \code{'Condition'} in the \code{pData()} table).
#' 
#' @title Distribution of CV of entities
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param ... arguments for palette
#' 
#' @return A density plot
#' 
#' @author Alexia Dorffer
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' wrapper.CVDistD(Exp1_R25_pept)
#' 
#' @importFrom Biobase exprs pData
#' 
#' @export
#' 
wrapper.CVDistD <- function(obj, ...){
  qData <- Biobase::exprs(obj)
  conds <- Biobase::pData(obj)[,"Condition"]
  CVDistD(qData, conds, ...)
}


#' Builds a densityplot of the CV of entities in the exprs() table. 
#' of an object \code{MSnSet}. The variance is calculated for each 
#' condition present
#' in the dataset (see the slot \code{'Condition'} in the \code{pData()} table).
#' Same as the function \code{\link{wrapper.CVDistD}} but uses the package \code{highcharter}
#' 
#' @title Distribution of CV of entities
#' 
#' @param obj An object of class \code{MSnSet}
#' 
#' @param ... arguments for palette.
#' 
#' @return A density plot
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' wrapper.CVDistD_HC(Exp1_R25_pept)
#' 
#' @importFrom Biobase exprs pData
#' 
#' 
#' @export
#' 
wrapper.CVDistD_HC <- function(obj, ...){
  qData <- Biobase::exprs(obj)
  conds <- Biobase::pData(obj)[,"Condition"]
  CVDistD_HC(qData, conds, ...)
}


#' Builds a densityplot of the CV of entities in the exprs() table
#' of a object. The CV is calculated for each condition present
#' in the dataset (see the slot \code{'Condition'} in the \code{pData()} table)
#' 
#' @title Distribution of CV of entities
#' 
#' @param qData A dataframe that contains quantitative data.
#' 
#' @param conds A vector of the conditions (one condition per sample).
#' 
#' @param palette xxx
#' 
#' @return A density plot
#' 
#' @author Florence Combes, Samuel Wieczorek
#' 
#' @seealso \code{\link{densityPlotD}}.
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
#' CVDistD(Biobase::exprs(Exp1_R25_pept), conds)
#' 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats density var
#' 
#' @export
#' 
CVDistD <- function(qData, conds=NULL, palette = NULL){
  
  if (is.null(conds)) {return(NULL)}
  if (is.null(palette)){
    #palette <- RColorBrewer::brewer.pal(length(unique(conds)),"Dark2")[1:length(unique(conds))]
    palette <- grDevices::colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(conds)))
  }else{
    if (length(palette) != ncol(qData)){
      warning("The color palette has not the same dimension as the number of samples. Set to default palette.")
      palette <- grDevices::colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(conds)))
      return(NULL)
    }
  }
  
  conditions <- unique(conds)
  n <- length(conditions)
  axis.limits <- matrix(data = 0, nrow = 4, ncol = n)
  for (i in conditions){
    if (length(which(conds == i)) > 1){
      t <- density(apply(qData[,which(conds == i)], 1, 
                         function(x) 100*var(x, na.rm=TRUE)/mean(x, na.rm=TRUE)), 
                   na.rm=TRUE)
      
      axis.limits[,which(conditions == i)]<- c(min(t$x), max(t$x), min(t$y),max(t$y))
    }
  }
  
  lim.x <- range(min(axis.limits[1,]), max(axis.limits[2,]))
  lim.y <- range(min(axis.limits[3,]), max(axis.limits[4,]))
  
  #par(mar = c(5, 5, 6, 3))
  plot(x = NULL
       , ylab ="Density"
       , xlab = "CV( log (intensity) )"
       , xlim = lim.x
       , ylim = lim.y
       , las=1
  )
  
  # density by condition
  conditions <- unique(conds)
  col.density = c(1:length(conditions))
  for (i in conditions){
    if (length(which(conds == i)) > 1){
      t <- apply(qData[,which(conds == i)], 1, 
                 function(x) 100*var(x, na.rm=TRUE)/mean(x, na.rm=TRUE))
      lines(density(t, na.rm = TRUE)
            , xlab=""
            , ylab=""
            , col=col.density[which(conditions == i)]
      )
    }
  }
  
  legend("topright"         
         , legend = conditions
         , col = col.density
         , pch = 15
         , bty = "n"
         , pt.cex = 2
         , cex = 1
         , horiz = FALSE
         , inset = c(0,0)
  )
  
}



#' Builds a densityplot of the CV of entities in the exprs() table
#' of a object. The CV is calculated for each condition present
#' in the dataset (see the slot \code{'Condition'} in the \code{pData()} table)
#' Same as the function \code{CVDistD} but uses the package \code{highcharter}
#' 
#' @title Distribution of CV of entities
#' 
#' @param qData A dataframe that contains quantitative data.
#' 
#' @param conds A vector of the conditions (one condition per sample).
#' 
#' @param palette xxx
#' 
#' @return A density plot
#' 
#' @author Samuel Wieczorek
#' 
#' @seealso \code{\link{densityPlotD}}.
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
#' CVDistD_HC(Biobase::exprs(Exp1_R25_pept), conds)
#' 
#' @import highcharter
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats density var
#' 
#' @export
#' 
CVDistD_HC <- function(qData, conds=NULL, palette = NULL){
  
  if (is.null(conds)) {
    warning("The vector of conditions is empty. The plot cannot be drawn.")
    return(NULL)}
  
  conditions <- unique(conds)
  n <- length(conditions)
  
  if (is.null(palette)){
    #palette <- RColorBrewer::brewer.pal(length(unique(conds)),"Dark2")[1:n]
    palette <- grDevices::colorRampPalette(brewer.pal(8, "Dark2"))(length(conditions))
  }else{
    if (length(palette) != ncol(qData)){
      warning("The color palette has not the same dimension as the number of samples")
      palette <- grDevices::colorRampPalette(brewer.pal(8, "Dark2"))(length(conditions))
      #return(NULL)
    }
  }
  
  
  
  # nbSeries = n
  # series <- list()
  # for (i in 1:length(conditions)){
  #     if (length(which(conds == conditions[i])) > 1){
  #         t <- apply(qData[,which(conds == conditions[i])], 1, 
  #                    function(x) 100*var(x, na.rm=TRUE)/mean(x, na.rm=TRUE))
  #         tmp <- data.frame(x = density(t, na.rm = TRUE)$x,
  #                           y = density(t, na.rm = TRUE)$y)
  #         series[[i]] <- list(name = conditions[i],
  #                             data = list_parse(tmp))
  #     }
  # }
  
  h1 <-  highchart() %>% 
    my_hc_chart(chartType = "spline", zoomType="x") %>%
    hc_colors(unique(palette)) %>%
    hc_legend(enabled = TRUE) %>%
    hc_xAxis(title = list(text = "CV(log(Intensity))")) %>%
    hc_yAxis(title = list(text = "Density")) %>%
    hc_tooltip(headerFormat= '',
               pointFormat = "<b>{series.name}</b>: {point.y} ",
               valueDecimals = 2) %>%
    my_hc_ExportMenu(filename = "logIntensity") %>%
    hc_plotOptions(
      series=list(
        connectNulls= TRUE,
        marker=list(
          enabled = FALSE)
      )
    )
  
  minX <- maxX <- 0
  maxY <- 0
  for (i in 1:n){
    if (length(which(conds == conditions[i])) > 1){
      t <- apply(qData[,which(conds == conditions[i])], 1, 
                 function(x) 100*var(x, na.rm=TRUE)/mean(x, na.rm=TRUE))
      tmp <- data.frame(x = density(t, na.rm = TRUE)$x,
                        y = density(t, na.rm = TRUE)$y)
      
      ymaxY <- max(maxY,tmp$y)
      xmaxY <- tmp$x[which(tmp$y==max(tmp$y))]
      minX <- min(minX, tmp$x)
      maxX <- max(maxX, 10*(xmaxY-minX))
      
      
      h1 <- h1 %>% hc_add_series(data=tmp, name=conditions[i]) }
  }
  
  h1 <- h1 %>%
    hc_chart(
      events = list(
        load = JS(paste0("function(){
                         var chart = this;
                         this.xAxis[0].setExtremes(",minX,",",maxX, ");
                         this.showResetZoom();}"))
      )
    )
  
  return(h1)
  
}
