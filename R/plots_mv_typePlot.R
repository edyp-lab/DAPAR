



#' This method is a wrapper for the function \code{\link{hc_mvTypePlot2}} adapted to objects
#' of class \code{MSnSet}).

#' @title Distribution of observed values with respect to intensity values 
#' from a \code{MSnSet} object
#' @param obj An object of class \code{MSnSet}.
#' @param ... See \code{\link{hc_mvTypePlot2}} 
#' @return A scatter plot
#' @author Florence Combes, Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' wrapper.hc_mvTypePlot2(Exp1_R25_pept)
#' @export
wrapper.hc_mvTypePlot2 <- function(obj,...){
    qData <- Biobase::exprs(obj)
    conds <- Biobase::pData(obj)[,"Condition"]
    hc_mvTypePlot2(qData, conds = conds,...)
}




#' This method shows density plots which represents the repartition of
#' Partial Observed Values for each replicate in the dataset.
#' The colors correspond to the different conditions (slot Condition in in the
#' dataset of class \code{MSnSet}).
#' The x-axis represent the mean of intensity for one condition and one
#' entity in the dataset (i. e. a protein) 
#' whereas the y-axis count the number of observed values for this entity
#' and the considered condition.
#' 
#' @title Distribution of Observed values with respect to intensity values
#' @param qData A dataframe that contains quantitative data.
#' @param conds A vector of the conditions (one condition per sample).
#' @param palette The different colors for conditions
#' @param typeofMV xxx
#' @param title The title of the plot
#' @return Density plots
#' @author Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' qData <- Biobase::exprs(Exp1_R25_pept)
#' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
#' hc_mvTypePlot2(qData, conds, title="POV distribution")
#' @export
#' @import highcharter
hc_mvTypePlot2 <- function(qData, conds, palette = NULL, typeofMV=NULL, title=NULL){
  if (is.null(conds)){return(NULL)}
  
  if (is.null(palette)){
    palette <- RColorBrewer::brewer.pal(length(unique(conds)),"Dark2")[1:length(unique(conds))]
  }else{
    if (length(palette) != length(unique(conds))){
      warning("The color palette has not the same dimension as the number of conditions")
      return(NULL)
    }
  }
  
  conditions <- conds
  mTemp <- nbNA <- nbValues <- matrix(rep(0,nrow(qData)*length(unique(conditions))), nrow=nrow(qData),
                                      dimnames=list(NULL,unique(conditions)))
  dataCond <- data.frame()
  ymax <- 0
  series <- list()
  myColors <- NULL
  j <- 1 
  
  for (iCond in unique(conditions)){
    if (length(which(conditions==iCond)) == 1){
      
      mTemp[,iCond] <- qData[,which(conditions==iCond)]
      nbNA[,iCond] <- as.integer(is.OfType(qData[,which(conditions==iCond)]))
      nbValues[,iCond] <- length(which(conditions==iCond)) - nbNA[,iCond]
    } else {
      mTemp[,iCond] <- apply(qData[,which(conditions==iCond)], 1, mean, na.rm=TRUE)
      nbNA[,iCond] <- apply(qData[,which(conditions==iCond)],1,function(x) length(which(is.na(x) == TRUE)))
      nbValues[,iCond] <- length(which(conditions==iCond)) - nbNA[,iCond]
    }
    
    
    for (i in 1:length(which(conditions==iCond))){
      data <- mTemp[which(nbValues[, iCond] == i), iCond]
      tmp <- NULL    
      if (length(data) >= 2)
      {
        tmp <- density(mTemp[which(nbValues[,iCond]==i),iCond])
        tmp$y <- tmp$y + i
        if (max(tmp$y) > ymax) { ymax <- max(tmp$y)}
      }
      series[[j]] <- tmp
      myColors <- c(myColors, palette[which(unique(conditions)==iCond)])
      j <- j+1
    }
    
  }
  
  
  hc <-  highchart() %>%
    hc_title(text = title) %>%
    dapar_hc_chart(chartType = "spline", zoomType="xy") %>%
    
    hc_legend(align = "left", verticalAlign = "top",
              layout = "vertical") %>%
    hc_xAxis(title = list(text = "Mean of intensities")) %>%
    hc_yAxis(title = list(text = "Number of quantity values per condition"),
             #categories = c(-1:3)
             #min = 1, 
             # max = ymax,
             tickInterval= 0.5
    ) %>%
    # hc_colors(palette) %>%
    hc_tooltip(headerFormat= '',
               pointFormat = "<b> {series.name} </b>: {point.y} ",
               valueDecimals = 2) %>%
    dapar_hc_ExportMenu(filename = "POV_distribution") %>%
    hc_plotOptions(
      series=list(
        showInLegend = TRUE,
        animation=list(
          duration = 100
        ),
        connectNulls= TRUE,
        marker=list(
          enabled = FALSE)
        
      )
    )
  
  for (i in 1:length(series)){
    hc <- hc_add_series(hc,
                        data = list_parse(data.frame(cbind(x = series[[i]]$x, 
                                                           y = series[[i]]$y))), 
                        showInLegend=FALSE,
                        color = myColors[i],
                        name=conds[i])
  }
  
  # add three empty series for the legend entries. Change color and marker symbol
  for (c in 1:length(unique(conds))){
    hc <-  hc_add_series(hc,data = data.frame(),
                         name = unique(conds)[c],
                         color = palette[c],
                         marker = list(symbol = "circle"),
                         type = "line")
  }
  
  hc
  return(hc)
}
