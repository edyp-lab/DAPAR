#' This function show the density plots of Fold Change (the same as calculated 
#' by limma) for a list of the comparisons of conditions in a differential 
#' analysis.
#' 
#' @title Density plots of logFC values
#' 
#' @param df_logFC A dataframe that contains the logFC values
#' 
#' @param threshold_LogFC The threshold on log(Fold Change) to
#' distinguish between differential and non-differential data 
#' 
#' @param pal xxx
#' 
#' @return A highcharts density plot
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_prot, package='DAPARdata')
#' obj <- Exp1_R25_prot[1:1000]
#' level <- obj@experimentData@other$typeOfData
#' metacell.mask <- match.metacell(GetMetacell(obj), 'missing', level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op='>=', th=1)
#' obj <- MetaCellFiltering(obj, indices, cmd='delete')
#' qData <- Biobase::exprs(obj$new)
#' sTab <- Biobase::pData(obj$new)
#' res <- limmaCompleteTest(qData, sTab, comp.type = "OnevsAll")
#' hc_logFC_DensityPlot(res$logFC, threshold_LogFC=1)
#' pal <- ExtendPalette(2, 'Dark2')
#' hc_logFC_DensityPlot(res$logFC, threshold_LogFC=1, pal=pal)
#' 
#' @export
#' 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' 
hc_logFC_DensityPlot <-function(df_logFC, 
                                threshold_LogFC = 0, 
                                pal = NULL){
  if (is.null(df_logFC) || threshold_LogFC < 0){
    hc <- NULL
    return(NULL)
  }
  
  myColors <- NULL
  if (is.null(pal)){
    warning("Color palette set to default.")
    myColors <- ExtendPalette(ncol(df_logFC), "Paired")
  } else {
    if (length(pal) != ncol(df_logFC)){
      warning("The color palette has not the same dimension as the number of samples")
      myColors <- ExtendPalette(pal, "Paired")
    }
    myColors <- pal
  }
  
  nValues <- nrow(df_logFC)*ncol(df_logFC)
  nInf <- length(which(df_logFC <= -threshold_LogFC))
  nSup <- length(which(df_logFC >= threshold_LogFC))
  nInside <- length(which(abs(df_logFC) < threshold_LogFC))
  hc <- NULL
  hc <-  highcharter::highchart() %>% 
    hc_title(text = "log(FC) repartition") %>% 
    my_hc_chart(chartType = "spline", zoomType="x") %>%
    hc_legend(enabled = TRUE) %>%
    hc_colors(myColors) %>%
    hc_xAxis(title = list(text = "log(FC)"),
             plotBands = list(list(from= -threshold_LogFC, to = threshold_LogFC, color = "lightgrey")),
             plotLines=list(list(color= "grey" , width = 2, value = 0, zIndex = 5))
    )%>%
    hc_yAxis(title = list(text="Density")) %>%
    hc_tooltip(headerFormat= '',
               pointFormat = "<b> {series.name} </b>: {point.y} ",
               valueDecimals = 2) %>%
    my_hc_ExportMenu(filename = "densityplot") %>%
    hc_plotOptions(
      series=list(
        animation=list(duration = 100),
        connectNulls= TRUE,
        marker=list(enabled = FALSE)
      )
    )
  
  maxY.inf <- NULL
  maxY.inside <- NULL
  maxY.sup <- NULL
  minX <- NULL
  maxX<- NULL
  
  
  for (i in 1:ncol(df_logFC)){
    tmp <- density(df_logFC[,i])
    ind <- tmp$y[which(tmp$x <= -threshold_LogFC)]
    maxY.inf <- max(maxY.inf, ifelse(length(ind)==0,0,ind))
    maxY.inside <- max(maxY.inf, tmp$y[intersect(which(tmp$x > -threshold_LogFC),which(tmp$x < threshold_LogFC))])
    ind <- tmp$y[which(tmp$x > threshold_LogFC)]
    maxY.sup <- max(maxY.sup, ifelse(length(ind)==0,tmp$y[length(tmp$y)],ind))
    minX <- min(minX,tmp$x)
    maxX <- max(maxX,tmp$x)
    
    
    hc <- hc_add_series(hc,
                        data.frame(x = tmp$x,  y = tmp$y), 
                        name=colnames(df_logFC)[i])
  }
  
  ## add annotations
  if(threshold_LogFC > 0) {
    hc <- hc %>% hc_add_annotation(
      labelOptions = list(
        shape='connector',
        backgroundColor = 'lightgrey',
        #verticalAlign = 'bottom',
        align='left',
        #distance=0,
        style=list(
          fontSize= '1.5em',
          textOutline= '1px white'
        ),
        borderWidth = 0,
        x = 20
      ),
      labels = list(
        list(
          point = list(
            xAxis = 0,
            yAxis = 0,
            x = 0,
            y = maxY.inside
          ),
          text = paste0("n Filtered out = ",nInside, "<br>(", round(100*nInside/nValues, digits=2), "%)")
        )
      )
    )
  }
  if (threshold_LogFC >= minX){
    hc <- hc %>%
      hc_add_annotation(
        labelOptions = list(
          shape='connector',
          backgroundColor = 'rgba(255,255,255,0.5)',
          verticalAlign = 'top',
          borderWidth = 0,
          crop=TRUE,
          style=list(
            color = 'blue',
            fontSize= '1.5em',
            textOutline= '1px white'
          ),
          y = -10
        ),
        labels = list(
          list(
            point = list(
              xAxis = 0,
              yAxis = 0,
              x = mean(c(minX,-threshold_LogFC)),
              y = maxY.inf
            ),
            text = paste0("nInf = ",nInf, "<br>(", round(100*nInf/nValues, digits=2), ")%")
          )
        )
      )
  }
  
  if (threshold_LogFC <= maxX){
    hc <- hc %>% hc_add_annotation(
      labelOptions = list(
        shape='connector',
        backgroundColor = 'blue',
        verticalAlign = 'top',
        borderWidth = 0,
        style=list(
          color = 'blue',
          fontSize= '1.5em',
          textOutline= '1px white'
        ),
        y = -5
      ),
      labels = list(
        list(
          point = list(
            xAxis = 0,
            yAxis = 0,
            x =  mean(c(maxX,threshold_LogFC)),
            y = maxY.sup
          ),
          text = paste0("nSup = ",nSup, "<br>(", round(100*nSup/nValues, digits=2), ")%")
        )
      )
    )
    
  }
  
  
  
  return(hc)
  
}


