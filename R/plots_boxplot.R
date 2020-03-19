#' Boxplot for quantitative proteomics data
#' 
#' @title Builds a boxplot from a dataframe
#' @param obj xxx
#' @param conds xxx
#' @param legend A vector of the conditions (one string per sample).
#' @param palette xxx
#' @return A boxplot
#' @author Florence Combes, Samuel Wieczorek
#' @seealso \code{\link{densityPlotD}}
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
#' boxPlotD(Exp1_R25_pept, conds)
#' @importFrom Biobase exprs
#' @importFrom RColorBrewer brewer.pal
#' @import graphics
#' @export
boxPlotD <- function(obj,conds, legend=NULL,palette=NULL){
  qData <- Biobase::exprs(obj)
  if (is.null(palette)){
    pal <- RColorBrewer::brewer.pal(length(unique(conds)),"Dark2")[1:length(unique(conds))]
    
    for (i in 1:ncol(qData)){
      palette[i] <- pal[ which(conds[i] == unique(conds))]
    }
    
  }else{
    if (length(palette) != ncol(qData)){
      warning("The color palette has not the same dimension as the number of samples")
      return(NULL)
    }
  }
  
  
  boxplot(qData
          ,las = 1
          , col = palette
          , cex = 2
          , axes=TRUE
          , xaxt = "n"
          , ylab = "Log (intensity)"
          , pt.cex = 4
          , horizontal = FALSE
  )
  
  
  if( is.null(legend)){legend <- Biobase::pData(obj)$Condition}
  axis(side=1,at = 1:ncol(qData), label = legend)
  #mtext("Samples", side=1,  line=(6+ncol(legend)), cex.lab=1, las=1)
  
  abline(h=0) 
  
}




#' Boxplot for quantitative proteomics data using the library \code{highcharter}
#' 
#' @title Builds a boxplot from a dataframe using the library \code{highcharter}
#' @param obj xxx
#' @param legend A vector of the conditions (one condition per sample).
#' @param palette xxx
#' @param subset.view A vector of index indicating rows to highlight
#' @return A boxplot
#' @author Samuel Wieczorek, Anais Courtier
#' @seealso \code{\link{densityPlotD_HC}}
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' legend <- Biobase::pData(Exp1_R25_pept)[,"Sample.name"]
#' boxPlotD_HC(Exp1_R25_pept, legend, subset.view=1:10)
#' @importFrom Biobase exprs pData
#' @import highcharter
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export
boxPlotD_HC <- function(obj, legend=NULL, palette = NULL,subset.view=NULL){
  
  if (is.null(obj)){
    warning('The dataset in NULL and cannot be shown')
    return(NULL)
  }
  
  qData <- Biobase::exprs(obj)
  if( is.null(legend)){legend <- Biobase::pData(obj)[,"Sample.name"]}
  
  palette <- BuildPalette(Biobase::pData(obj)$Condition, palette)
  
  bx <-boxplot(qData, na.rm=TRUE)
  df_outlier <- data.frame(x=bx$group-1,
                           y = bx$out)
  
  #tmp <- NULL
  #for (i in 1:ncol(qData)){
  #  tmp <- c(tmp, rep(paste("S",i,legend[i], collapse="_"),nrow(qData)))
  #}
  tmp <- NULL
  for (i in 1:ncol(qData)){
    tmp <- c(tmp, rep(paste(paste0(rep("A", i), collapse=""),legend[i], sep='_'),nrow(qData)))
  }
  
  df <- data.frame(values = as.vector(qData,mode='numeric'),samples = tmp, stringsAsFactors = FALSE)
  
  
  hc <- highcharter::hcboxplot(x=df$values, var = df$samples, colorByPoint = TRUE, outliers = TRUE) %>%
    highcharter::hc_chart(type="column") %>%
    highcharter::hc_yAxis(title = list(text = "Log (intensity)")) %>%
    highcharter::hc_xAxis(title = list(text = "Samples"), categories=legend) %>%
    highcharter::hc_colors(palette) %>%
    highcharter::hc_add_series(type= "scatter",df_outlier,name="Outliers",tooltip=list(enabled=F,headerFormat ="",pointFormat="{point.y: .2f} ")) %>%  
    highcharter::hc_plotOptions(
      
      boxplot= list(
        
        fillColor= c('lightgrey'),
        lineWidth= 3,
        medianColor= 'grey',
        medianWidth= 3,
        stemColor= '#A63400',
        stemDashStyle= 'dot',
        stemWidth= 1,
        whiskerColor= '#3D9200',
        whiskerLength= '20%',
        whiskerWidth= 3
      ),
      scatter = list(
        marker=list(
          fillColor = 'white',
          lineWidth = 0.5,
          lineColor = 'grey'
        )
      )
    )
  
  # Display of rows to highlight (index of row in subset.view) 
  if(!is.null(subset.view)){
    idVector <- Biobase::fData(obj)[,keyId(obj)]
    pal=grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(subset.view))    
    n=0
    for(i in subset.view){
      n=n+1
      dfSubset <- data.frame(y = as.vector(qData[i,],mode='numeric'), x = as.numeric(factor(names(qData[i,])))-1, stringsAsFactors = FALSE)
      hc<-hc %>%
        highcharter::hc_add_series(type= "line",data=dfSubset,color=pal[n], dashStyle = "shortdot",name=idVector[i],
                      tooltip=list(enabled=T,headerFormat ="",pointFormat="{point.series.name} : {point.y: .2f} ") )
      
    }
  }  
  
  hc
  
}


