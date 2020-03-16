



#' ViolinPlot for quantitative proteomics data
#' 
#' @title Builds a violinplot from a dataframe
#' @param obj xxx
#' @param legend A vector of the conditions (one condition per sample).
#' @param palette xxx
#' @param subset.view A vector of index indicating rows to highlight
#' @return A violinplot
#' @author Samuel Wieczorek, Anais Courtier
#' @seealso \code{\link{densityPlotD}}
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' library(vioplot)
#' legend <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
#' violinPlotD(Exp1_R25_pept, legend=legend,subset.view=20:30)
#' @importFrom Biobase exprs
#' @importFrom vioplot vioplot
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics plot.window
#' @export
violinPlotD <- function(obj, legend=NULL, palette = NULL,subset.view=NULL){
  plot.new()
  qData <- Biobase::exprs(obj)
  
  if (!is.null(palette)) {
    if (length(palette) != ncol(qData)){
      warning("The color palette has not the same dimension as the number of samples")
      return(NULL)
    }
  } else {
    palette <- rep('#FFFFFF',ncol(qData))
  }
  
  
  graphics::plot.window(xlim=c(0,ncol(qData)+1),
              ylim=c(min(na.omit(qData)),max(na.omit(qData))))
  title( ylab="Log (intensity)")
  
  for (i in 1:ncol(qData)) {
    vioplot::vioplot(na.omit(qData[,i]), col = palette[i], add=TRUE, at=i)}
  
  
  axis(2, yaxp = c(floor(min(na.omit(qData))), 
                   floor(max(na.omit(qData))), 5), las=1)
  
  if( !is.null(legend))
  {
    if (is.vector(legend) ){
      N <- 1} else{ N <- ncol(legend)}
    
    for (i in 1:N){
      axis(side=1,
           at = 1:ncol(qData),
           label = if (is.vector(legend) ) 
           {legend} else {legend[,i]},
           line= 2*i-1
      )
    }
    
    mtext("Samples",side=1,line=6+length(colnames(legend)), cex.lab=1, las=1)
  }
  # Display of rows to highlight (index of row in subset.view) 
  if(!is.null(subset.view)){
    idColName<-obj@experimentData@other$proteinId
    idVector=obj@featureData@data[,idColName]
    pal=grDevices::colorRampPalette(brewer.pal(8, "Set1"))(length(subset.view))
    
    n=0
    for (i in subset.view) {
      n=n+1
      for (c in 1:(ncol(qData)-1)) {
        segments(y0=qData[i,c],y1=qData[i,c+1],x0=c,x1=c+1,pch=16,col=pal[n],lwd=2)
        points(y=qData[i,c],x=c,pch=16,col=pal[n])
      }
      points(y=qData[i,ncol(qData)],x=ncol(qData),pch=16,col=pal[n])
    }
    legend("topleft",legend=idVector[subset.view],lty=1,lwd=2,col=pal,pch=16,bg="transparent",bty="n")
  }
  
}



