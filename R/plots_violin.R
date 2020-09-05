
#' @title Builds a violinplot from a dataframe
#' 
#' @param obj xxx
#' 
#' @param legend A vector of the conditions (one condition per sample).
#' 
#' @param palette xxx
#' 
#' @return A violinplot
#' 
#' @author Samuel Wieczorek
#' 
#' @seealso \code{\link{densityPlotD}}
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' library(vioplot)
#' legend <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
#' conds <- Biobase::pData(Exp1_R25_pept)[["Condition"]]
#' key <- "Protein_group_IDs"
#' violinPlotD(exprs(Exp1_R25_pept), conds=legend, keyId=key, palette=c(rep('blue',3), rep('green',3)),subset.view=1:10)
#' 
#' @importFrom vioplot vioplot
#' 
#' @importFrom grDevices colorRampPalette
#' 
#' @importFrom graphics plot.window axis mtext legend points segments plot.new
#' 
#' @importFrom RColorBrewer brewer.pal
#' 
#' @export
#' 
violinPlotD <- function(qData, conds, keyId, legend=NULL, palette = NULL, subset.view=NULL){
  
  
  graphics::plot.new()
  
  if (is.null(qData)) {
    warning('The dataset in NULL and cannot be shown')
    return(NULL)
  }
  
  
  if(missing(conds))
    stop("'conds' is missing.")
  
  if (is.null(legend)) {
    legend <- conds
    for (i in unique(conds))
      legend[which(conds==i)] <- paste0(i,'_', 1:length(which(conds==i)))
  }
  
  if (!is.null(subset.view)) {
    if (is.null(keyId)|| missing(keyId))
      stop("'keyId' is missing.")
  }
  
  
  #palette <- BuildPalette(conds, palette)
  if (is.null(palette)){palette <- rep("#FFFFFF", ncol(qData))
  } else {
    if (length(palette) != ncol(qData)){
      warning("The color palette has not the same dimension as the number of samples")
      return(NULL)
    }
  }
  
  graphics::plot.window(xlim=c(0,ncol(qData)+1),
                        ylim=c(min(na.omit(qData)),max(na.omit(qData)))
  )
  title( ylab="Log (intensity)")
  
  for (i in 1:ncol(qData)) {
    vioplot::vioplot(na.omit(qData[,i]), col = palette[i], add=TRUE, at=i)
  }
  
  
  graphics::axis(2, yaxp = c(floor(min(na.omit(qData))), 
                             floor(max(na.omit(qData))), 5), las=1)
  
  if( !is.null(legend) ) {
    if ( is.vector(legend) ) {
      N <- 1
    } else { 
      N <- ncol(legend)
    }
    
    for (i in 1:N) {
      graphics::axis(side=1,
                     at = 1:ncol(qData),
                     labels = if (is.vector(legend) ) {legend} else {legend[,i]},
                     line= 2*i-1
      )
    }
    
    graphics::mtext("Samples",side=1,line=6+length(colnames(legend)), cex.lab=1, las=1)
  }
  
  # Display of rows to highlight (index of row in subset.view) 
  if(!is.null(subset.view)){
    idVector <- keyId
    pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(subset.view))
    
    n=0
    for (i in subset.view) {
      n=n+1
      for (c in 1:(ncol(qData)-1)) {
        graphics::segments(y0=qData[i,c],y1=qData[i,c+1],x0=c,x1=c+1,pch=16,col=pal[n],lwd=2)
        graphics::points(y=qData[i,c],x=c,pch=16,col=pal[n])
      }
      graphics::points(y=qData[i,ncol(qData)],x=ncol(qData),pch=16,col=pal[n])
    }
    graphics::legend("topleft",legend=idVector[subset.view],lty=1,lwd=2,col=pal,pch=16,bg="transparent",bty="n")
  }
  
}



# 
# violinPlotD <- function(obj, legend=NULL, palette = NULL){
#   plot.new()
#   qData <- Biobase::exprs(obj)
#   
#   if (!is.null(palette)) {
#     if (length(palette) != ncol(qData)){
#       warning("The color palette has not the same dimension as the number of samples")
#       return(NULL)
#     }
#   } else {
#     palette <- rep('#FFFFFF',ncol(qData))
#   }
#   
#   
#   plot.window(xlim=c(0,ncol(qData)+1),
#               ylim=c(min(na.omit(qData)),max(na.omit(qData))))
#   title( ylab="Log (intensity)")
#   
#   for (i in 1:ncol(qData)) {
#     vioplot(na.omit(qData[,i]), col = palette[i], add=TRUE, at=i)}
#   
#   
#   axis(2, yaxp = c(floor(min(na.omit(qData))), 
#                    floor(max(na.omit(qData))), 5), las=1)
#   
#   if( !is.null(legend)) {
#     if (is.vector(legend) ){
#       N <- 1} else{ N <- ncol(legend)}
#     
#     for (i in 1:N){
#       axis(side=1,
#            at = 1:ncol(qData),
#            label = if (is.vector(legend) ) 
#            {legend} else {legend[,i]},
#            line= 2*i-1
#       )
#     }
#     
#     mtext("Samples",side=1,line=6+length(colnames(legend)), cex.lab=1, las=1)
#   }
#   
# }
