##' This function is a wrapper for using the boxPlotD function with objects of class \code{\link{MSnSet}}
##' 
##' @title Wrapper to the boxplotD function on an object \code{\link{MSnSet}}
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param dataForXAxis A vector of strings containing the names of columns 
##' in \code{pData()} to print labels on X-axis (Default is "Label").
##' @return A boxplot
##' @author Florence Combes, Samuel Wieczorek
##' @seealso \code{\link{wrapper.densityPlotD}}
##' @examples
##' data(UPSprotx2)
##' types <- c("Label","Analyt.Rep")
##' wrapper.boxPlotD(UPSprotx2, types)
wrapper.boxPlotD <- function(obj, dataForXAxis="Label"){
  
  qData <- exprs(obj)
  dataForXAxis <- as.matrix(pData(obj)[,dataForXAxis])
  labels <- pData(obj)[,"Label"]
  
  boxPlotD(qData, dataForXAxis, labels)
  
}



##' Boxplot for quantitative proteomics data
##' 
##' @title Builds a boxplot from a dataframe
##' @param qData A dataframe that contains quantitative data.
##' @param dataForXAxis A vector containing the types of replicates 
##' to use as X-axis. Available values are: Label, Analyt.Rep,
##' Bio.Rep and Tech.Rep. Default is "Label".
##' @param labels A vector of the conditions (labels) (one label per sample).
##' @return A boxplot
##' @author Florence Combes, Samuel Wieczorek
##' @seealso \code{\link{densityPlotD}}
##' @examples
##' data(UPSprotx2)
##' qData <- exprs(UPSprotx2)
##' types <- c("Label","Analyt.Rep")
##' dataForXAxis <- pData(UPSprotx2)[,types]
##' labels <- pData(UPSprotx2)[,"Label"]
##' boxPlotD(qData, dataForXAxis, labels)
boxPlotD <- function(qData, dataForXAxis=NULL, labels=NULL){
  
  pal <- getPaletteForLabels(labels)
  par(oma = c(2+length(colnames(dataForXAxis)), 0, 0, 0))
  
  boxplot(qData
          , las = 1
          , col = pal
          , cex = 2
          , axes=TRUE
          , xaxt = "n"
          , ylab = "Log (intensity)"
          , pt.cex = 4
          , horizontal = FALSE
  )
  
  if( !is.null(dataForXAxis))
    {
    for (i in 1:ncol(dataForXAxis)){
    axis(side=1,
         at = 1:ncol(qData),
         labels = dataForXAxis[,i],
         line= 2*i-1
    )
  }
  
 mtext("Samples", side=1, line=6+length(colnames(dataForXAxis)), cex.lab=1, las=1)
}
  
 abline(h=0) 
  palette("default")
}

##' Plot to compare the quantitative proteomics data before and after normalization
##' 
##' @title Builds a plot from a dataframe
##' @param qDataBefore A dataframe that contains quantitative data before normalization.
##' @param qDataAfter A dataframe that contains quantitative data after normalization.
##' @param dataForXAxis A vector containing the types of replicates 
##' to use as X-axis. Available values are: Label, Analyt.Rep,
##' Bio.Rep and Tech.Rep. Default is "Label".
##' @param labels A vector of the conditions (labels) (one label per sample).
##' @return A plot
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' data(UPSprotx2)
##' qDataBefore <- exprs(UPSprotx2)
##' labels <- pData(UPSprotx2)[,"Label"]
##' qDataAfter <- normalizeD(qDataBefore, labels, "Median Centering", "within conditions")
##' types <- c("Label","Analyt.Rep")
##' dataForXAxis <- pData(UPSprotx2)[,types]
##' compareNormalizationD(qDataBefore, qDataAfter, dataForXAxis, labels)
compareNormalizationD <- function(qDataBefore, qDataAfter, dataForXAxis=NULL, labels=NULL){
  #requireNamespace(scales)
  pal <- getPaletteForLabels(labels)
  par(oma = c(2+length(colnames(dataForXAxis)), 0, 0, 0))
  
  plot(qDataBefore, 
       qDataAfter,
          , las = 1
          , col = alpha(pal, 0.5)
          , cex = 0.5
          , axes=TRUE
       , pch = 16
  )
  

  palette("default")
}




##' This function is a wrapper for using the densityPlotD function with objects of class \code{\link{MSnSet}}
##' 
##' @title Builds a densityplot from an object of class \code{\link{MSnSet}}
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param lab2Show A vector of labels to show in densityplot.
##' @param highLightLabel The name of the Label to highlight
##' in the density plot.
##' @return A density plot
##' @author Alexia Dorffer
##' @seealso \code{\link{wrapper.boxPlotD}}, \code{\link{wrapper.varianceDistD}}
##' @examples data(UPSprotx2)
##' wrapper.densityPlotD(UPSprotx2)
wrapper.densityPlotD <- function(obj, lab2Show=NULL,  highLightLabel=NULL){
  qData <- exprs(obj)
  labels <- pData(obj)[,"Label"]
  densityPlotD(qData, labels, lab2Show,  highLightLabel)
}




##' Densityplot of quantitative proteomics data over samples.
##' 
##' @title Builds a densityplot from a dataframe
##' @param qData A dataframe that contains quantitative data.
##' @param labels A vector of the conditions (labels) (one label per sample).
##' @param lab2Show A vector of labels to show in densityplot. If NULL, then all labels are displayed.
##' @param highLightLabel The name of the Label to highlight
##' in the density plot.
##' @return A density plot
##' @author Florence Combes, Samuel Wieczorek
##' @seealso \code{\link{boxPlotD}}, \code{\link{varianceDistD}}
##' @examples data(UPSprotx2)
##' qData <- exprs(UPSprotx2)
##' labels <- lab2Show <- pData(UPSprotx2)[,"Label"]
##' densityPlotD(qData, labels, lab2Show)
densityPlotD <- function(qData, labels=NULL,lab2Show=NULL, highLightLabel=NULL){
    
  if (is.null(lab2Show)){
    lab2Show <- unique(labels)
  }
  
  if (is.null(labels)) return(NULL)
  
  indices <- which(labels %in% lab2Show)
  n <- length(indices)
  
  axis.limits <- matrix(data = 0, nrow = 4, ncol = ncol(qData))
  
  for (i in 1:ncol(qData)){
    dens <- density(qData[,i], na.rm = TRUE)
    axis.limits[,i] <- c(min(dens$x), max(dens$x), min(dens$y), max(dens$y))
  }
  lim.x <- range(min(axis.limits[1,]), max(axis.limits[2,]))
  lim.y <- range(min(axis.limits[3,]), max(axis.limits[4,]))
  
  ##Colors definition
  nColors <- length(unique(labels))
  col.density <- c(1:nColors)
  pal <- getPaletteForLabels(labels)
  
  lineWD <- NULL
  lineWD <- c(rep(1, length(colnames(qData))))
  if (!is.null(highLightLabel)) {
    lineWD[which(labels == highLightLabel)] <- 3
  }
  
  plot(x=NULL
       , ylab ="Density"
       , xlab = "log(intensity)"
       ,xlim = lim.x
       ,ylim = lim.y
       ,las = 1
       ,cex.lab = 1
       ,cex.axis = 1
       ,cex.main = 3)
  
  for (i in indices){
    lines(density(qData[,i], na.rm=TRUE)
          , col = col.density[which(labels[i] == unique(labels))]
          , lwd = lineWD[i])
  }
  
  legend("topright"         
         , legend = unique(labels)
         , col = col.density
         , pch = 15 
         , bty = "n"
         , pt.cex = 2
         , cex = 1
         , horiz = FALSE
         , inset=c(0,0)
  )
}



##' Builds a densityplot of the variance of entities in the exprs() table
##' of an object \code{\link{MSnSet}}. The variance is calculated for each condition (Label) present
##' in the dataset (see the slot \code{'Label'} in the \code{pData()} table).
##' 
##' @title Distribution of variance of proteins
##' @param obj An object of class \code{\link{MSnSet}}.
##' @return A density plot
##' @author Alexia Dorffer
##' @seealso \code{\link{wrapper.densityPlotD}}
##' @examples data(UPSprotx2)
##' wrapper.varianceDistD(UPSprotx2)
wrapper.varianceDistD <- function(obj){
  qData <- exprs(obj)
  labels <- pData(obj)[,"Label"]
  varianceDistD(qData, labels)
}




##' Builds a densityplot of the variance of entities in the exprs() table
##' of a object. The variance is calculated for each condition (Label) present
##' in the dataset (see the slot \code{'Label'} in the \code{pData()} table)
##' 
##' @title Distribution of variance of proteins
##' @param qData A dataframe that contains quantitative data.
##' @param labels A vector of the conditions (labels) (one label per sample).
##' @return A density plot
##' @author Florence Combes, Samuel Wieczorek
##' @seealso \code{\link{densityPlotD}}.
##' @examples data(UPSprotx2)
##' varianceDistD(UPSprotx2)
varianceDistD <- function(qData, labels=NULL){
    
  if (is.null(labels)) {return(NULL)}
  conditions <- unique(labels)
  n <- length(conditions)
  axis.limits <- matrix(data = 0, nrow = 4, ncol = n)
  for (i in conditions){
    t <- density(apply(qData[,which(labels == i)], 1, 
                       function(x) var(x, na.rm=TRUE)), na.rm=TRUE)
    axis.limits[,which(conditions == i)]<- c(min(t$x), max(t$x), min(t$y),
                                             max(t$y))
  }
  
  lim.x <- range(min(axis.limits[1,]), max(axis.limits[2,]))
  lim.y <- range(min(axis.limits[3,]), max(axis.limits[4,]))
  
  par(mar = c(5, 5, 6, 3))
  plot(x = NULL
       , ylab ="Density"
       , xlab = "Variance( log (intensity) )"
       , xlim = lim.x
       , ylim = lim.y
       , las=1
  )
  
  # density by condition
  pal <- getPaletteForLabels(labels)
  conditions <- unique(labels)
  col.density = c(1:length(conditions))
  for (i in conditions){
    t <- apply(qData[,which(labels == i)], 1, 
               function(x) var(x, na.rm = TRUE))
    lines(density(t, na.rm = TRUE)
          , xlab=""
          , ylab=""
          , col=col.density[which(conditions == i)]
    )
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


##' Builds a correlation matrix based on a \code{\link{MSnSet}} object.
##' 
##' @title Displays a correlation matrix of the quantitative data of the
##' \code{exprs()} table
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param rate xxxx
##' @param indLegend A vector of indices in the columns in
##' the \code{pData()} table choosen for the labels in the axes.
##' @return A colored correlation matrix
##' @author Alexia Dorffer
##' @examples data(UPSprotx2)
##' wrapper.corrMatrixD(UPSprotx2)
wrapper.corrMatrixD <- function(obj, rate=5, indLegend=NULL){
  qData <- exprs(obj)
  samplesData <- pData(obj)
  corrMatrixD(qData, samplesData, indLegend, rate)
}




##' Correlation matrix based on a \code{\link{MSnSet}} object
##' 
##' @title Displays a correlation matrix of the quantitative data of the
##' \code{exprs()} table.
##' @param qData A dataframe of quantitative data.
##' @param samplesData A dataframe where lines correspond to samples and columns to the meta-data for those samples.
##' @param indLegend A vector of indices in the columns of
##' the samplesData table choosen for the labels in the axes.
##' @param rate The rate parameter to control de exponential law for the gradient of colors
##' @return A colored correlation matrix
##' @author Florence Combes, Samuel Wieczorek
##' @examples data(UPSprotx2)
##' qData <- exprs(UPSprotx2)
##' samplesData <- pData(UPSprotx2)
##' corrMatrixD(qData, samplesData)
corrMatrixD <- function(qData, samplesData, indLegend=NULL, gradientRate = 5){
  Var1 <- Var2 <- value <- NULL
  if (is.null(indLegend)) {indLegend <- colnames(samplesData[-1])} 
  
  nb <- length(colnames(qData))
  for (j in 1:nb){
    noms <- NULL
    
    for (i in indLegend){
      #if (pData(obj)[i,j] != "_"){
      noms <- paste(noms, samplesData[j,i], sep=" ")
      #}     
    }
    colnames(qData)[j] <- noms
  }
  
  z <- cor(qData,use = 'pairwise.complete.obs')
  text <- element_text(colour="black", size = 16, face = "bold")
  d <- qplot(x = Var1, 
             y = Var2, 
             data = melt(z), 
             fill = value, 
             geom = "tile") +
    theme(axis.text = element_text(size=16),
          axis.title = element_text(size=20, face="bold"),
          axis.text.x = element_text(angle=30, vjust=1, hjust=1),
          legend.text = text,
          legend.title = text) +
    labs(x = "", y = "") +
    #scale_fill_gradient2(midpoint=0.5, space="Lab",low = "white",mid = "white", high = "steelblue", limits=c(0.1,1))
  #scale_fill_gradientn(colours=cscale(seq(0.1,1,0.1), seq_gradient_pal("grey80", "black")),trans=exp_trans(), limits=c(0.1,1))
    scale_fill_gradientn (
    colours=colorRampPalette (c ("white", "lightblue","darkblue")) (101),
    values = c(pexp(seq(0,1,0.01), rate=gradientRate),1), limits=c(0,1))
  plot(d)
}


##' Builds a heatmap of the quantitative proteomic data of a \code{\link{MSnSet}} object.
##' 
##' @title This function is a wrapper to \code{\link{heatmap.2}} that displays 
##' quantitative data in the \code{exprs()} table of an object of
##' class \code{\link{MSnSet}}
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param distance The distance used by the clustering algorithm to compute 
##' the dendrogram. See \code{help(heatmap.2)}.
##' @param cluster the clustering algorithm used to build the dendrogram.
##' See \code{help(heatmap.2)}
##' @return A heatmap
##' @author Alexia Dorffer
##' @examples data(testWithoutNA)
##' wrapper.heatmapD(testWithoutNA)
wrapper.heatmapD  <- function(obj, distance="euclidean", cluster="average"){
  qData <- exprs(obj)
  heatmapD(qData, distance, cluster)
}




##' Heatmap of the quantitative proteomic data of a \code{\link{MSnSet}} object
##' 
##' @title This function is a wrapper to \code{\link{heatmap.2}} that displays 
##' quantitative data in the \code{exprs()} table of an object of
##' class \code{\link{MSnSet}}
##' @param qData A dataframe that contains quantitative data.
##' @param distance The distance used by the clustering algorithm to compute 
##' the dendrogram. See \code{help(heatmap.2)}
##' @param cluster the clustering algorithm used to build the dendrogram.
##' See \code{help(heatmap.2)}
##' @return A heatmap
##' @author Florence Combes, Samuel Wieczorek
##' @examples data(testWithoutNA)
##' qData <- exprs(testWithoutNA)
##' heatmapD(qData)
heatmapD <- function(qData, distance="euclidean", cluster="average"){
  ##Check parameters
  paramdist <- c("euclidean", "manhattan") 
  if (!(distance %in% paramdist)){
    stop("Param distance is not correct.")
    return (NULL)
  }
  
  paramcluster <- c("ward.D", "average")
  if (!(cluster %in%  paramcluster)){
    stop("Param clustering is not correct.")
    return (NULL)
  }
  
  if (getNumberOfEmptyLines(qData) != 0)  {
    stop("Your dataset contains empty lines.
         Please filter or impute missing values before.")
    return (NULL)
  }
  else {
    .data <- matrix(qData, 
                    ncol = ncol(qData), 
                    byrow = FALSE,
                    dimnames = list(rownames(qData), colnames(qData))
    )
    colors = c(seq(-3, -2, length=100),
               seq(-2, 0.5, length=100),
               seq(0.5, 6, length=100))
    heatmap.color <- colorRampPalette(c("green", "black", "red"))(n = 1000)
    
    p <- heatmap.2(
      x=t(.data),
      distfun = function(x) {
        x[is.na(x)] <- -1e5
        dist(x, method=distance)
      },
      hclustfun = function(x) {
        x[is.na(x)] <- -1e5
        hclust(x, method=cluster)
      },
      dendrogram="row",
      Rowv=TRUE,
      col=heatmap.color ,
      density.info='none',
      key=TRUE,
      trace="none",
      scale="none",
      #srtCol=45,
      labCol="",
      margins=c(4,9),
      cexRow=1
    )
  }
}
