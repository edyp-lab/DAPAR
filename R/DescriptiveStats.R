##' Boxplot for quantitative proteomics data
##' 
##' @title Builds a boxplot from an object of class \code{\link{MSnSet}}
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param Xaxis A vector containing the indices of columns 
##' in \code{pData()} to use as X-axis (Default is "Label").
##' @return A boxplot
##' @author Florence Combes, Samuel Wieczorek
##' @seealso \code{\link{densityPlotD}}
##' @examples
##' data(UPSprotx2)
##' boxPlotD(UPSprotx2)
boxPlotD <-function(obj, Xaxis="Label"){
    pal<-getPaletteForLabels(obj)
    par(oma = c(2+length(Xaxis), 0, 0, 0))
    .data <- exprs(obj)
    .pData <- pData(obj)
    boxplot(.data
            , las=1
            , col = pal
            , cex = 2
            , axes=TRUE
            , xaxt = "n"
            , ylab = "Log (intensity)"
            , pt.cex = 4
            , horizontal = FALSE
    )

    for (i in Xaxis){
    axis(side=1,
        at = 1:ncol(.data),
        labels = .pData[,i],
        line= 2 * which(i == Xaxis) - 1
    )

    mtext(colnames(.pData)[i],
            1,
            line = 2 * which(i == Xaxis) - 1,
            at=0
            )
    }
    
    mtext("Samples", side=1, line=6+length(Xaxis), cex.lab=1,las=1)

    abline(h=0) 
    palette("default")
}


##' Densityplot of quantitative proteomics data over samples.
##' 
##' @title Builds a densityplot from an object of class \code{\link{MSnSet}}
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param highLightLabel The name of the Label to highlight
##' in the density plot.
##' @param lab2Show A vector of labels to show in densityplot.
##' @return A density plot
##' @author Florence Combes, Samuel Wieczorek
##' @seealso \code{\link{boxPlotD}}, \code{\link{varianceDistD}}
##' @examples data(UPSprotx2)
##' densityPlotD(UPSprotx2)
densityPlotD <- function(obj, highLightLabel=NULL,lab2Show=NULL){
    
    if (is.null(lab2Show)){
        lab2Show <- unique(pData(obj)$Label)
    }
    
    indices <- which(pData(obj)[,"Label"] %in% lab2Show)
    n <- length(indices)
    
    axis.limits <- matrix(data = 0, nrow = 4, ncol = ncol(exprs(obj)))
    
    for (i in 1:ncol(exprs(obj))){
        dens <- density(exprs(obj)[,i], na.rm = TRUE)
        axis.limits[,i]<- c(min(dens$x), max(dens$x), min(dens$y), max(dens$y))
    }
    lim.x <- range(min(axis.limits[1,]), max(axis.limits[2,]))
    lim.y <- range(min(axis.limits[3,]), max(axis.limits[4,]))
    
    ##Colors definition
    nColors <- length(unique(pData(obj)$Label))
    col.density <- c(1:nColors)
    pal <- getPaletteForLabels(obj)
    
    lineWD <- NULL
    lineWD <- c(rep(1, length(colnames(exprs(obj)))))
    if (!is.null(highLightLabel)) {
        lineWD[which(pData(obj)[, "Label"]==highLightLabel)] <- 3
    }
    
    plot(x=NULL
        , ylab="Density"
        , xlab = "log(intensity)"
        ,xlim = lim.x
        ,ylim = lim.y
        ,las = 1
        ,cex.lab = 1
        ,cex.axis = 1
        ,cex.main = 3)
    
    for (i in indices){
        lines(density(exprs(obj)[,i], na.rm=TRUE)
                , col = col.density[which(pData(obj)[i, "Label"]
                                        == unique(pData(obj)$Label))]
                , lwd = lineWD[i])
    }
    
    legend("topright"         
            , legend = unique(pData(obj)$Label)
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
##' of a object. The variance is calculated for each condition (Label) present
##' in the dataset (see the slot \code{'Label'} in the \code{pData()} table)
##' 
##' @title Distribution of variance of proteins
##' @param obj An object of class \code{\link{MSnSet}}.
##' @return A density plot
##' @author Florence Combes, Samuel Wieczorek
##' @seealso \code{\link{densityPlotD}}.
##' @examples data(UPSprotx2)
##' varianceDistD(UPSprotx2)
varianceDistD <-function(obj){
    
    conditions <- unique(pData(obj)[,"Label"])
    n <- length(conditions)
    axis.limits <- matrix(data = 0, nrow = 4, ncol = n)
    for (i in conditions){
        t <- density(apply(exprs(obj)[,which(pData(obj)[,"Label"] == i)], 1, 
                    function(x) var(x, na.rm=TRUE)), na.rm=TRUE)
        axis.limits[,which(conditions == i)]<- c(min(t$x), max(t$x), min(t$y),
                                                max(t$y) )
    }
    
    lim.x <- range(min(axis.limits[1,]), max(axis.limits[2,]))
    lim.y <- range(min(axis.limits[3,]), max(axis.limits[4,]))
    
    par(mar = c(5, 5, 6, 3))
    plot(x=NULL
        , ylab="Density"
        , xlab = "variance(log(intensity))"
        , xlim = lim.x
        , ylim = lim.y
        , las=1
    )
    
    # density by condition
    pal <- getPaletteForLabels(obj)
    conditions <- unique(pData(obj)[,"Label"])
    col.density = c(1:length(conditions))
    for (i in conditions){
        t <- apply(exprs(obj)[,which(pData(obj)[,"Label"] == i)], 1, 
                function(x) var(x, na.rm=TRUE))
        lines(density(t,na.rm=TRUE)
            , xlab=""
            , ylab=""
            , col=col.density[which(conditions == i)]
        )
    }
    
    legend("topright"         
            , legend = conditions
            , col = col.density
            , pch = 15, 
            bty = "n"
            , pt.cex = 2
            , cex = 1
            , horiz = FALSE
            , inset=c(0,0)
    )
}


##' Correlation matrix based on a \code{\link{MSnSet}} object
##' 
##' @title Displays a correlation matrix of the quantitative data of the
##' \code{exprs()} table.
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param indLegend A vector of indices in the columns in
##' the \code{pData()} table choosen for the labels in the axes.
##' @return A colored correlation matrix
##' @author Florence Combes, Samuel Wieczorek
##' @examples data(UPSprotx2)
##' corrMatrixD(UPSprotx2)
corrMatrixD <- function(obj, indLegend=NULL){
    Var1 <- Var2 <- value <- NULL
    if (is.null(indLegend)) {indLegend <- colnames(pData(obj)[-1])} 
    .data <- exprs(obj)
    nb <- length(colnames(.data))
    for (j in 1:nb){
        noms <- NULL
        
        for (i in indLegend){
            #if (pData(obj)[i,j] != "_"){
            noms <- paste(noms, pData(obj)[j,i], sep=" ")
            #}     
        }
        colnames(.data)[j] <- noms
    }
    
    z <- cor(.data,use = 'pairwise.complete.obs')
    text <- element_text(colour="black", size = 16, face = "bold")
    d <- qplot(x=Var1, y=Var2, data= melt(z), fill=value, geom="tile") +
        theme(axis.text = element_text(size=16),
                axis.title = element_text(size=20, face="bold"),
                axis.text.x = element_text(angle=30, vjust=1, hjust=1),
                legend.text = text,
                legend.title = text) +
        labs(x = "", y = "") +
        scale_fill_gradient(low = "lightblue", high = "steelblue")
    plot(d)
}


##' Heatmap of the quantitative proteomic data of a \code{\link{MSnSet}} object
##' 
##' @title This function is a wrapper to \code{\link{heatmap.2}} that displays 
##' quantitative data in the \code{exprs()} table of an object of
##' class \code{\link{MSnSet}}
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param distance The distance used by the clustering algorithm to compute 
##' the dendrogram. See \code{help(heatmap.2)}
##' @param cluster the clustering algorithm used to build the dendrogram.
##' See \code{help(heatmap.2)}
##' @return A heatmap
##' @author Florence Combes, Samuel Wieczorek
##' @examples data(testWithoutNA)
##' heatmapD(testWithoutNA)
heatmapD <- function(obj, distance="euclidean", cluster="average"){
        ##Check parameters
        paramdist<-c("euclidean", "manhattan") 
        if (!(distance %in% paramdist)){
            stop("Param distance is not correct.")
                    return (NULL)
        }
    
    paramcluster<-c("ward.D", "average")
    if (!(cluster %in%  paramcluster)){
        stop("Param clustering is not correct.")
        return (NULL)
    }
    
    if (getNumberOfEmptyLines(obj) != 0)  {
        stop("Your dataset contains empty lines.
            Please filter or impute missing values before.")
        return (NULL)
    }
    else {
    .data <- matrix(exprs(obj), 
                    ncol=ncol(exprs(obj)), 
                    byrow= FALSE,
                    dimnames=list(rownames(exprs(obj)), colnames(exprs(obj)))
    )
    colors = c(seq(-3,-2,length=100),
                seq(-2,0.5,length=100),
                seq(0.5,6,length=100))
    heatmap.color <- colorRampPalette(c("green", "black", "red"))(n = 1000)
    
    p <- heatmap.2(
        x=t(.data),
        distfun=function(x) {
            x[is.na(x)] <- -1e5
            dist(x,method=distance)
        },
        hclustfun=function(x) {
            x[is.na(x)] <- -1e5
            hclust(x,method=cluster)
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
