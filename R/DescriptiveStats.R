##' This function is a wrapper for using the boxPlotD function with objects of 
##' class \code{MSnSet}
##' 
##' @title Wrapper to the boxplotD function on an object \code{MSnSet}
##' @param obj An object of class \code{MSnSet}.
##' @param dataForXAxis A vector of strings containing the names of columns 
##' in \code{pData()} to print labels on X-axis (Default is "Label").
##' @param group2Color A string that indicates how to color the replicates: one
##' color per condition (value "Condition") or one 
##' color per replicate (value "Replicate"). Default value is by Condition.
##' @return A boxplot
##' @author Florence Combes, Samuel Wieczorek
##' @seealso \code{\link{wrapper.densityPlotD}}
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' types <- c("Label","Analyt.Rep")
##' wrapper.boxPlotD(Exp1_R25_pept, types)
wrapper.boxPlotD <- function(obj, 
                            dataForXAxis="Label", 
                            group2Color="Condition"){

qData <- Biobase::exprs(obj)
#dataForXAxis <- as.matrix(Biobase::pData(obj)[,dataForXAxis])
dataForXAxis <- as.matrix(Biobase::pData(obj)[,dataForXAxis])

labels <- Biobase::pData(obj)[,"Label"]

boxPlotD(qData, dataForXAxis, labels, group2Color)

}







##' This function is a wrapper for using the boxPlotD_HC function with objects of 
##' class \code{MSnSet}
##' 
##' @title Wrapper to the boxplotD_HC function on an object \code{MSnSet}
##' @param obj An object of class \code{MSnSet}.
##' @param dataForXAxis A vector of strings containing the names of columns 
##' in \code{pData()} to print labels on X-axis (Default is "Label").
##' @param group2Color A string that indicates how to color the replicates: one
##' color per condition (value "Condition") or one 
##' color per replicate (value "Replicate"). Default value is by Condition.
##' @return A boxplot
##' @author Samuel Wieczorek
##' @seealso \code{\link{wrapper.densityPlotD}}
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' types <- c("Label","Analyt.Rep")
##' wrapper.boxPlotD_HC(Exp1_R25_pept, types)
wrapper.boxPlotD_HC <- function(obj, 
                             dataForXAxis="Label", 
                             group2Color="Condition"){
    
    qData <- Biobase::exprs(obj)
    #dataForXAxis <- as.matrix(Biobase::pData(obj)[,dataForXAxis])
    dataForXAxis <- as.matrix(Biobase::pData(obj)[,dataForXAxis])
    
    labels <- Biobase::pData(obj)[,"Label"]
    
    boxPlotD_HC(qData, dataForXAxis, labels, group2Color)
    
}








##' Boxplot for quantitative proteomics data
##' 
##' @title Builds a boxplot from a dataframe
##' @param qData A dataframe that contains quantitative data.
##' @param dataForXAxis A vector containing the types of replicates 
##' to use as X-axis. Available values are: Label, Analyt.Rep,
##' Bio.Rep and Tech.Rep. Default is "Label".
##' @param labels A vector of the conditions (labels) (one label per sample).
##' @param group2Color A string that indicates how to color the replicates: 
##' one color per condition (value "Condition") or one 
##' color per replicate (value "Replicate"). Default value is by Condition.
##' @return A boxplot
##' @author Florence Combes, Samuel Wieczorek
##' @seealso \code{\link{densityPlotD}}
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' types <- c("Label","Analyt.Rep")
##' dataForXAxis <- Biobase::pData(Exp1_R25_pept)[,types]
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
##' boxPlotD(qData, dataForXAxis, labels)
boxPlotD <- function(qData, 
                    dataForXAxis=NULL, 
                    labels=NULL, 
                    group2Color="Condition"){

if (group2Color == "Condition") {
    pal <- getPaletteForLabels(labels)
}else { 
        pal <- getPaletteForReplicates(ncol(qData))}

    if (is.null(labels)){size <- 2}else{size <- 2*length(labels)}
   #par(mar=c(size,4,2,2))
        
boxplot(qData
        ,las = 1
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
if (is.vector(dataForXAxis) ){
    xAxisLegend <- dataForXAxis
    } 
    else{ 
        xAxisLegend <- NULL
        N <- nrow(dataForXAxis)
        for (i in 1:N){
            xAxisLegend <- c(xAxisLegend, paste(dataForXAxis[i,], collapse="_"))
        }
    }
    axis(side=1,
         at = 1:ncol(qData),
         labels = xAxisLegend)

mtext("Samples", 
    side=1, 
    line=6+length(colnames(dataForXAxis)), 
    cex.lab=1, las=1)
}

abline(h=0) 
palette("default")
}




##' Boxplot for quantitative proteomics data using the library \code{highcharter}
##' 
##' @title Builds a boxplot from a dataframe using the library \code{highcharter}
##' @param qData A dataframe that contains quantitative data.
##' @param dataForXAxis A vector containing the types of replicates 
##' to use as X-axis. Available values are: Label, Analyt.Rep,
##' Bio.Rep and Tech.Rep. Default is "Label".
##' @param labels A vector of the conditions (labels) (one label per sample).
##' @param group2Color A string that indicates how to color the replicates: 
##' one color per condition (value "Condition") or one 
##' color per replicate (value "Replicate"). Default value is by Condition.
##' @return A boxplot
##' @author Samuel Wieczorek
##' @seealso \code{\link{densityPlotD_HC}}
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' types <- c("Label","Analyt.Rep")
##' dataForXAxis <- Biobase::pData(Exp1_R25_pept)[,types]
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
##' boxPlotD_HC(qData, dataForXAxis, labels)
boxPlotD_HC <- function(qData, 
                     dataForXAxis="Label", 
                     labels=NULL, 
                     group2Color="Condition"){

    
    if (is.null(labels) ) {
        labels <- rep("",ncol(qData))
    }

    
bx <-boxplot(qData, na.rm=TRUE)

df_outlier <- data.frame(x=bx$group-1,y = bx$out)

titi <- NULL
nc <- ncol(qData)
for (i in 1:nc){
    titi <- c(titi, rep(colnames(qData)[i],nrow(qData)))
}

df <- data.frame(values = as.vector(qData,mode='numeric'),
                 samples = titi)

pal <- c("#002F80", "#002F80","#002F80","#002F80","#F9AF38","#F9AF38","#F9AF38","#F9AF38")

myColors <- NULL
##Colors definition
if (group2Color == "Condition") {
    myColors <- getPaletteForLabels_HC(labels)
}else { 
    myColors <- getPaletteForReplicates_HC(ncol(qData))}



hcboxplot(x=df$values, var = df$samples, 
          colorByPoint = TRUE, 
          outliers = FALSE) %>%
    
    hc_chart(type="column") %>%
    hc_yAxis(title = list(text = "Log (intensity)")) %>%
    hc_xAxis(title = list(text = "Samples")) %>%
    hc_colors(myColors) %>%
    hc_add_series(type= "scatter",df_outlier) %>%
    hc_tooltip(enabled = FALSE) %>%
    hc_plotOptions(
        boxplot= list(
            fillColor= myColors,
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





}



##' This function is a wrapper for using the violinPlotD function with objects 
##' of class \code{MSnSet}
##' 
##' @title Wrapper to the violinPlotD function on an object \code{MSnSet}
##' @param obj An object of class \code{MSnSet}.
##' @param dataForXAxis A vector of strings containing the names of columns 
##' in \code{pData()} to print labels on X-axis (Default is "Label").
##' @param group2Color A string that indicates how to color the replicates: one
##' color per condition (value "Condition") or one 
##' color per replicate (value "Replicate"). Default value is by Condition.
##' @return A violin plot
##' @author Samuel Wieczorek
##' @seealso \code{\link{wrapper.densityPlotD}}, \code{\link{wrapper.boxPlotD}}
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' library(vioplot)
##' types <- c("Label","Analyt.Rep")
##' wrapper.violinPlotD(Exp1_R25_pept, types)
wrapper.violinPlotD <- function(obj, 
                             dataForXAxis="Label", 
                             group2Color="Condition"){
    
    qData <- Biobase::exprs(obj)
    dataForXAxis <- as.matrix(Biobase::pData(obj)[,dataForXAxis])
    labels <- Biobase::pData(obj)[,"Label"]
    
    violinPlotD(qData, dataForXAxis, labels, group2Color)
    
}




##' ViolinPlot for quantitative proteomics data
##' 
##' @title Builds a violinplot from a dataframe
##' @param qData A dataframe that contains quantitative data.
##' @param dataForXAxis A vector containing the types of replicates 
##' to use as X-axis. Available values are: Label, Analyt.Rep,
##' Bio.Rep and Tech.Rep. Default is "Label".
##' @param labels A vector of the conditions (labels) (one label per sample).
##' @param group2Color A string that indicates how to color the replicates: 
##' one color per condition (value "Condition") or one 
##' color per replicate (value "Replicate"). Default value is by Condition.
##' @return A violinplot
##' @author Florence Combes, Samuel Wieczorek
##' @seealso \code{\link{densityPlotD}}
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' library(vioplot)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' types <- c("Label","Analyt.Rep")
##' dataForXAxis <- Biobase::pData(Exp1_R25_pept)[,types]
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
##' violinPlotD(qData, dataForXAxis, labels)
violinPlotD <- function(qData, 
                     dataForXAxis=NULL, 
                     labels=NULL, 
                     group2Color="Condition"){

    plot.new()
    if (group2Color == "Condition") {
        pal <- getPaletteForLabels(labels)
    }else { pal <- getPaletteForReplicates(ncol(qData))}
    
    plot.window(xlim=c(0,ncol(qData)+1),
                ylim=c(min(na.omit(qData)),max(na.omit(qData))))
    title( ylab="Log (intensity)")
    for (i in 1:ncol(qData)) {
        vioplot(na.omit(qData[,i]), col = pal[i], add=TRUE, at=i)}
    
    
    axis(2, yaxp = c(floor(min(na.omit(qData))), 
                     floor(max(na.omit(qData))), 5), las=1)
    
     if( !is.null(dataForXAxis))
     {
        if (is.vector(dataForXAxis) ){
            N <- 1} else{ N <- ncol(dataForXAxis)}
         
        for (i in 1:N){
            axis(side=1,
                 at = 1:ncol(qData),
                 labels = if (is.vector(dataForXAxis) ) 
                     {dataForXAxis} else {dataForXAxis[,i]},
                 line= 2*i-1
            )
        }

        mtext("Samples",
              side=1,
              line=6+length(colnames(dataForXAxis)),
              cex.lab=1, las=1)
    }


     palette("default")
}



##' Wrapper to the function that plot to compare the quantitative proteomics 
##' data before and after normalization
##' 
##' @title Builds a plot from a dataframe
##' @param objBefore A dataframe that contains quantitative data before 
##' normalization.
##' @param objAfter A dataframe that contains quantitative data after 
##' normalization.
##' @param labelsForLegend A vector of the conditions (labels) (one label 
##' per sample).
##' @param indData2Show A vector of the indices of the columns to show in the 
##' plot. The indices are those of indices of 
##' the columns int the data.frame qDataBefore.
##' @param group2Color A string that indicates how to color the replicates: 
##' one color per condition (value "Condition") or one 
##' color per replicate (value "Replicate"). Default value is by Condition.
##' @return A plot
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
##' objAfter <- wrapper.normalizeD(Exp1_R25_pept, "Median Centering", 
##' "within conditions")
##' wrapper.compareNormalizationD(Exp1_R25_pept, objAfter, labels)
wrapper.compareNormalizationD <- function(objBefore, objAfter, 
                                        labelsForLegend=NULL,
                                        indData2Show=NULL,
                                        group2Color="Condition"){

qDataBefore <- Biobase::exprs(objBefore)
qDataAfter <- Biobase::exprs(objAfter)
if (is.null(labelsForLegend)){
  labelsForLegend <- Biobase::pData(objBefore)[,"Label"]}

compareNormalizationD(qDataBefore, qDataAfter, labelsForLegend, indData2Show, 
                    group2Color)
}

##' Wrapper to the function that plot to compare the quantitative proteomics 
##' data before and after normalization. Same as the function \link{wrapper.compareNormalizationD}
##' but uses the package \code{highcharter}
##' 
##' @title Builds a plot from a dataframe
##' @param objBefore A dataframe that contains quantitative data before 
##' normalization.
##' @param objAfter A dataframe that contains quantitative data after 
##' normalization.
##' @param labelsForLegend A vector of the conditions (labels) (one label 
##' per sample).
##' @param indData2Show A vector of the indices of the columns to show in the 
##' plot. The indices are those of indices of 
##' the columns int the data.frame qDataBefore.
##' @param group2Color A string that indicates how to color the replicates: 
##' one color per condition (value "Condition") or one 
##' color per replicate (value "Replicate"). Default value is by Condition.
##' @return A plot
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
##' objAfter <- wrapper.normalizeD(Exp1_R25_pept, "Median Centering", 
##' "within conditions")
##' wrapper.compareNormalizationD_HC(Exp1_R25_pept, objAfter, labels)
wrapper.compareNormalizationD_HC <- function(objBefore, objAfter, 
                                          labelsForLegend=NULL,
                                          indData2Show=NULL,
                                          group2Color="Condition"){
    
    qDataBefore <- Biobase::exprs(objBefore)
    qDataAfter <- Biobase::exprs(objAfter)
    
    compareNormalizationD_HC(qDataBefore, qDataAfter, labelsForLegend, indData2Show, 
                          group2Color)
}

##' Plot to compare the quantitative proteomics data before and after 
##' normalization
##' 
##' @title Builds a plot from a dataframe
##' @param qDataBefore A dataframe that contains quantitative data before 
##' normalization.
##' @param qDataAfter A dataframe that contains quantitative data after 
##' normalization.
##' @param labelsForLegend A vector of the conditions (labels) (one label 
##' per sample).
##' @param indData2Show A vector of the indices of the columns to show in 
##' the plot. The indices are those of indices of 
##' the columns int the data.frame qDataBefore.
##' @param group2Color A string that indicates how to color the replicates:
##' one color per condition (value "Condition") or one 
##' color per replicate (value "Replicate"). Default value is by Condition.
##' @return A plot
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qDataBefore <- Biobase::exprs(Exp1_R25_pept)
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
##' qDataAfter <- normalizeD(qDataBefore,labels,"Median Centering",
##' "within conditions")
##' compareNormalizationD(qDataBefore, qDataAfter, labels)
compareNormalizationD <- function(qDataBefore,
                                qDataAfter,
                                labelsForLegend=NULL,
                                indData2Show=NULL,
                                group2Color="Condition"){

if (is.null(labelsForLegend)) return(NULL)
if (is.null(indData2Show)) {indData2Show <- c(1:ncol(qDataAfter)) }

x <- qDataBefore
y <- qDataAfter/qDataBefore

lim.x <- range(min(x, na.rm=TRUE), max(x, na.rm=TRUE))
lim.y <- range(min(y, na.rm=TRUE), max(y, na.rm=TRUE))


##Colors definition
if (group2Color == "Condition") {
    pal <- getPaletteForLabels(labelsForLegend)
    legendColor <- unique(pal)
    txtLegend <- unique(labelsForLegend)
}else { 
    pal <- getPaletteForReplicates(ncol(x))
    legendColor <- pal[indData2Show]
    txtLegend <- paste("Replicate", seq(1,ncol(x)), labelsForLegend,sep=" ")
    txtLegend <- txtLegend[indData2Show]
}


plot(x=NULL
    ,xlim = lim.x
    ,ylim = lim.y
    , cex = 1
    , axes=TRUE
    , xlab = "Intensities before normalization"
    , ylab = "Intensities after normalization / Intensities before 
    normalization"
    ,cex.lab = 1
    ,cex.axis = 1
    ,cex.main = 3)


for (i in indData2Show){
    points(x[,i], y[,i], col = pal[i], cex = 1,pch=16)
}

legend("topleft"
        , legend = txtLegend
        , col = legendColor
        , pch = 15 
        , bty = "n"
        , pt.cex = 2
        , cex = 1
        , horiz = FALSE
        , inset=c(0,0)
)


palette("default")
}



##' Plot to compare the quantitative proteomics data before and after 
##' normalization using the library \code{highcharter}
##' 
##' @title Builds a plot from a dataframe. Same as compareNormalizationD but 
##' uses the library \code{highcharter}
##' @param qDataBefore A dataframe that contains quantitative data before 
##' normalization.
##' @param qDataAfter A dataframe that contains quantitative data after 
##' normalization.
##' @param labelsForLegend A vector of the conditions (labels) (one label 
##' per sample).
##' @param indData2Show A vector of the indices of the columns to show in 
##' the plot. The indices are those of indices of 
##' the columns int the data.frame qDataBefore.
##' @param group2Color A string that indicates how to color the replicates:
##' one color per condition (value "Condition") or one 
##' color per replicate (value "Replicate"). Default value is by Condition.
##' @return A plot
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qDataBefore <- Biobase::exprs(Exp1_R25_pept)
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
##' qDataAfter <- normalizeD(qDataBefore,labels,"Median Centering",
##' "within conditions")
##' compareNormalizationD_HC(qDataBefore, qDataAfter, labels)
compareNormalizationD_HC <- function(qDataBefore,
                                  qDataAfter,
                                  labelsForLegend=NULL,
                                  indData2Show=NULL,
                                  group2Color="Condition"){
    
    if (is.null(labelsForLegend)) return(NULL)
    if (is.null(indData2Show)) {indData2Show <- c(1:ncol(qDataAfter)) }
    
    x <- qDataBefore
    y <- qDataAfter/qDataBefore
   
    ##Colors definition
    if (group2Color == "Condition") {
        pal <- getPaletteForLabels(labelsForLegend)
        legendColor <- unique(pal)
        txtLegend <- unique(labelsForLegend)
    }else { 
        pal <- getPaletteForReplicates(ncol(x))
        legendColor <- pal[indData2Show]
        txtLegend <- paste("Replicate", seq(1,ncol(x)), labelsForLegend,sep=" ")
        txtLegend <- txtLegend[indData2Show]
    }
    
    
    series <- list()
    for (i in 1:length(indData2Show)){
        tmp <- list(name=labelsForLegend[i], data =list_parse(data.frame(x=x[,indData2Show[i]],y=y[,indData2Show[i]])))
        series[[i]] <- tmp
    }
   
    h1 <-  highchart() %>% 
        my_hc_chart( chartType = "scatter") %>%
        hc_add_series_list(series) %>%
        hc_tooltip(enabled= "false" ) %>%
        my_hc_ExportMenu(filename = "compareNormalization")
    h1
    return(h1)
    
    
    
    # 
    # plot(x=NULL
    #      ,xlim = lim.x
    #      ,ylim = lim.y
    #      , cex = 1
    #      , axes=TRUE
    #      , xlab = "Intensities before normalization"
    #      , ylab = "Intensities after normalization / Intensities before 
    #      normalization"
    #      ,cex.lab = 1
    #      ,cex.axis = 1
    #      ,cex.main = 3)
    # 
    # 
    # for (i in indData2Show){
    #     points(x[,i], y[,i], col = pal[i], cex = 1,pch=16)
    # }
    # 
    # legend("topleft"
    #        , legend = txtLegend
    #        , col = legendColor
    #        , pch = 15 
    #        , bty = "n"
    #        , pt.cex = 2
    #        , cex = 1
    #        , horiz = FALSE
    #        , inset=c(0,0)
    # )
    # 
    # 
    # palette("default")
}





##' This function is a wrapper for using the densityPlotD function with 
##' objects of class \code{MSnSet}
##' 
##' @title Builds a densityplot from an object of class \code{MSnSet}
##' @param obj An object of class \code{MSnSet}.
##' @param labelsForLegend A vector of labels to show in densityplot.
##' @param indData2Show A vector of the indices of the columns to show in 
##' the plot. The indices are those of indices of the columns int the data
##' frame qDataBefore in the density plot.
##' @param group2Color A string that indicates how to color the replicates: 
##' one color per condition (value "Condition") or one color per replicate
##' (value "Replicate"). Default value is by Condition.
##' @return A density plot
##' @author Alexia Dorffer
##' @seealso \code{\link{wrapper.boxPlotD}}, 
##' \code{\link{wrapper.CVDistD}}
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
##' wrapper.densityPlotD(Exp1_R25_pept, labels)
wrapper.densityPlotD <- function(obj, labelsForLegend=NULL,  indData2Show=NULL,
                                group2Color = "Condition"){
qData <- Biobase::exprs(obj)
if (is.null(labelsForLegend) ) {
  labelsForLegend <- Biobase::pData(obj)[,"Label"]}

densityPlotD(qData, labelsForLegend, indData2Show,group2Color)
}

##' This function is a wrapper for using the densityPlotD function with 
##' objects of class \code{MSnSet}. Same as the function \code{\link{wrapper.densityPlotD}}
##' but uses the package \code{highcharter}
##' 
##' @title Builds a densityplot from an object of class \code{MSnSet}
##' @param obj An object of class \code{MSnSet}.
##' @param labelsForLegend A vector of labels to show in densityplot.
##' @param indData2Show A vector of the indices of the columns to show in 
##' the plot. The indices are those of indices of the columns int the data
##' frame qDataBefore in the density plot.
##' @param group2Color A string that indicates how to color the replicates: 
##' one color per condition (value "Condition") or one color per replicate
##' (value "Replicate"). Default value is by Condition.
##' @return A density plot
##' @author Samuel Wieczorek
##' @seealso \code{\link{wrapper.boxPlotD}}, 
##' \code{\link{wrapper.CVDistD}}
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
##' wrapper.densityPlotD_HC(Exp1_R25_pept, labels)
wrapper.densityPlotD_HC <- function(obj, labelsForLegend=NULL,  indData2Show=NULL,
                                 group2Color = "Condition"){
    qData <- Biobase::exprs(obj)
    
    if (is.null(labelsForLegend) ) {
      labelsForLegend <- Biobase::pData(obj)[,"Label"]}
    
    densityPlotD_HC(qData, labelsForLegend, indData2Show,group2Color)
}



##' Densityplot of quantitative proteomics data over samples.
##' 
##' @title Builds a densityplot from a dataframe
##' @param qData A dataframe that contains quantitative data.
##' @param labelsForLegend A vector of the conditions (labels) (one label 
##' per sample).
##' @param indData2Show A vector of indices to show in densityplot. If NULL, 
##' then all labels are displayed.
##' @param group2Color A string that indicates how to color the replicates: 
##' one color per condition (value "Condition") or one 
##' color per replicate (value "Replicate"). Default value is by Condition.
##' @return A density plot
##' @author Florence Combes, Samuel Wieczorek
##' @seealso \code{\link{boxPlotD}}, \code{\link{CVDistD}}
##' @examples 
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' labels <- lab2Show <- Biobase::pData(Exp1_R25_pept)[,"Label"]
##' densityPlotD(qData, labels)
densityPlotD <- function(qData, labelsForLegend=NULL,indData2Show=NULL,
                        group2Color = "Condition"){
    
#if (is.null(labelsForLegend) ) {
#  labelsForLegend <- Biobase::pData(rv$current.obj)[,"Label"]}

if (is.null(indData2Show)) {indData2Show <- c(1:ncol(qData)) }

### Range of axis definition
axis.limits <- matrix(data = 0, nrow = 4, ncol = ncol(qData))
for (i in 1:ncol(qData)){
    dens <- density(qData[,i], na.rm = TRUE)
    axis.limits[,i] <- c(min(dens$x), max(dens$x), min(dens$y), max(dens$y))
    }
lim.x <- range(min(axis.limits[1,]), max(axis.limits[2,]))
lim.y <- range(min(axis.limits[3,]), max(axis.limits[4,]))


##Colors definition
if (group2Color == "Condition") {
    pal <- getPaletteForLabels(labelsForLegend)
    legendColor <- unique(pal)
    txtLegend <- unique(labelsForLegend)
}else { 
    pal <- getPaletteForReplicates(ncol(qData))
    legendColor <- pal[indData2Show]
    txtLegend <- paste("Replicate", seq(1,ncol(qData)), 
                        labelsForLegend,sep=" ")
    txtLegend <- txtLegend[indData2Show]
}

###Erase data not to show (color in white)
# lineWD <- NULL
# lineWD <- c(rep(1, length(colnames(qData))))
# if (!is.null(highLightLabel)) {
#   lineWD[which(labels == highLightLabel)] <- 3
# }
# 
plot(x =NULL
    , ylab ="Density"
    , xlab = "log(intensity)"
    , col = pal
    ,xlim = lim.x
    ,ylim = lim.y
    ,las = 1
    ,cex.lab = 1
    ,cex.axis = 1
    ,cex.main = 3)

for (i in indData2Show){
    lines(density(qData[,i], na.rm=TRUE), col = pal[i])
}


legend("topleft"         
        , legend = txtLegend
        , col = legendColor
        , pch = 15 
        , bty = "n"
        , pt.cex = 2
        , cex = 1
        , horiz = FALSE
        , inset=c(0,0)
)
}




##' Densityplot of quantitative proteomics data over samples. Same as the function \code{\link{densityPlotD}}
##' but uses the package \code{highcharter}
##' 
##' @title Builds a densityplot from a dataframe
##' @param qData A dataframe that contains quantitative data.
##' @param labelsForLegend A vector of the conditions (labels) (one label 
##' per sample).
##' @param indData2Show A vector of indices to show in densityplot. If NULL, 
##' then all labels are displayed.
##' @param group2Color A string that indicates how to color the replicates: 
##' one color per condition (value "Condition") or one 
##' color per replicate (value "Replicate"). Default value is by Condition.
##' @return A density plot
##' @author Samuel Wieczorek
##' @seealso \code{\link{boxPlotD}}, \code{\link{CVDistD}}
##' @examples 
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' labels <- lab2Show <- Biobase::pData(Exp1_R25_pept)[,"Label"]
##' densityPlotD_HC(qData, labels)
densityPlotD_HC <- function(qData, labelsForLegend=NULL,indData2Show=NULL,
                            group2Color = "Condition"){
    
    
    if (is.null(labelsForLegend) ) {
      labelsForLegend <- rep("",ncol(qData))
      }
    
    if (is.null(indData2Show)) {
      indData2Show <- c(1:ncol(qData)) 
      }
    
    pal <- c("#002F80", "#002F80","#002F80","#002F80","#F9AF38","#F9AF38","#F9AF38","#F9AF38")
    
    myColors <- NULL
    ##Colors definition
    if (group2Color == "Condition") {
      myColors <- getPaletteForLabels_HC(labelsForLegend)[indData2Show]
    } else { 
      myColors<- getPaletteForReplicates_HC(ncol(qData))[indData2Show]
    }

    series <- list()
    for (i in 1:length(indData2Show)){
  
    tmp <- data.frame(x = density(qData[,indData2Show[i]], na.rm = TRUE)$x, 
                y = density(qData[,indData2Show[i]], na.rm = TRUE)$y)
    series[[i]] <- list(name = labelsForLegend[i],
                              data = list_parse(tmp))
}


    h1 <-  highchart() %>% 
        hc_title(text = "Density plot") %>% 
        my_hc_chart(chartType = "spline", zoomType="x") %>%
        hc_add_series_list(series) %>%
        hc_legend(enabled = TRUE) %>%
        hc_xAxis(title = list(text = "log(Intensity)")) %>%
        hc_yAxis(title = list(text = "Density")) %>%
        hc_colors(myColors) %>%
        hc_tooltip(headerFormat= '',
                   pointFormat = "<b> {series.name} </b>: {point.y} ",
                   valueDecimals = 2) %>%
      my_hc_ExportMenu(filename = "densityplot") %>%
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
    
    
    return(h1)

}




##' Builds a densityplot of the CV of entities in the exprs() table
##' of an object \code{MSnSet}. The variance is calculated for each 
##' condition (Label) present
##' in the dataset (see the slot \code{'Label'} in the \code{pData()} table).
##' 
##' @title Distribution of CV of entities
##' @param obj An object of class \code{MSnSet}.
##' @return A density plot
##' @author Alexia Dorffer
##' @seealso \code{\link{wrapper.densityPlotD}}
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.CVDistD(Exp1_R25_pept)
wrapper.CVDistD <- function(obj){
qData <- Biobase::exprs(obj)
labels <- Biobase::pData(obj)[,"Label"]
CVDistD(qData, labels)
}


##' Builds a densityplot of the CV of entities in the exprs() table. 
##' of an object \code{MSnSet}. The variance is calculated for each 
##' condition (Label) present
##' in the dataset (see the slot \code{'Label'} in the \code{pData()} table).
##' Same as the function \code{\link{wrapper.CVDistD}} but uses the package \code{highcharter}
##' 
##' @title Distribution of CV of entities
##' @param obj An object of class \code{MSnSet}.
##' @return A density plot
##' @author Samuel Wieczorek
##' @seealso \code{\link{wrapper.densityPlotD}}
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.CVDistD_HC(Exp1_R25_pept)
wrapper.CVDistD_HC <- function(obj){
    qData <- Biobase::exprs(obj)
    labels <- Biobase::pData(obj)[,"Label"]
    CVDistD_HC(qData, labels)
}


##' Builds a densityplot of the CV of entities in the exprs() table
##' of a object. The CV is calculated for each condition (Label) present
##' in the dataset (see the slot \code{'Label'} in the \code{pData()} table)
##' 
##' @title Distribution of CV of entities
##' @param qData A dataframe that contains quantitative data.
##' @param labels A vector of the conditions (labels) (one label per sample).
##' @return A density plot
##' @author Florence Combes, Samuel Wieczorek
##' @seealso \code{\link{densityPlotD}}.
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
##' CVDistD(Biobase::exprs(Exp1_R25_pept), labels)
CVDistD <- function(qData, labels=NULL){
    
if (is.null(labels)) {return(NULL)}
conditions <- unique(labels)
n <- length(conditions)
axis.limits <- matrix(data = 0, nrow = 4, ncol = n)
for (i in conditions){
    if (length(which(labels == i)) > 1){
    t <- density(apply(qData[,which(labels == i)], 1, 
                    function(x) 100*var(x, na.rm=TRUE)/mean(x, na.rm=TRUE)), 
                 na.rm=TRUE)

    axis.limits[,which(conditions == i)]<- c(min(t$x), max(t$x), min(t$y),
                                            max(t$y))
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
pal <- getPaletteForLabels(labels)
conditions <- unique(labels)
col.density = c(1:length(conditions))
for (i in conditions){
    if (length(which(labels == i)) > 1){
        t <- apply(qData[,which(labels == i)], 1, 
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



##' Builds a densityplot of the CV of entities in the exprs() table
##' of a object. The CV is calculated for each condition (Label) present
##' in the dataset (see the slot \code{'Label'} in the \code{pData()} table)
##' Same as the function \code{CVDistD} but uses the package \code{highcharter}
##' @title Distribution of CV of entities
##' @param qData A dataframe that contains quantitative data.
##' @param labels A vector of the conditions (labels) (one label per sample).
##' @return A density plot
##' @author Samuel Wieczorek
##' @seealso \code{\link{densityPlotD}}.
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
##' CVDistD_HC(Biobase::exprs(Exp1_R25_pept), labels)
CVDistD_HC <- function(qData, labels=NULL){
    
    if (is.null(labels)) {return(NULL)}
    conditions <- unique(labels)
    n <- length(conditions)
    
    nbSeries = n
    series <- list()
    for (i in 1:length(conditions)){
        if (length(which(labels == conditions[i])) > 1){
            t <- apply(qData[,which(labels == conditions[i])], 1, 
                       function(x) 100*var(x, na.rm=TRUE)/mean(x, na.rm=TRUE))
            tmp <- data.frame(x = density(t, na.rm = TRUE)$x,
                              y = density(t, na.rm = TRUE)$y)
            series[[i]] <- list(name = conditions[i],
                                data = list_parse(tmp))
        }
    }

    h1 <-  highchart() %>% 
        my_hc_chart(chartType = "spline", zoomType="x") %>%
        hc_add_series_list(series) %>%
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
    
    return(h1)

}






##' Builds a correlation matrix based on a \code{MSnSet} object.
##' 
##' @title Displays a correlation matrix of the quantitative data of the
##' \code{exprs()} table
##' @param obj An object of class \code{MSnSet}.
##' @param rate A float that defines the gradient of colors.
##' @return A colored correlation matrix
##' @author Alexia Dorffer
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.corrMatrixD(Exp1_R25_pept)
wrapper.corrMatrixD <- function(obj, rate=5){
qData <- Biobase::exprs(obj)
samplesData <- Biobase::pData(obj)
corrMatrixD(qData, samplesData, rate)
}

##' Builds a correlation matrix based on a \code{MSnSet} object. 
##' Same as the function \code{\link{wrapper.corrMatrixD}} but uses the package \code{highcharter}
##' 
##' @title Displays a correlation matrix of the quantitative data of the
##' \code{exprs()} table
##' @param obj An object of class \code{MSnSet}.
##' @param rate A float that defines the gradient of colors.
##' @return A colored correlation matrix
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.corrMatrixD_HC(Exp1_R25_pept)
wrapper.corrMatrixD_HC <- function(obj, rate=0.5){
    qData <- Biobase::exprs(obj)
    samplesData <- Biobase::pData(obj)
    data <- cor(qData,use = 'pairwise.complete.obs')
    corrMatrixD_HC(data,samplesData, rate)
}



##' Correlation matrix based on a \code{MSnSet} object
##' 
##' @title Displays a correlation matrix of the quantitative data of the
##' \code{exprs()} table.
##' @param qData A dataframe of quantitative data.
##' @param samplesData A dataframe where lines correspond to samples and 
##' columns to the meta-data for those samples.
##' @param gradientRate The rate parameter to control the exponential law for 
##' the gradient of colors
##' @return A colored correlation matrix
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' samplesData <- Biobase::pData(Exp1_R25_pept)
##' corrMatrixD(qData, samplesData)
corrMatrixD <- function(qData, samplesData, gradientRate = 5){
Var1 <- Var2 <- value <- NULL

for (j in 1:length(colnames(qData))){
    colnames(qData)[j] <- paste(as.character(samplesData[j,2:5]), 
                                collapse =" ")
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

    scale_fill_gradientn (
    colours=colorRampPalette (c ("white", "lightblue","darkblue")) (101),
    values = c(pexp(seq(0,1,0.01), rate=gradientRate),1), limits=c(0,1))

plot(d)
}




##' Correlation matrix based on a \code{MSnSet} object. Same as the 
##' function \link{corrMatrixD} but uses the package \code{highcharter}
##' 
##' @title Displays a correlation matrix of the quantitative data of the
##' \code{exprs()} table.
##' @param object The result of the \code{cor} function.
##' @param samplesData A dataframe in which lines correspond to samples and 
##' columns to the meta-data for those samples.
##' @param rate The rate parameter to control the exponential law for 
##' the gradient of colors
##' @return A colored correlation matrix
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' samplesData <- Biobase::pData(Exp1_R25_pept)
##' res <- cor(qData,use = 'pairwise.complete.obs')
##' corrMatrixD_HC(res, samplesData)
corrMatrixD_HC <- function(object,samplesData = NULL, rate = 0.5) {
    
    df <- as.data.frame(object)
    
    if (!is.null(samplesData)){
        for (j in 1:ncol(df)){
            names(df)[j] <- paste(as.character(samplesData[j,2:5]), 
                                        collapse =" ")
        }
        }
    is.num <- sapply(df, is.numeric)
    df[is.num] <- lapply(df[is.num], round, 2)
    dist <- NULL
    
    x <- y <- names(df)
    
    df <- tbl_df(cbind(x = y, df)) %>% 
        gather(y, dist, -x) %>% 
        mutate(x = as.character(x),
               y = as.character(y)) %>% 
        left_join(data_frame(x = y,
                             xid = seq(length(y)) - 1), by = "x") %>% 
        left_join(data_frame(y = y,
                             yid = seq(length(y)) - 1), by = "y")
    
    ds <- df %>% 
        select_("xid", "yid", "dist") %>% 
        list_parse2()
    
    fntltp <- JS("function(){
                  return this.series.xAxis.categories[this.point.x] + ' ~ ' +
                         this.series.yAxis.categories[this.point.y] + ': <b>' +
                         Highcharts.numberFormat(this.point.value, 2)+'</b>';
               ; }")
    cor_colr <- list( list(0, '#FF5733'),
                      list(0.5, '#F8F5F5'),
                      list(1, '#2E86C1')
    )
    highchart() %>% 
        my_hc_chart(chartType = "heatmap") %>% 
        hc_xAxis(categories = y, title = NULL) %>% 
        hc_yAxis(categories = y, title = NULL) %>% 
        hc_add_series(data = ds) %>% 
        hc_plotOptions(
            series = list(
                boderWidth = 0,
                dataLabels = list(enabled = TRUE)
            )) %>% 
        hc_tooltip(formatter = fntltp) %>% 
        hc_legend(align = "right", layout = "vertical",
                  verticalAlign="middle") %>% 
        hc_colorAxis(  stops= cor_colr,min=rate,max=1) %>%
      my_hc_ExportMenu(filename = "corrMatrix")
}




##' Builds a heatmap of the quantitative proteomic data of a 
##' \code{MSnSet} object.
##' 
##' @title This function is a wrapper to \code{\link{heatmap.2}} that displays 
##' quantitative data in the \code{exprs()} table of an object of
##' class \code{MSnSet}
##' @param obj An object of class \code{MSnSet}.
##' @param distance The distance used by the clustering algorithm to compute 
##' the dendrogram. See \code{help(heatmap.2)}.
##' @param cluster the clustering algorithm used to build the dendrogram.
##' See \code{help(heatmap.2)}
##' @param dendro A boolean to indicate fi the dendrogram has to be displayed
##' @return A heatmap
##' @author Alexia Dorffer
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- mvFilter(Exp1_R25_pept[1:1000], "wholeMatrix", 6)
##' wrapper.heatmapD(obj)
wrapper.heatmapD  <- function(obj, distance="euclidean", cluster="complete", 
                            dendro = FALSE){
qData <- Biobase::exprs(obj)
for (j in 1:length(colnames(qData))){
    colnames(qData)[j] <- paste(as.character(Biobase::pData(obj)[j,2:5]), 
                                collapse =" ")
}

heatmapD(qData, distance, cluster, dendro)
}




##' Heatmap of the quantitative proteomic data of a \code{MSnSet} object
##' 
##' @title This function is a wrapper to \code{\link{heatmap.2}} that displays 
##' quantitative data in the \code{exprs()} table of an object of
##' class \code{MSnSet}
##' @param qData A dataframe that contains quantitative data.
##' @param distance The distance used by the clustering algorithm to compute 
##' the dendrogram. See \code{help(heatmap.2)}
##' @param cluster the clustering algorithm used to build the dendrogram.
##' See \code{help(heatmap.2)}
##' @param dendro A boolean to indicate fi the dendrogram has to be displayed
##' @return A heatmap
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- mvFilter(Exp1_R25_pept[1:1000], "wholeMatrix", 6)
##' qData <- Biobase::exprs(obj)
##' heatmapD(qData)
heatmapD <- function(qData, distance="euclidean", cluster="complete", 
                    dendro = FALSE){
##Check parameters
# paramdist <- c("euclidean", "manhattan") 
# if (!(distance %in% paramdist)){
#     stop("Param distance is not correct.")
#     return (NULL)
# }
# 
# paramcluster <- c("ward.D", "average")
# if (!(cluster %in%  paramcluster)){
#     stop("Param clustering is not correct.")
#     return (NULL)
# }


# if (isTRUE(dendro) && getNumberOfEmptyLines(qData) != 0)  {
#     stop("Your dataset contains empty lines: the dendrogram cannot 
# be computed.
#         Please filter or impute missing values before.")
#     return (NULL)
# }
# else {
    .data <- matrix(qData, 
                    ncol = ncol(qData), 
                    byrow = FALSE,
                    dimnames = list(rownames(qData), colnames(qData))
    )
    colors = c(seq(-3, -2, length=100),
                seq(-2, 0.5, length=100),
                seq(0.5, 6, length=100))
    heatmap.color <- colorRampPalette(c("green", "red"))(n = 1000)
    
    
    if (dendro){ .dendro = "row"} else {.dendro = "none"}
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
        dendrogram =.dendro,
        Rowv=TRUE,
        col=heatmap.color ,
        density.info='none',
        key=TRUE,
        trace="none",
        scale="none",
        #srtCol=45,
        labCol="",
        margins=c(4,12),
        cexRow=1.5
    )
#    }
}


##' Heatmap inspired by the heatmap.2 function.
##' 
##' @title This function is inspired from the function \code{\link{heatmap.2}} 
##' that displays quantitative data in the \code{exprs()} table of an object of
##' class \code{MSnSet}. For more information, please refer to the help 
##' of the heatmap.2 function.
##' @param x A dataframe that contains quantitative data.
##' @param col colors used for the image. Defaults to heat colors (heat.colors).
##' @param srtCol angle of column labels, in degrees from horizontal 
##' @param labCol character vectors with column labels to use.
##' @param labRow character vectors with row labels to use.
##' @param key logical indicating whether a color-key should be shown.
##' @param key.title main title of the color key. If set to NA no title will 
##' be plotted.
##' @param main main title; default to none.
##' @param ylab y-axis title; default to none.
##' @return A heatmap
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- mvFilter(Exp1_R25_pept, "wholeMatrix", 6)
##' qData <- Biobase::exprs(obj)
##' heatmap.DAPAR(qData)
heatmap.DAPAR <- 
    function (x, 
              col = heat.colors(100),
              srtCol=NULL,
              labCol = NULL,
              labRow = NULL,
              key = TRUE, 
              key.title = NULL,
              main = NULL,  
              ylab = NULL) 
    {
        scale01 <- function(x, low = min(x), high = max(x)) {
            x <- (x - low)/(high - low)
            x
        }
        
        offsetCol <- 0.5
        offsetRow = 0.5
        srtRow = NULL
        colRow = NULL
        colCol = NULL 
        xlab = NULL
        key.par = list()
        margins = c(5, 5)
        sepcolor = "white"
        na.color = "white"
        keysize = 1.5
        breaks <- NULL
        na.rm = TRUE
        
        if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
            stop("`x' must be a numeric matrix")
        nr <- di[1]
        nc <- di[2]
        if (nr <= 1 || nc <= 1) 
            stop("`x' must have at least 2 rows and 2 columns")
        x <- x[nr:1,]
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
        cexCol = 0.2 + 1/log10(nc)
        cexRow = 0.2 + 1/log10(nr)
        iy <- 1:nr
        breaks <- length(col) + 1
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                      length = breaks)
        
        nbr <- length(breaks)
        ncol <- length(breaks) - 1
        
        min.breaks <- min(breaks)
        max.breaks <- max(breaks)
        x[x < min.breaks] <- min.breaks
        x[x > max.breaks] <- max.breaks
        lhei <- c(keysize, 4)
        lwid <- c(keysize, 4)
        lmat <- rbind(4:3, 2:1)
        lmat[is.na(lmat)] <- 0
        
        op <- par(no.readonly = TRUE)
        on.exit(par(op))
        layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
        
        par(mar = c(margins[1], 0, 0, margins[2]))
        x <- t(x)
        
        
        image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
                  c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
              breaks = breaks)
        
        
        if (!is.null(labCol)) 
        {
            axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 + 
                     offsetCol, tick = 0, cex.axis = cexCol, hadj = NA, 
                 padj = 0)
        }
        else {
                adjCol = c(1, NA)
                xpd.orig <- par("xpd")
                par(xpd = NA)
                xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, 
                             tick = 0)
                text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * 
                         strheight("M"), labels = labCol, adj = adjCol, 
                     cex = cexCol, srt = srtCol, col = colCol)
                par(xpd = xpd.orig)
        }
        
        
        if (!is.null(labRow) ) {
            axis(4, iy, labels = labRow, las = 5, line = -0.5 + offsetRow, 
                 tick = 0, cex.axis = cexRow, hadj = 0, padj = NA)
        }
        else {
                xpd.orig <- par("xpd")
                par(xpd = NA)
                ypos <- axis(4, iy, labels = rep("", nr), las = 2, 
                             line = -0.5, tick = 0)
                text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"), 
                     y = ypos, labels = labRow, adj = c(0,NA), cex = cexRow, 
                     srt = srtRow, col = colRow)
                par(xpd = xpd.orig)
        }
        
        
        
        
        par(mar = c(margins[1], 0, 0, 0))
        plot.new()
        par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
        
        plot.new()
        if (!is.null(main)) 
            title(main, cex.main = 1.5 * op[["cex.main"]])
        
        
        if (key) {
            mar <- c(5, 4, 2, 1)
            par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
            if (length(key.par) > 0) 
                do.call(par, key.par)
            
            tmpbreaks <- breaks
            min.raw <- min.breaks
            max.raw <- max.breaks
            
            z <- seq(min.raw, max.raw, by = min(diff(breaks)/100))
            image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
                  xaxt = "n", yaxt = "n")
            par(usr = c(0, 1, 0, 1))
            lv <- pretty(breaks)
            xv <- scale01(as.numeric(lv), min.raw, max.raw)
            xargs <- list(at = xv, labels = lv)
            
            xargs$side <- 1
            do.call(axis, xargs)
            key.xlab <- "Intensity value"
            
            mtext(side = 1, key.xlab, line = par("mgp")[1], padj = 0.5, 
                  cex = par("cex") * par("cex.lab"))
            
            if (is.null(key.title)) 
                title("Color Key")
        }
        
    }
