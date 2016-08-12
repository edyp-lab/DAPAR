##' This method is a wrapper to plots from a \code{\link{MSnSet}} object a 
##' histogram which represents the distribution of the 
##' number of missing values (NA) per lines (ie proteins).
##' 
##' @title Histogram of missing values per lines from an object 
##' \code{\link{MSnSet}}
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param indLegend The indice of the column name's in \code{pData()} tab .
##' @param showValues A logical that indicates wether numeric values should be
##' drawn above the bars.
##' @return A histogram
##' @author Alexia Dorffer
##' @examples
##' data(UPSpep25)
##' wrapper.mvPerLinesHisto(UPSpep25)
wrapper.mvPerLinesHisto <- function(obj, indLegend="auto", showValues=FALSE){
qData <- Biobase::exprs(obj)
samplesData <- Biobase::pData(obj)
mvPerLinesHisto(qData, samplesData, indLegend, showValues)
}

##' This method plots a bar plot which represents the distribution of the 
##' number of missing values (NA) per lines (ie proteins).
##' 
##' @title Bar plot of missing values per lines
##' @param qData A dataframe that contains the data to plot.
##' @param samplesData A dataframe which contains informations about 
##' the replicates.
##' @param indLegend The indice of the column name's in \code{pData()} tab 
##' @param showValues A logical that indicates wether numeric values should be
##' drawn above the bars.
##' @return A bar plot
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' qData <- Biobase::exprs(UPSpep25)
##' samplesData <- Biobase::pData(UPSpep25)
##' mvPerLinesHisto(qData, samplesData)
mvPerLinesHisto <- function(qData, samplesData, indLegend="auto", 
                        showValues=FALSE){

if (identical(indLegend,"auto")) { indLegend <- c(2:length(colnames(samplesData)))}
    
for (j in 1:length(colnames(qData))){
    noms <- NULL
    for (i in 1:length(indLegend)){
    noms <- paste(noms, samplesData[j,indLegend[i]], sep=" ")
    }
    colnames(qData)[j] <- noms
}

coeffMax <- .1

NbNAPerCol <- colSums(is.na(qData))
NbNAPerRow <- rowSums(is.na(qData))
#par(mar = c(10,3, 3, 3))

nb.col <- dim(qData)[2] 
nb.na <- NbNAPerRow
temp <- table(NbNAPerRow)
nb.na2barplot <- c(temp, rep(0,1+ncol(qData)-length(temp)))

if (sum(NbNAPerRow) == 0){
    nb.na2barplot <- rep(0,1+ncol(qData))
}
x <- barplot(nb.na2barplot[-1], 
                main = "# lines by # of NA",
                xlab = "# NA per lines",
                names.arg = as.character(c(1:(ncol(qData)))), 
                col = c(rep("lightgrey",nb.col-1), "red"),
                ylim = c(0, 1.2*max(1,nb.na2barplot[-1])), 
                las=1,
                cex.names=1.5,
                cex.axis=1.5
)

# if (showValues) {
#   graphics::text(x, 
#        nb.na2barplot[-1], 
#        labels = nb.na2barplot[-1],
#        cex= .75,
#        pos = 3,
#        col = rep("black",nb.col)) 
#   
# }
}



##' This method is a wrapper to plots from a \code{\link{MSnSet}} object a 
##' bar plot which represents the distribution of the 
##' number of missing values (NA) per lines (ie proteins) and per conditions.
##' 
##' @title Bar plot of missing values per lines and per conditions from an 
##' object \code{\link{MSnSet}}
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param indLegend The indice of the column name's in \code{pData()} tab .
##' @param showValues A logical that indicates wether numeric values should be
##' drawn above the bars.
##' @return A bar plot
##' @author Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' wrapper.mvPerLinesHistoPerCondition(UPSpep25)
wrapper.mvPerLinesHistoPerCondition <- function(obj, indLegend="auto", 
                                            showValues=FALSE){
qData <- Biobase::exprs(obj)
samplesData <- Biobase::pData(obj)
mvPerLinesHistoPerCondition(qData, samplesData, indLegend, showValues)
}

##' This method plots a bar plot which represents the distribution of the 
##' number of missing values (NA) per lines (ie proteins) and per conditions.
##' 
##' @title Bar plot of missing values per lines and per condition
##' @param qData A dataframe that contains quantitative data.
##' @param samplesData A dataframe where lines correspond to samples and 
##' columns to the meta-data for those samples.
##' @param indLegend The indice of the column name's in \code{pData()} tab 
##' @param showValues A logical that indicates wether numeric values should be
##' drawn above the bars.
##' @return A bar plot
##' @author Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' qData <- Biobase::exprs(UPSpep25)
##' samplesData <- Biobase::pData(UPSpep25)
##' mvPerLinesHistoPerCondition(qData, samplesData)
mvPerLinesHistoPerCondition <- function(qData, samplesData, indLegend="auto", 
                                        showValues=FALSE){

pal <- getPaletteForLabels(samplesData[,"Label"])
if (identical(indLegend,"auto")) { indLegend <- c(2:length(colnames(samplesData)))}

nbLabels <- length(unique(samplesData[,"Label"]))

m <- matrix(rep(0,nbLabels*(1+(nrow(samplesData)/nbLabels ))), 
            ncol = nbLabels, 
            dimnames=list(seq(0:(nrow(samplesData)/nbLabels)),
                            unique(samplesData[,"Label"])))

for (i in unique(samplesData[,"Label"]))
{
    nSample <- length(which(samplesData[,"Label"] == i))
    t <- NULL
    if (nSample == 1) {
        t <- table(as.integer(is.na(qData[,which(samplesData[,"Label"] == i)])))
    } else {t <- table(rowSums(is.na(qData[,which(samplesData[,"Label"] == i)])))}
    
    m[as.integer(names(t))+1,i] <- t
}

m <- t(m)

x <- barplot(m, 
                main = "# lines by # of NA",
                xlab = "# NA per lines",
                names.arg = as.character(
                    seq(0:(nrow(samplesData)/nbLabels))-1), 
                col = unique(pal),
                ylim = c(0, 1.2*max(m)), 
                xpd = FALSE,
                las=1,
                cex.names=1.5,
                cex.axis=1.5,
                beside=TRUE
)

# if (showValues) {
#   text(x, 
#        m, 
#        labels = m,
#        cex= .75,
#        pos = 3) 
# }

}


##' This method plots from a \code{\link{MSnSet}} object a histogram of 
##' missing values.
##' 
##' @title Histogram of missing values from a \code{\link{MSnSet}} object
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param indLegend The indices of the column name's in \code{pData()} tab.
##' @param showValues A logical that indicates wether numeric values should be
##' drawn above the bars.
##' @return A histogram
##' @author Alexia Dorffer
##' @examples
##' data(UPSpep25)
##' wrapper.mvHisto(UPSpep25, showValues=TRUE)
wrapper.mvHisto <- function(obj, indLegend="auto", showValues=FALSE){
qData <- Biobase::exprs(obj)
samplesData <- Biobase::pData(obj)
labels <- samplesData[,"Label"]
mvHisto(qData, samplesData, labels, indLegend, showValues)
}



##' This method plots a histogram of missing values.
##' 
##' @title Histogram of missing values
##' @param qData A dataframe that contains quantitative data.
##' @param samplesData A dataframe where lines correspond to samples and 
##' columns to the meta-data for those samples.
##' @param labels A vector of the conditions (labels) (one label per sample).
##' @param indLegend The indices of the column name's in \code{pData()} tab
##' @param showValues A logical that indicates wether numeric values should be
##' drawn above the bars.
##' @return A histogram
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' qData <- Biobase::exprs(UPSpep25)
##' samplesData <- Biobase::pData(UPSpep25)
##' labels <- Biobase::pData(UPSpep25)[,"Label"]
##' mvHisto(qData, samplesData, labels, indLegend="auto", showValues=TRUE)
mvHisto <- function(qData, samplesData, labels, indLegend="auto", 
                    showValues=FALSE){

if (identical(indLegend,"auto")) { 
    indLegend <- c(2:length(colnames(samplesData)))
}


colnames(qData) <- samplesData[,"Label"]

coeffMax <- .1
pal <- getPaletteForLabels(labels)

NbNAPerCol <- colSums(is.na(qData))
NbNAPerRow <- rowSums(is.na(qData))

if (sum(NbNAPerCol) == 0) {if (sum(NbNAPerCol) == 0){
    NbNAPerCol <- rep(0,1+ncol(qData))
}} 
x <- barplot(NbNAPerCol, 
                main = "# NA per columns",
                col=pal,
                las=1,
                ylim = c(0, 1.2*max(1,NbNAPerCol)),
                names.arg = c(1:18), 
                cex.names=1.5,
                cex.axis=1.5,
                axisnames = FALSE
)

par(xpd = TRUE)
graphics::text(x, -3,
        labels = colnames(qData),
        srt = 45,
        adj=1,
        cex=1.4)
# 
#   
#   if (showValues) {
#     graphics::text(x, 
#          NbNAPerCol, 
#          labels = NbNAPerCol, 
#          cex= .75,
#          pos = 3) 
#   }

}


##' Plots a heatmap of the quantitative data. Each column represent one of
##' the conditions in the object of class \code{\link{MSnSet}} and 
##' the color is proportional to the mean of intensity for each line of
##' the dataset.
##' The lines have been sorted in order to vizualize easily the different
##' number of missing values. A white square is plotted for missing values.
##' 
##' @title Heatmap of missing values from a \code{\link{MSnSet}} object
##' @param obj An object of class \code{\link{MSnSet}}.
##' @return A heatmap
##' @author Alexia Dorffer
##' @examples
##' data(UPSpep25)
##' wrapper.mvImage(UPSpep25)
wrapper.mvImage <- function(obj){
qData <- Biobase::exprs(obj)
labels <- Biobase::pData(obj)[,"Label"]
mvImage(qData, labels)
}



##' Plots a heatmap of the quantitative data. Each column represent one of
##' the conditions in the object of class \code{\link{MSnSet}} and 
##' the color is proportional to the mean of intensity for each line of
##' the dataset.
##' The lines have been sorted in order to vizualize easily the different
##' number of missing values. A white square is plotted for missing values.
##' 
##' @title Heatmap of missing values
##' @param qData A dataframe that contains quantitative data.
##' @param labels A vector of the conditions (labels) (one label per sample).
##' @return A heatmap
##' @author Samuel Wieczorek, Thomas Burger
##' @examples
##' data(UPSpep25)
##' qData <- Biobase::exprs(UPSpep25)
##' labels <- Biobase::pData(UPSpep25)[,"Label"]
##' mvImage(qData, labels)
mvImage <- function(qData, labels){
### build indices of conditions
indCond <- list()
LabelNames <- unique(labels)
for (i in LabelNames) {
    indCond <- append(indCond, list(which(i == labels)))
}
indCond <- setNames(indCond, as.list(c("cond1", "cond2")))

nNA1 = apply(as.matrix(qData[,indCond$cond1]), 1, function(x) sum(is.na(x)))
nNA2 = apply(as.matrix(qData[,indCond$cond2]), 1, function(x) sum(is.na(x)))
o <- order(((nNA1 +1)^2) / (nNA2 +1))
exprso <- qData[o,]

for (i in 1:nrow(exprso)){
    k <- order(exprso[i,indCond$cond1])
    exprso[i,rev(indCond$cond1)] <- exprso[i, k]
    .temp <- mean(exprso[i,rev(indCond$cond1)], na.rm = TRUE)
    exprso[i,which(!is.na(exprso[i,indCond$cond1]))] <- .temp
    
    k <- order(exprso[i,indCond$cond2])
    exprso[i,indCond$cond2] <- exprso[i, k+length(indCond$cond1)]
    .temp <- mean(exprso[i,indCond$cond2], na.rm = TRUE)
    exprso[i,length(indCond$cond1) + 
            which(!is.na(exprso[i,indCond$cond2]))] <- .temp
}

colfunc <- colorRampPalette(c("yellow", "red"))

heatmap.DAPAR(exprso,
                col = colfunc(100),
                key=TRUE,
                srtCol= 0,
                labCol=labels,
                ylab = "Peptides / proteins",
                main = "Missing values heatmap"
) 
}




##' This method plots a scatter plot which represents the distribution of
##' missing values.
##' The colors correspond to the different conditions (slot Label in in the
##' dataset of class \code{\link{MSnSet}}).
##' The x-axis represent the mean of intensity for one condition and one
##' entity in the dataset (i. e. a protein) 
##' whereas the y-axis count the number of missing values for this entity
##' and the considered condition.
##' The data have been jittered for an easier vizualisation.
##' 
##' @title Distribution of missing values with respect to intensity values 
##' from a \code{\link{MSnSet}} object
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param threshold An integer for the intensity that delimits MNAR and 
##' MCAR missing values.
##' @return A scatter plot
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' wrapper.mvTypePlot(UPSpep25)
wrapper.mvTypePlot <- function(obj, threshold=0){
qData <- Biobase::exprs(obj)
labels <- Biobase::pData(obj)[,"Label"]
mvTypePlot(qData, labels, threshold)
}




##' This method plots a scatter plot which represents the distribution of
##' missing values.
##' The colors correspond to the different conditions (slot Label in in the
##' dataset of class \code{\link{MSnSet}}).
##' The x-axis represent the mean of intensity for one condition and one
##' entity in the dataset (i. e. a protein) 
##' whereas the y-axis count the number of missing values for this entity
##' and the considered condition.
##' The data have been jittered for an easier vizualisation.
##' 
##' @title Distribution of missing values with respect to intensity values
##' @param qData A dataframe that contains quantitative data.
##' @param labels A vector of the conditions (labels) (one label per sample).
##' @param threshold An integer for the intensity that delimits MNAR and 
##' MCAR missing values.
##' @return A scatter plot
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' qData <- Biobase::exprs(UPSpep25)
##' labels <- Biobase::pData(UPSpep25)[,"Label"]
##' mvTypePlot(qData, labels, threshold=0)
mvTypePlot <- function(qData, labels, threshold=0){
#require(scales)
pal <- unique(getPaletteForLabels(labels))
color <- NULL
col.legend <- c(1:length(pal))


conditions <- labels
mTemp <- colorTemp <- nbNA <- matrix(rep(0,nrow(qData)*length(unique(conditions))), nrow=nrow(qData),
                dimnames=list(NULL,unique(conditions)))

for (iCond in unique(conditions)){
    if (length(which(conditions==iCond)) == 1){
        mTemp[,iCond] <- qData[,which(conditions==iCond)]
        nbNA[,iCond] <- as.integer(is.na(qData[,which(conditions==iCond)]))
    }else {
        mTemp[,iCond] <- apply(qData[,which(conditions==iCond)], 1, mean, na.rm=TRUE)
        nbNA[,iCond] <- apply(qData[,which(conditions==iCond)],1,function(x) length(which(is.na(x) == TRUE)))
    }
    colorTemp[,iCond] <- pal[which( unique(conditions) ==iCond)]
}


plot(c(mTemp),
        rep(-1,length(c(mTemp))),
        xlim = range(c(mTemp), na.rm = TRUE),
        ylim = c(0, ncol(qData)/length(unique(conditions))),
        xlab = "Mean of quantity values", 
        ylab = "Number of missing values",
        main =  "Missing values repartition")


    points(c(mTemp),
            jitter(c(nbNA), 0.3), 
            col = alpha(c(colorTemp), 0.5),
            pch = 16,
            cex=0.8)

  #  abline(v=threshold, col="blue", lwd=3)


legend("topright"         
        , legend = unique(labels)
        , col = col.legend
        , pch = 15 
        , bty = "n"
        , pt.cex = 2
        , cex = 1
        , horiz = FALSE
        , inset=c(0,0)
)
}
