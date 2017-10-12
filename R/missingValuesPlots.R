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
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.mvPerLinesHisto(Exp1_R25_pept)
wrapper.mvPerLinesHisto <- function(obj, indLegend="auto", showValues=FALSE){
qData <- Biobase::exprs(obj)
samplesData <- Biobase::pData(obj)
mvPerLinesHisto(qData, samplesData, indLegend, showValues)
}


##' This method is a wrapper to plots from a \code{\link{MSnSet}} object a 
##' histogram which represents the distribution of the 
##' number of missing values (NA) per lines (ie proteins).
##' 
##' @title Histogram of missing values per lines from an object using highcharter
##' \code{\link{MSnSet}}
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param indLegend The indice of the column name's in \code{pData()} tab .
##' @param showValues A logical that indicates wether numeric values should be
##' drawn above the bars.
##' @return A histogram
##' @author Alexia Dorffer
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.mvPerLinesHisto(Exp1_R25_pept)
wrapper.mvPerLinesHisto_HC <- function(obj, indLegend="auto", showValues=FALSE){
    qData <- Biobase::exprs(obj)
    samplesData <- Biobase::pData(obj)
    hc <- mvPerLinesHisto_HC(qData, samplesData, indLegend, showValues)
    return(hc)
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
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' samplesData <- Biobase::pData(Exp1_R25_pept)
##' mvPerLinesHisto(qData, samplesData)
mvPerLinesHisto <- function(qData, samplesData, indLegend="auto", showValues=FALSE){

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

print(nb.na2barplot[-1])

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


##' This method plots a bar plot which represents the distribution of the 
##' number of missing values (NA) per lines (ie proteins).
##' 
##' @title Bar plot of missing values per lines using highcharter
##' @param qData A dataframe that contains the data to plot.
##' @param samplesData A dataframe which contains informations about 
##' the replicates.
##' @param indLegend The indice of the column name's in \code{pData()} tab 
##' @param showValues A logical that indicates wether numeric values should be
##' drawn above the bars.
##' @return A bar plot
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' samplesData <- Biobase::pData(Exp1_R25_pept)
##' mvPerLinesHisto_HC(qData, samplesData)
mvPerLinesHisto_HC <- function(qData, samplesData, indLegend="auto", showValues=FALSE){
    
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
    
    df <- data.frame(y=nb.na2barplot[-1])
    myColors = rep("lightgrey",nrow(df))
    myColors[nrow(df)] <- "red"
    
    #df1 <- df2 <- df
    #df2[1:(nrow(df)-1),] <- 0
    #df1 [nrow(df),] <- 0
    
    
    #, series = list( pointWidth = 50)
    
    h1 <-  highchart() %>% 
        hc_title(text = "#[lines] with X NA values") %>% 
        hc_add_series(data = df, type="column", colorByPoint = TRUE) %>%
        hc_colors(myColors) %>%
        hc_plotOptions( column = list(stacking = "normal"),
                        animation=list(duration = 100)) %>%
        hc_legend(enabled = FALSE) %>%
        hc_xAxis(categories = row.names(df), title = list(text = "#[NA values] per line")) %>%
      my_hc_ExportMenu(filename = "missingValuesPlot1") %>%
        hc_tooltip(enabled = FALSE)
    
    return(h1)
 
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
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.mvPerLinesHistoPerCondition(Exp1_R25_pept)
wrapper.mvPerLinesHistoPerCondition <- function(obj, indLegend="auto", 
                                            showValues=FALSE){
qData <- Biobase::exprs(obj)
samplesData <- Biobase::pData(obj)
mvPerLinesHistoPerCondition(qData, samplesData, indLegend, showValues)
}


##' This method is a wrapper to plots (using highcharts) from a \code{\link{MSnSet}} object a 
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
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.mvPerLinesHistoPerCondition(Exp1_R25_pept)
wrapper.mvPerLinesHistoPerCondition_HC <- function(obj, indLegend="auto", 
                                                showValues=FALSE){
    qData <- Biobase::exprs(obj)
    samplesData <- Biobase::pData(obj)
    mvPerLinesHistoPerCondition_HC(qData, samplesData, indLegend, showValues)
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
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' samplesData <- Biobase::pData(Exp1_R25_pept)
##' mvPerLinesHistoPerCondition(qData, samplesData)
mvPerLinesHistoPerCondition <- function(qData, samplesData, indLegend="auto", 
                                        showValues=FALSE){

pal <- getPaletteForLabels(samplesData[,"Label"])
if (identical(indLegend,"auto")) { indLegend <- c(2:length(colnames(samplesData)))}

nbLabels <- length(unique(samplesData[,"Label"]))

ncolMatrix <- max(unlist(lapply(unique(samplesData[,"Label"]), function(x){length(which(samplesData[,"Label"]==x))})))
m <- matrix(rep(0, nbLabels*(1+ncolMatrix)), 
            ncol = nbLabels, 
            dimnames=list(seq(0:(ncolMatrix)),unique(samplesData[,"Label"])))

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
                names.arg = as.character(0:ncolMatrix), 
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




##' This method plots a bar plot which represents the distribution of the 
##' number of missing values (NA) per lines (ie proteins) and per conditions.
##' Same as the function \link{mvPerLinesHistoPerCondition} but uses the package
##' \CRANpkg{highcharter}.
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
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' samplesData <- Biobase::pData(Exp1_R25_pept)
##' mvPerLinesHistoPerCondition_HC(qData, samplesData)
mvPerLinesHistoPerCondition_HC <- function(qData, samplesData, indLegend="auto", 
                                        showValues=FALSE){
    
    labels <- samplesData[,"Label"]
        pal <- getPaletteForLabels(samplesData[,"Label"])
    if (identical(indLegend,"auto")) { indLegend <- c(2:length(colnames(samplesData)))}
    
    nbLabels <- length(unique(samplesData[,"Label"]))
    
    ncolMatrix <- max(unlist(lapply(unique(samplesData[,"Label"]), function(x){length(which(samplesData[,"Label"]==x))})))
    m <- matrix(rep(0, nbLabels*(1+ncolMatrix)), 
                ncol = nbLabels, 
                dimnames=list(seq(0:(ncolMatrix)),unique(samplesData[,"Label"])))
    
    for (i in unique(samplesData[,"Label"]))
    {
        nSample <- length(which(samplesData[,"Label"] == i))
        t <- NULL
        if (nSample == 1) {
            t <- table(as.integer(is.na(qData[,which(samplesData[,"Label"] == i)])))
        } else {t <- table(rowSums(is.na(qData[,which(samplesData[,"Label"] == i)])))}
        
        m[as.integer(names(t))+1,i] <- t
    }
    
    
    colnames(m) <- unique(labels)
    rownames(m) <- 0:(nrow(m)-1)
    
    nbSeries = length(unique(labels))
    series <- list()
    for (i in 1:nbSeries){
        tmp <- list(data = m[,i], name=(unique(labels))[i])
        names(tmp$data) <- c()
        series[[i]] <- tmp
    }
    
    h1 <-  highchart() %>% 
        hc_title(text = "#[lines] with X NA values (condition-wise)") %>% 
        my_hc_chart(chartType = "column") %>%
        hc_plotOptions( column = list(stacking = ""),
                        animation=list(duration = 100)) %>%
        hc_add_series_list(series) %>%
        hc_legend(enabled = FALSE) %>%
        hc_xAxis(categories = row.names(m), title = list(text = "#[NA values] per line (condition-wise)")) %>%
      my_hc_ExportMenu(filename = "missingValuesPlot_2") %>%
        hc_tooltip(headerFormat= '',
                   pointFormat = "{point.y} ")
    
    return(h1)
    

    
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
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.mvHisto(Exp1_R25_pept, showValues=TRUE)
wrapper.mvHisto <- function(obj, indLegend="auto", showValues=FALSE){
qData <- Biobase::exprs(obj)
samplesData <- Biobase::pData(obj)
labels <- samplesData[,"Label"]
mvHisto(qData, samplesData, labels, indLegend, showValues)
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
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.mvHisto(Exp1_R25_pept, showValues=TRUE)
wrapper.mvHisto_HC <- function(obj, indLegend="auto", showValues=FALSE){
    qData <- Biobase::exprs(obj)
    samplesData <- Biobase::pData(obj)
    labels <- samplesData[,"Label"]
    mvHisto_HC(qData, samplesData, labels, indLegend, showValues)
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
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' samplesData <- Biobase::pData(Exp1_R25_pept)
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
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



##' This method plots a histogram of missing values. Same as the function \code{mvHisto}
##' but uses the package \CRANpkg{highcharter}
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
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' samplesData <- Biobase::pData(Exp1_R25_pept)
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
##' mvHisto_HC(qData, samplesData, labels, indLegend="auto", showValues=TRUE)
mvHisto_HC <- function(qData, samplesData, labels, indLegend="auto", 
                    showValues=FALSE){
    
    if (identical(indLegend,"auto")) { 
        indLegend <- c(2:length(colnames(samplesData)))
    }
    
    
    colnames(qData) <- samplesData[,"Label"]
    
    coeffMax <- .1
    pal <- getPaletteForLabels(labels)
    
    NbNAPerCol <- colSums(is.na(qData))
    NbNAPerRow <- rowSums(is.na(qData))
    
    # if (sum(NbNAPerCol) == 0){
    #     NbNAPerCol <- rep(0,1+ncol(qData))
    # }
    
    
    nbSeries = length(unique(labels))
    series <- list()
    for (i in 1:nbSeries){
        data <- rep(0, length(NbNAPerCol))
        data[which((unique(labels))[i] == labels)] <- NbNAPerCol[which(labels == (unique(labels))[i])]
        
        tmp <- list(data = data, name=(unique(labels))[i])
        names(tmp$data) <- c()
        series[[i]] <- tmp
    }
    
    h1 <-  highchart() %>%
         my_hc_chart(chartType = "column") %>%
         hc_title(text = "#[non-NA values] by replicate") %>%
        hc_add_series_list(series) %>%
        hc_plotOptions( column = list(stacking = "normal"),
                        animation=list(duration = 100)) %>%
        hc_legend(enabled = FALSE) %>%
        hc_xAxis(categories = labels, title = list(text = "Replicates")) %>%
      my_hc_ExportMenu(filename = "missingValuesPlot_3") %>%
        hc_tooltip(headerFormat= '',
                   pointFormat = "{point.y}")
    

    return(h1)
    
    
    

    
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
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.mvImage(Exp1_R25_pept)
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
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
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
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.mvTypePlot(Exp1_R25_pept)
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
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
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




