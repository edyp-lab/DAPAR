##' Returns the contents of the slot processing of an object of class \code{\link{MSnSet}}
##' 
##' @title Returns the contains of the slot processing  of an object of 
##' class \code{\link{MSnSet}}
##' @param  obj An object (peptides) of class \code{\link{MSnSet}}.
##' @return The slot processing of obj@processingData
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' getProcessingInfo(Exp1_R25_pept)
getProcessingInfo <- function(obj){
return(obj@processingData@processing)
}

##' Returns the number of empty lines in a matrix.
##' 
##' @title Returns the number of empty lines in the data
##' @param qData A matrix corresponding to the quantitative data.
##' @return An integer
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' getNumberOfEmptyLines(qData)
getNumberOfEmptyLines <- function(qData){
n <- sum(apply(is.na(as.matrix(qData)), 1, all))
return (n)
}




##' Returns a list for the two conditions where each slot is a vector of 
##' indices for the samples.
##' 
##' @title Gets the conditions indices.
##' @param labels A vector of strings containing the column "Label" of 
##' the \code{pData()}.
##' @param cond1 A vector of Labels (a slot in the \code{pData()} table) for
##' the condition 1.
##' @param cond2 A vector of Labels (a slot in the \code{pData()} table) for
##' the condition 2.
##' @return A list with two slots \code{iCond1} and \code{iCond2} containing
##' respectively the indices of samples in the \code{pData()} table of the
##' dataset. 
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
##' getIndicesConditions(labels, "25fmol", "10fmol")
getIndicesConditions <- function(labels, cond1, cond2){
indCondition1 <- indCondition2 <- NULL

for(i in 1:length(cond1)){
    indCondition1 <- c(indCondition1,
                        which(labels == cond1[i]))
}
for(i in 1:length(cond2)){
    indCondition2 <- c(indCondition2,
                        which(labels == cond2[i]))
}

return(list(iCond1 = indCondition1, iCond2 = indCondition2))
}



##' Selects colors for the plots in \code{\link{DAPAR}} based on the different conditions in
##' the dataset. The palette is derived from
##' the brewer palette "Dark2" (see \code{\link{RColorBrewer}}).
##' 
##' @title Palette for plots in \code{\link{DAPAR}}
##' @param labels A vector of labels (strings).
##' @return A palette designed for the data manipulated in DAPAR
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
##' getPaletteForLabels(labels)
getPaletteForLabels <- function(labels){
    nColors <- 8
    col <- c(1:nColors)
    palette(brewer.pal(nColors,"Dark2"))
    
    ## Define one color per label/condition
    col.boxplot <- NULL
    for (i in 1:length(labels)){
        col.boxplot[which(labels == unique(labels)[i])] <- col[i]
    }
    
    return (col.boxplot)
}


##' Selects colors for the highcharter plots in \code{\link{DAPAR}} based on the different conditions in
##' the dataset. The palette is derived from
##' the brewer palette "Dark2" (see \code{\link{RColorBrewer}}).
##' 
##' @title Palette for highcharter plots in \code{\link{DAPAR}}
##' @param labels A vector of labels (strings).
##' @return A palette designed for the data manipulated in \code{\link{DAPAR}}
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' labels <- Biobase::pData(Exp1_R25_pept)[,"Label"]
##' getPaletteForLabels_HC(labels)
getPaletteForLabels_HC <- function(labels){
  nColors <- 8
  #col <- c(1:nColors)
  #pal <- c('#7cb5ec', '#434348', '#90ed7d', '#f7a35c', '#8085e9', '#f15c80', '#e4d354', '#2b908f', '#f45b5b', '#91e8e1')
  #pal <- c("#002F80", "#002F80","#002F80","#002F80","#F9AF38","#F9AF38","#F9AF38","#F9AF38")
  pal <- brewer.pal(nColors, "Dark2")
  ## Define one color per label/condition
  col.boxplot <- NULL
  for (i in 1:length(labels)){
    col.boxplot[which(labels == unique(labels)[i])] <- pal[i]
  }
  
  return (col.boxplot)
}

##' Selects colors for the highcharter plots in \code{\link{DAPAR}} based on the replicates in
##' the dataset. The palette is derived from the brewer palette "Dark2" (see \code{\link{RColorBrewer}}).
##' 
##' @title Palette for highcharter plot the replicates in \code{\link{DAPAR}}
##' @param nColors The desired number of colors
##' @return A palette designed for the data manipulated in \code{\link{DAPAR}}
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' n <- nrow(Biobase::pData(Exp1_R25_pept))
##' getPaletteForReplicates_HC(n)
getPaletteForReplicates_HC <- function(nColors){
  #col <- c(1:nColors)
  #pal <- c("#002F80", "#002F80","#002F80","#002F80","#F9AF38","#F9AF38","#F9AF38","#F9AF38")
  pal <- brewer.pal(nColors, "Dark2")
  return(pal)
}



##' Selects colors for the plots in DAPAR based on the replicates in
##' the dataset. The palette is derived from
##' the brewer palette "Dark2" (see \code{\link{RColorBrewer}}).
##' 
##' @title Palette for plot the replicates in DAPAR
##' @param nColors The desired number of colors
##' @return A palette designed for the data manipulated in DAPAR
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' n <- nrow(Biobase::pData(Exp1_R25_pept))
##' getPaletteForReplicates(n)
getPaletteForReplicates <- function(nColors){
#col <- c(1:nColors)
#getPalette <- colorRampPalette(brewer.pal(8, "Dark2"))
#return(col)
pal <- brewer.pal(nColors, "Dark2")
return(pal)
}



##' Customise the contextual menu of highcharts plots.
##' 
##' @title Customised contextual menu of highcharts plots
##' @param hc A highcharter object
##' @param filename The filename under which the plot has to be saved
##' @return A contextual menu for highcharts plots
##' @author Samuel Wieczorek
##' @examples
##' library("highcharter")
##' hc <- highchart() 
##' hc_chart(hc,type = "line") 
##' hc_add_series(hc,data = c(29, 71, 40))
##' my_hc_ExportMenu(hc,filename='foo')
my_hc_ExportMenu <- function(hc, filename){
  hc_exporting(hc, enabled=T,
               filename = filename,
               buttons= list(
                 contextButton= list(
                   menuItems= list('downloadPNG', 'downloadSVG','downloadPDF')
                 )
               )
  )
}



##' Customise the resetZoomButton of highcharts plots.
##' 
##' @title Customised resetZoomButton of highcharts plots
##' @param hc A highcharter object
##' @param chartType The type of the plot
##' @param zoomType The type of the zoom (one of "x", "y", "xy", "None")
##' @return A highchart plot
##' @author Samuel Wieczorek
##' @examples
##' library("highcharter")
##' hc <- highchart() 
##' hc_chart(hc,type = "line") 
##' hc_add_series(hc,data = c(29, 71, 40))
##' my_hc_ExportMenu(hc,filename='foo')
my_hc_chart <- function(hc,  chartType,zoomType="None"){
  hc %>% 
    hc_chart(type = chartType, 
           zoomType=zoomType,
           resetZoomButton= list(
             position = list(
               align= 'left',
               verticalAlign = 'top')
           ))
}
